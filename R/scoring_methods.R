# =========================================================
# scoring_methods.R
# Metodos de ordenamento de SNPs usados no artigo (main.tex, Secao 2):
#   1) Correlacao Marginal        -> score_MAR
#   2) Correlacao Parcial         -> score_PAR
#   3) CAR Score                  -> score_CAR
#   4) I Score                    -> score_I_Score
#   5) Importancia por Permutacao -> score_RF_Permutation (%IncMSE)
#   6) Importancia por Impureza   -> score_RF_Impureza (IncNodePurity)
# (score_ValorP_* tambem e calculado apenas para conferir, no texto do artigo
# (Secao 2.1), a equivalencia entre ordenamento por valor-p e por correlacao
# marginal em regressao linear simples; nao gera coluna propria nas tabelas.)
#
# Este arquivo NAO inclui Kemeny/BARD/BARDHM: o artigo trata somente dos
# metodos de ordenamento acima, sem agregacao bayesiana de ranks.
# =========================================================

# ---------------------------------------------------------
# Secao: Calculo do I-Score de um SNP
# Objetivo: Medir o quanto um SNP individual ajuda a explicar o fenotipo.
# Entrada: Vetor genotipico de um SNP e vetor resposta `y`.
# Saida: Valor numerico do I-Score.
# Para que isso existe: O I-Score e um dos metodos base usados para ordenar SNPs.
# ---------------------------------------------------------
# I-Score marginal para um único SNP.
compute_i_score_single <- function(snp_vector, y) {
  # n = tamanho da amostra.
  n <- length(y)
  # Média global do fenótipo.
  y_bar <- mean(y)
  # Variância amostral na parametrização usada pelo I-score.
  s2 <- sum((y - y_bar)^2) / n
  # Se não há variabilidade no fenótipo, I-score é zero.
  if (s2 <= 0) return(0)
  # Níveis genotípicos observados no SNP.
  genotype_levels <- sort(unique(snp_vector))
  # Acumulador do I-score.
  i_score_total <- 0
  for (g in genotype_levels) {
    idx <- which(snp_vector == g)
    n_g <- length(idx)
    y_g_bar <- mean(y[idx])
    # Formas equivalentes do I-Score:
    # I = sum_j (n_j/n) * ((Ybar_j - Ybar)^2 / (s^2 / n_j))
    #   = [sum_j n_j^2 * (Ybar_j - Ybar)^2] / [sum_i (Y_i - Ybar)^2]
    contribution_g <- (n_g / n) * ((y_g_bar - y_bar)^2 / (s2 / n_g))
    i_score_total <- i_score_total + contribution_g
  }
  return(i_score_total)
}

# ---------------------------------------------------------
# Secao: Calculo de valor-p por SNP
# Objetivo: Ajustar regressao simples SNP a SNP e extrair valor-p bruto e ajustado.
# Entrada: Data frame com colunas de SNPs e a ultima coluna como fenotipo.
# Saida: Tabela com valor-p bruto e valor-p ajustado por Bonferroni.
# Para que isso existe: Permite comparar os metodos de score com uma referencia estatistica classica.
# ---------------------------------------------------------
# Função para cálculo do valor-p bruto e ajustado (Bonferroni), conforme solicitado.
valor.p <- function(genotipo_fenotipo) {
  # Vetor de p-valores brutos (um por SNP).
  valor_p_bruto <- vector()
  # Vetor de p-valores ajustados (Bonferroni).
  valor_p_ajustado <- vector()
  # Número de testes (todas as colunas exceto fenótipo).
  m <- ncol(genotipo_fenotipo) - 1
  # Guarda modelos lineares ajustados (opcional para inspeção).
  model_regression <- list()

  for (i in seq_len(ncol(genotipo_fenotipo) - 1)) {
    x <- genotipo_fenotipo[, i]  # Genótipo do SNP i.
    # Ajusta modelo linear simples fenótipo ~ SNP_i.
    model_regression[[i]] <- lm(genotipo_fenotipo[, ncol(genotipo_fenotipo)] ~ x)
    # Se SNP é constante, força p-valor = 1 para evitar artefatos numéricos.
    valor_p_bruto[i] <- ifelse(
      all(x == 0) || all(x == 1) || all(x == 2) || all(x == 3),
      1,
      summary(model_regression[[i]])[[4]][2, 4]
    )
  }
  # Ajuste de Bonferroni.
  valor_p_ajustado <- m * valor_p_bruto
  # Cria tabela de saída.
  valor_p <- data.frame()
  valor_p <- cbind(valor_p_bruto, valor_p_ajustado)
  colnames(valor_p) <- c("Valor p bruto", "Valor p ajustado")
  rownames(valor_p) <- names(genotipo_fenotipo)[1:(ncol(genotipo_fenotipo) - 1)]
  return(valor_p)
}

# ---------------------------------------------------------
# Secao: Calculo de todos os scores base
# Objetivo: Executar todos os metodos de score do projeto em um unico ponto.
# Entrada: Matriz de genotipos `X` e fenotipo `y`.
# Saida: Tabela com uma coluna `score_*` para cada metodo.
# Para que isso existe: Centraliza o calculo dos metodos usados depois para gerar ranks.
# ---------------------------------------------------------
# Calcula scores CAR, MAR, PAR e importância RF. Retorna data.frame com uma coluna por método.
compute_all_scores <- function(X, y) {
  message("[scores] Iniciando cálculo de scores...")
  # Nomes dos SNPs vindos das colunas de X.
  snp_names <- colnames(X)
  # Se não houver nomes, cria nomes padrão SNP1..SNPp.
  if (is.null(snp_names)) snp_names <- paste0("SNP", seq_len(ncol(X)))
  # n = indivíduos; p = número de SNPs.
  n <- nrow(X)
  p <- ncol(X)
  # Converte para matriz numérica para funções estatísticas.
  Xmat <- as.matrix(X)
  # Garante vetor numérico da resposta.
  y <- as.numeric(y)

  # Data.frame base com coluna identificadora.
  out <- data.frame(SNP = snp_names, stringsAsFactors = FALSE)
  # Vetor nomeado para armazenar tempos (em segundos) de cada método de score.
  timing_seconds <- numeric(0)
  # Política automática de shrinkage:
  # se p >= 0.5*n, usa shrinkage automático para CAR/PAR (e MAR por consistência de implementação).
  use_auto_shrink <- (p >= 0.5 * n)
  # Tabela com lambda por método (quando aplicável).
  lambda_info <- data.frame(
    metodo = c("CAR", "MAR", "PAR"),
    lambda = c(NA_real_, NA_real_, NA_real_),
    modo = c("manual", "manual", "manual"),
    regra_aplicada = rep(ifelse(use_auto_shrink, "auto_if_p_ge_0.5n", "manual_lambda_0_if_p_lt_0.5n"), 3),
    n = rep(n, 3),
    p = rep(p, 3),
    limiar = rep(0.5 * n, 3),
    stringsAsFactors = FALSE
  )

  # Extrai lambda de objetos que eventualmente trazem essa informação em atributos/listas.
  extract_lambda <- function(obj) {
    if (is.null(obj)) return(NA_real_)
    if (is.list(obj) && !is.null(obj$lambda)) {
      v <- suppressWarnings(as.numeric(obj$lambda))
      if (length(v) > 0 && is.finite(v[1])) return(v[1])
    }
    at <- attributes(obj)
    if (!is.null(at) && "lambda" %in% names(at)) {
      v <- suppressWarnings(as.numeric(at$lambda))
      if (length(v) > 0 && is.finite(v[1])) return(v[1])
    }
    NA_real_
  }

  # Estimativa fallback de lambda de correlação para X (quando método não expõe lambda diretamente).
  estimate_lambda_x <- function(Xmat_local) {
    lam <- tryCatch(
      suppressWarnings(as.numeric(corpcor::estimate.lambda(Xmat_local))),
      error = function(e) NA_real_
    )
    if (length(lam) > 0 && is.finite(lam[1])) lam[1] else NA_real_
  }

  # CAR (Correlation-Adjusted marginal coRelation)
  message("[scores] Calculando CAR...")
  t0 <- Sys.time()
  tryCatch({
    # carscore com diagonal=FALSE usa ajuste de correlação entre SNPs.
    # Política: em alta dimensão usa shrinkage automático.
    if (use_auto_shrink) {
      lam_car <- estimate_lambda_x(Xmat)
      if (!is.finite(lam_car)) lam_car <- 0
      car <- care::carscore(Xmat, y, lambda = lam_car, diagonal = FALSE, verbose = FALSE)
      lambda_info$modo[lambda_info$metodo == "CAR"] <- "auto_estimado"
      lambda_info$lambda[lambda_info$metodo == "CAR"] <- lam_car
    } else {
      car <- care::carscore(Xmat, y, lambda = 0, diagonal = FALSE, verbose = FALSE)
      lambda_info$modo[lambda_info$metodo == "CAR"] <- "manual_lambda_0"
      lambda_info$lambda[lambda_info$metodo == "CAR"] <- 0
    }
    # Armazena CAR^2 como score.
    out$score_CAR <- as.numeric(car^2)
  }, error = function(e) {
    # Em caso de falha, preenche com NA para preservar pipeline.
    out$score_CAR <- NA_real_
    lambda_info$modo[lambda_info$metodo == "CAR"] <- "erro"
  })
  timing_seconds["score_CAR"] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # MAR (correlação marginal)
  message("[scores] Calculando MAR...")
  t0 <- Sys.time()
  tryCatch({
    # carscore com diagonal=TRUE equivale à versão marginal.
    # Mantido consistente com a política do CAR.
    if (use_auto_shrink) {
      lam_mar <- estimate_lambda_x(Xmat)
      if (!is.finite(lam_mar)) lam_mar <- 0
      mar <- care::carscore(Xmat, y, lambda = lam_mar, diagonal = TRUE, verbose = FALSE)
      lambda_info$modo[lambda_info$metodo == "MAR"] <- "auto_estimado"
      lambda_info$lambda[lambda_info$metodo == "MAR"] <- lam_mar
    } else {
      mar <- care::carscore(Xmat, y, lambda = 0, diagonal = TRUE, verbose = FALSE)
      lambda_info$modo[lambda_info$metodo == "MAR"] <- "manual_lambda_0"
      lambda_info$lambda[lambda_info$metodo == "MAR"] <- 0
    }
    out$score_MAR <- as.numeric(mar^2)
  }, error = function(e) {
    out$score_MAR <- NA_real_
    lambda_info$modo[lambda_info$metodo == "MAR"] <- "erro"
  })
  timing_seconds["score_MAR"] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # PAR (correlação parcial)
  message("[scores] Calculando PAR...")
  t0 <- Sys.time()
  tryCatch({
    # Política:
    # - se p >= 0.5*n -> shrinkage automático;
    # - senão -> lambda = 0.
    if (use_auto_shrink) {
      lam_par <- estimate_lambda_x(cbind(y, Xmat))
      if (!is.finite(lam_par)) {
        pc <- corpcor::pcor.shrink(cbind(y, Xmat), verbose = FALSE)
      } else {
        pc <- corpcor::pcor.shrink(cbind(y, Xmat), lambda = lam_par, verbose = FALSE)
      }
      par <- pc[-1, 1]
      lambda_info$modo[lambda_info$metodo == "PAR"] <- "auto_estimado"
      lambda_info$lambda[lambda_info$metodo == "PAR"] <- if (is.finite(lam_par)) lam_par else extract_lambda(pc)
    } else {
      pc <- corpcor::pcor.shrink(cbind(y, Xmat), lambda = 0, verbose = FALSE)
      par <- pc[-1, 1]
      lambda_info$modo[lambda_info$metodo == "PAR"] <- "manual_lambda_0"
      lambda_info$lambda[lambda_info$metodo == "PAR"] <- 0
    }
    out$score_PAR <- as.numeric(par^2)
  }, error = function(e) {
    out$score_PAR <- NA_real_
    lambda_info$modo[lambda_info$metodo == "PAR"] <- "erro"
  })
  timing_seconds["score_PAR"] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # I-Score marginal
  message("[scores] Calculando I-Score...")
  t0 <- Sys.time()
  tryCatch({
    # Calcula I-score SNP a SNP.
    i_scores <- vapply(seq_len(p), function(j) compute_i_score_single(Xmat[, j], y), numeric(1))
    out$score_I_Score <- as.numeric(i_scores)
  }, error = function(e) {
    out$score_I_Score <<- NA_real_
  })
  timing_seconds["score_I_Score"] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # valor-p bruto e ajustado por Bonferroni
  message("[scores] Calculando valor-p (bruto/Bonferroni)...")
  t0 <- Sys.time()
  tryCatch({
    # Monta tabela de entrada esperada por valor.p.
    genotipo_fenotipo <- data.frame(X, fenotipo = y, check.names = FALSE)
    # Calcula p-valores.
    vp <- valor.p(genotipo_fenotipo)
    out$score_ValorP_Bruto <- as.numeric(vp[, "Valor p bruto"])
    out$score_ValorP_Ajustado <- as.numeric(vp[, "Valor p ajustado"])
  }, error = function(e) {
    out$score_ValorP_Bruto <<- NA_real_
    out$score_ValorP_Ajustado <<- NA_real_
  })
  timing_seconds["score_ValorP"] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # Random Forest (uma única floresta, com duas medidas de importancia extraidas dela):
  # - score_RF_Permutation ("Importancia por Permutacao", %IncMSE): Secao 2.5 do artigo.
  # - score_RF_Impureza ("Importancia por Impureza"/MDI, IncNodePurity): Secao 2.6 do artigo.
  # Observacao (gVI): como o fenotipo aqui e continuo, a impureza de cada no e medida por
  # soma de quadrados dos residuos (RSS), nao pelo indice de Gini. Por isso nao existe um
  # "gVI" (Mean Decrease in Gini) neste pipeline -- a medida equivalente para fenotipo
  # continuo e o IncNodePurity (baseado em RSS), tratado aqui como "Importancia por Impureza".
  message("[scores] Calculando Random Forest (mtry padrao de regressao, ntree = 4000)...")
  t0 <- Sys.time()
  tryCatch({
    rf_ntree <- 4000
    rf_mtry <- max(1L, floor(ncol(X) / 3))
    rf_replace <- TRUE
    rf_sampsize <- nrow(X)
    rf_nodesize <- 5
    rf_maxnodes <- NULL
    rf_importance <- TRUE

    rf_fit <- randomForest::randomForest(
      x = X,
      y = y,
      ntree = rf_ntree,
      mtry = rf_mtry,
      replace = rf_replace,
      sampsize = rf_sampsize,
      nodesize = rf_nodesize,
      maxnodes = rf_maxnodes,
      importance = rf_importance
    )
    imp_rf <- randomForest::importance(rf_fit)
    if ("%IncMSE" %in% colnames(imp_rf)) {
      out$score_RF_Permutation <- as.numeric(imp_rf[, "%IncMSE"])
    } else {
      out$score_RF_Permutation <- as.numeric(imp_rf[, 1])
    }
    if ("IncNodePurity" %in% colnames(imp_rf)) {
      out$score_RF_Impureza <- as.numeric(imp_rf[, "IncNodePurity"])
    } else {
      out$score_RF_Impureza <- as.numeric(imp_rf[, ncol(imp_rf)])
    }
  }, error = function(e) {
    out$score_RF_Permutation <<- NA_real_
    out$score_RF_Impureza <<- NA_real_
  })
  timing_seconds["score_RF"] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # Anexa tempos como atributo para consumo no pipeline principal.
  message("[scores] Cálculo de scores concluído.")
  attr(out, "timing_seconds") <- timing_seconds
  attr(out, "lambda_info") <- lambda_info
  out  # Retorna tabela com todos os scores calculados.
}
