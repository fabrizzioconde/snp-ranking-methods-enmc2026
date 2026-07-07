# =========================================================
# utils.R
# Utility functions and package management
# =========================================================

# ---------------------------------------------------------
# Secao: Preparacao inicial do ambiente
# Objetivo: Definir dependencias e garantir que o pipeline tenha o que precisa para rodar.
# Entrada: Nomes de pacotes obrigatorios e opcionais.
# Saida: Ambiente preparado para as etapas seguintes.
# Para que isso existe: Evita falhas por falta de bibliotecas logo no inicio.
# ---------------------------------------------------------

# Pacotes usados pelos métodos de ordenamento deste artigo:
# - obrigatórios: sem eles a simulação/análise não roda
# (Kemeny/BARD/BARDHM não fazem parte do escopo do artigo, por isso não
# aparecem aqui os pacotes ConsRank/rgl usados por eles no projeto original)
required_packages <- c("care", "corpcor", "randomForest", "scrime")
optional_packages <- character(0)

# ---------------------------------------------------------
# Secao: Verificacao e instalacao de pacotes
# Objetivo: Conferir se os pacotes do projeto estao disponiveis e instalar os ausentes.
# Entrada: Vetores com pacotes obrigatorios e opcionais.
# Saida: Lista invisivel com o status da verificacao.
# Para que isso existe: Reduz erros de execucao em maquinas diferentes.
# ---------------------------------------------------------
ensure_packages <- function(required_pkgs, optional_pkgs = character(0), repos = "https://cloud.r-project.org") {
  # Junta obrigatórios e opcionais removendo duplicados.
  all_pkgs <- unique(c(required_pkgs, optional_pkgs))
  # Lista quais obrigatórios estão faltando antes da tentativa de instalação.
  missing_required_before <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  # Lista quais opcionais estão faltando antes da tentativa de instalação.
  missing_optional_before <- optional_pkgs[!vapply(optional_pkgs, requireNamespace, logical(1), quietly = TRUE)]

  # Se faltou obrigatório, tenta instalar com dependências.
  if (length(missing_required_before) > 0) {
    message("Pacotes obrigatórios ausentes: ", paste(missing_required_before, collapse = ", "))
    install.packages(missing_required_before, repos = repos, dependencies = TRUE)
  }
  # Se faltou opcional, tenta instalar de forma leve.
  if (length(missing_optional_before) > 0) {
    message("Pacotes opcionais ausentes: ", paste(missing_optional_before, collapse = ", "))
    # Para opcionais, tenta instalação leve (evita árvore de dependências muito grande).
    install.packages(missing_optional_before, repos = repos, dependencies = FALSE)
  }
  # Mensagem amigável quando não há nada para instalar.
  if (length(missing_required_before) == 0 && length(missing_optional_before) == 0) {
    message("Todos os pacotes já estão instalados.")
  }

  # Revalida obrigatórios após tentativa de instalação.
  missing_required <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_required) > 0) {
    stop("Pacotes obrigatórios não disponíveis após tentativa de instalação: ",
         paste(missing_required, collapse = ", "))
  }

  # Revalida opcionais após tentativa de instalação.
  missing_optional <- optional_pkgs[!vapply(optional_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_optional) > 0) {
    warning("Pacotes opcionais não disponíveis: ", paste(missing_optional, collapse = ", "))
  }

  # Monta tabela de status final para inspeção no console.
  status_tbl <- data.frame(
    pacote = all_pkgs,
    instalado = vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE),
    tipo = ifelse(all_pkgs %in% required_pkgs, "obrigatorio", "opcional"),
    stringsAsFactors = FALSE
  )
  message("Status dos pacotes:\n", paste(capture.output(print(status_tbl, row.names = FALSE)), collapse = "\n"))

  invisible(list(
    required = required_pkgs,
    optional = optional_pkgs,
    missing_before = c(missing_required_before, missing_optional_before),
    missing_required = missing_required,
    missing_optional = missing_optional
  ))
}

ensure_packages(required_packages, optional_packages)

# Semente fixa para reprodutibilidade dos resultados do pipeline.
set.seed(123)

# ---------------------------------------------------------
# Secao: Conversao de score para rank inteiro
# Objetivo: Transformar qualquer vetor de score em um ranking completo de SNPs.
# Entrada: Vetor de scores, nomes dos SNPs e direcao de ordenacao.
# Saida: Vetor inteiro com posicoes 1, 2, 3, ...
# Para que isso existe: O restante do pipeline trabalha com ranks completos.
# ---------------------------------------------------------
build_integer_rank <- function(score_vector, snp_names, decreasing = TRUE) {
  # Garante tipo numérico para operação de ordenação.
  score_vector <- as.numeric(score_vector)
  
  if (length(score_vector) != length(snp_names)) {
    stop("O comprimento do vetor de scores difere do número de SNPs.")
  }
  
  # Define ordenação: decreasing=TRUE significa maior score => melhor rank.
  if (decreasing) {
    ord <- order(-score_vector, snp_names, na.last = TRUE)
  } else {
    ord <- order(score_vector, snp_names, na.last = TRUE)
  }
  
  # Obtém nomes dos SNPs já em ordem de prioridade.
  ordered_snps <- snp_names[ord]
  
  # Inicializa vetor de ranks inteiros.
  rank_vector <- integer(length(snp_names))
  # Nomeia vetor por SNP para facilitar indexação por nome.
  names(rank_vector) <- ordered_snps
  # Atribui posições 1,2,3,... na ordem calculada.
  rank_vector[ordered_snps] <- seq_along(ordered_snps)
  
  # Reordena para a ordem original dos SNPs de entrada.
  rank_vector <- rank_vector[snp_names]
  # Retorna ranks como inteiros.
  return(as.integer(rank_vector))
}

# ---------------------------------------------------------
# Secao: Identificacao de SNPs causais
# Objetivo: Extrair os nomes dos SNPs causais a partir da definicao do cenario.
# Entrada: Estrutura `list_snp` com os indices causais.
# Saida: Vetor com nomes como `SNP1`, `SNP2`, etc.
# Para que isso existe: Facilita comparar o resultado do metodo com a verdade simulada.
# ---------------------------------------------------------
# Extrai nomes dos SNPs causais a partir de list.snp (ex.: list(1,2,c(3,4)) -> SNP1, SNP2, SNP3, SNP4).
# Preserva a ordem de aparição no cenário.
extract_causal_snp_names <- function(list_snp) {
  idx <- unique(unlist(list_snp))
  paste0("SNP", idx)
}

# ---------------------------------------------------------
# Secao: Tabela formal dos SNPs causais
# Objetivo: Organizar a definicao causal em formato de tabela.
# Entrada: Estrutura `list_snp` do cenario.
# Saida: Data frame com termo causal, ordem e nome do SNP.
# Para que isso existe: Melhora a rastreabilidade e a exportacao dos cenarios simulados.
# ---------------------------------------------------------
# Constrói dataframe com definição dos SNPs causais por termo do cenário.
build_causal_definition_df <- function(list_snp) {
  parts <- vector("list", length(list_snp))
  for (k in seq_along(list_snp)) {
    idx <- as.integer(list_snp[[k]])
    parts[[k]] <- data.frame(
      termo_causal = k,
      ordem_no_termo = seq_along(idx),
      snp_idx = idx,
      SNP = paste0("SNP", idx),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, parts)
}

# ---------------------------------------------------------
# Secao: Exportacao de tabela em LaTeX
# Objetivo: Salvar um data frame no formato de tabela LaTeX.
# Entrada: Data frame, caminho do arquivo e metadados da tabela.
# Saida: Arquivo `.tex` pronto para ser usado em relatorios.
# Para que isso existe: Facilita levar os resultados para artigos e documentos.
# ---------------------------------------------------------
# Salva um data.frame como tabela LaTeX (tabular) em um arquivo .tex.
write_df_latex <- function(df, file_path, caption = "", label = "tab:ranks", digits = 0) {
  # Se não há linhas, não escreve arquivo.
  if (nrow(df) == 0) return(invisible(NULL))
  # Número de colunas.
  nc <- ncol(df)
  # Alinhamento: primeira coluna à esquerda, restantes à direita.
  col_align <- c("l", rep("r", nc - 1))
  # Cabeçalho da tabela.
  header <- paste(colnames(df), collapse = " & ")
  # Linhas iniciais do ambiente tabular.
  lines <- c(
    "\\begin{table}[ht]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    paste0("\\begin{tabular}{", paste(col_align, collapse = ""), "}"),
    "\\hline",
    paste0(header, " \\\\"),
    "\\hline"
  )
  # Percorre linhas do data.frame para escrever corpo da tabela.
  for (i in seq_len(nrow(df))) {
    # Vetor de strings formatadas para a linha i.
    row_vals <- character(nc)
    for (j in seq_len(nc)) {
      v <- df[i, j]  # Valor bruto da célula.
      # Formata números com arredondamento; texto segue como está.
      if (is.numeric(v)) row_vals[j] <- format(round(v, digits), nsmall = digits) else row_vals[j] <- as.character(v)
    }
    lines <- c(lines, paste0(paste(row_vals, collapse = " & "), " \\\\"))
  }
  lines <- c(lines, "\\hline", "\\end{tabular}", "\\end{table}")
  writeLines(lines, file_path)
}
