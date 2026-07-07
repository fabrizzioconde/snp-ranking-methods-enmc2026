# =========================================================
# gerar_tabelas_artigo.R
# Simula as Bases de Dados 1 e 2 do artigo, calcula os 6 metodos de
# ordenamento (Secao 2 do main.tex) e exporta CSVs, fragmentos LaTeX e
# figuras (histograma e boxplot) para cada base.
#
# Como executar (a partir desta pasta):
#   Rscript gerar_tabelas_artigo.R
# =========================================================

script_dir <- tryCatch(
  dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))),
  error = function(e) getwd()
)
if (is.null(script_dir) || !nzchar(script_dir) || !dir.exists(script_dir)) script_dir <- getwd()
setwd(script_dir)

source("utils.R")
source("scoring_methods.R")
source("simulacao_base_dados_1.R")
source("simulacao_base_dados_2.R")

dir_saida <- file.path(script_dir, "resultados")
dir_figuras <- file.path(dirname(script_dir), "figures")
if (!dir.exists(dir_saida)) dir.create(dir_saida, recursive = TRUE)
if (!dir.exists(dir_figuras)) dir.create(dir_figuras, recursive = TRUE)

metodos_tab1 <- c(
  "Correlação Marginal"        = "rank_Marginal",
  "Correlação Parcial"         = "rank_Parcial",
  "CAR Score"                  = "rank_CAR",
  "I Score"                    = "rank_IScore",
  "Importância por Permutação" = "rank_Permutation",
  "Importância por Impureza"   = "rank_Impureza"
)

fmt_snp <- function(snp, snps_causais) {
  if (snp %in% snps_causais) paste0("\\textbf{", snp, "}") else snp
}

configurar_par_grafico <- function() {
  graphics::par(
    mar = c(5.8, 5.8, 4.8, 2.2),
    cex.lab = 1.9,
    cex.axis = 1.6,
    cex.main = 1.7,
    mgp = c(3.8, 1.1, 0)
  )
}

plotar_boxplot_destacado <- function(y, titulo) {
  stats <- grDevices::boxplot.stats(y)
  outliers <- stats$out

  y_lim <- range(c(y, outliers))
  y_pad <- diff(y_lim) * 0.06

  configurar_par_grafico()
  graphics::boxplot(
    y,
    main = titulo,
    ylab = "Fenótipo simulado",
    col = "lightblue",
    border = "gray40",
    lwd = 1.3,
    outline = FALSE,
    staplelwd = 1.2,
    whisklty = 1,
    whisklwd = 1.2,
    medlwd = 2,
    medcol = "gray20",
    ylim = c(y_lim[1] - y_pad, y_lim[2] + y_pad),
    xlim = c(0.35, 1.65)
  )

  if (length(outliers) > 0) {
    n_out <- length(outliers)
    idx <- order(outliers)
    x_out <- 1 + (seq_len(n_out) - (n_out + 1) / 2) * 0.24
    pt_cex <- max(1.4, 2.4 - 0.18 * n_out)
    graphics::points(
      x_out,
      outliers[idx],
      pch = 21,
      col = "#7B241C",
      bg = "#E74C3C",
      cex = pt_cex,
      lwd = 1.8
    )
    graphics::mtext(
      paste0("Outliers: n = ", n_out),
      side = 3,
      line = 0.2,
      adj = 1,
      cex = 1.2,
      col = "#7B241C",
      font = 2
    )
  }

  invisible(stats)
}

gerar_figuras_fenotipo <- function(y, base_num) {
  titulo <- paste0("Base de Dados ", base_num)
  arquivo_hist <- file.path(dir_figuras, paste0("histograma_base", base_num, ".jpeg"))
  arquivo_box <- file.path(dir_figuras, paste0("boxplot_base", base_num, ".jpeg"))

  grDevices::jpeg(arquivo_hist, width = 900, height = 700, quality = 95)
  configurar_par_grafico()
  graphics::hist(
    y,
    main = titulo,
    xlab = "Fenótipo simulado",
    ylab = "Número de indivíduos",
    col = "lightblue",
    border = "white"
  )
  grDevices::dev.off()

  grDevices::jpeg(arquivo_box, width = 700, height = 900, quality = 95)
  plotar_boxplot_destacado(y, titulo)
  grDevices::dev.off()

  list(histograma = arquivo_hist, boxplot = arquivo_box)
}

processar_base <- function(sim_fun, base_num, tab_top20_label, tab_causais_label) {
  message("[artigo] Simulando Base de Dados ", base_num, "...")
  sim <- sim_fun(seed = 123)
  X <- sim$X
  y <- sim$y
  snps_causais <- sim$snps_causais
  message("[artigo] SNPs causais (Base ", base_num, "): ", paste(snps_causais, collapse = ", "))

  gerar_figuras_fenotipo(y, base_num)

  scores <- compute_all_scores(X, y)
  snp_names <- scores$SNP

  scores$rank_Marginal    <- build_integer_rank(scores$score_MAR, snp_names, decreasing = TRUE)
  scores$rank_Parcial     <- build_integer_rank(scores$score_PAR, snp_names, decreasing = TRUE)
  scores$rank_CAR         <- build_integer_rank(scores$score_CAR, snp_names, decreasing = TRUE)
  scores$rank_IScore      <- build_integer_rank(scores$score_I_Score, snp_names, decreasing = TRUE)
  scores$rank_Permutation <- build_integer_rank(scores$score_RF_Permutation, snp_names, decreasing = TRUE)
  scores$rank_Impureza    <- build_integer_rank(scores$score_RF_Impureza, snp_names, decreasing = TRUE)
  scores$causal <- scores$SNP %in% snps_causais

  prefixo <- paste0("base_dados_", base_num)
  write.csv(scores, file.path(dir_saida, paste0("scores_e_ranks_", prefixo, ".csv")), row.names = FALSE)

  ranks_causais <- scores[scores$causal,
    c("SNP", "rank_Marginal", "rank_Parcial", "rank_CAR", "rank_IScore", "rank_Permutation", "rank_Impureza")]
  ranks_causais <- ranks_causais[order(as.integer(sub("SNP", "", ranks_causais$SNP))), ]
  write.csv(ranks_causais, file.path(dir_saida, paste0("ranks_causais_", prefixo, ".csv")), row.names = FALSE)

  message("[artigo] Ranks dos SNPs causais (Base ", base_num, "):")
  print(ranks_causais, row.names = FALSE)

  tabela_top20 <- sapply(names(metodos_tab1), function(label) {
    col <- metodos_tab1[[label]]
    ordenado <- scores[order(scores[[col]]), ]
    head(ordenado$SNP, 20)
  })

  linhas_tab1 <- apply(tabela_top20, 1, function(linha) {
    paste(vapply(linha, fmt_snp, character(1), snps_causais = snps_causais), collapse = " & ")
  })

  tex_tab1 <- c(
    "\\begin{table}[H]",
    "\\centering",
    paste0("\\caption{Rank dos SNPs por seis métodos de ordenamento para os 20 primeiros SNPs na Base de Dados ", base_num, ".}"),
    paste0("\\label{", tab_top20_label, "}"),
    "\\small",
    "\\begin{tabular}{l|l|l|l|l|l}",
    "  \\hline",
    paste0(
      "\\begin{tabular}[c]{@{}c@{}}Correlação\\\\Marginal\\end{tabular} & ",
      "\\begin{tabular}[c]{@{}c@{}}Correlação\\\\Parcial\\end{tabular} & ",
      "CAR Score & I Score & ",
      "\\begin{tabular}[c]{@{}c@{}}Importância\\\\por Permutação\\end{tabular} & ",
      "\\begin{tabular}[c]{@{}c@{}}Importância\\\\por Impureza\\end{tabular} \\\\"
    ),
    "  \\hline",
    paste0("  ", linhas_tab1, " \\\\"),
    "   \\hline",
    "\\end{tabular}",
    "\\end{table}"
  )
  writeLines(tex_tab1, file.path(dir_saida, paste0("tabela", base_num, "_top20.tex")))

  linhas_tab2 <- apply(ranks_causais, 1, function(row) {
    paste(row["SNP"],
          sprintf("%3d", as.integer(row["rank_Marginal"])),
          sprintf("%3d", as.integer(row["rank_Parcial"])),
          sprintf("%3d", as.integer(row["rank_CAR"])),
          sprintf("%3d", as.integer(row["rank_IScore"])),
          sprintf("%3d", as.integer(row["rank_Permutation"])),
          sprintf("%3d", as.integer(row["rank_Impureza"])),
          sep = " & ")
  })

  tex_tab2 <- c(
    "\\begin{table}[H]",
    "\\centering",
    paste0("\\caption{Posição dos SNPs causais por seis métodos de ordenamento na Base de Dados ", base_num, ".}"),
    paste0("\\label{", tab_causais_label, "}"),
    "\\small",
    "\\begin{tabular}{l|r|r|r|r|r|r}",
    "  \\hline",
    paste0(
      "SNPs & ",
      "\\begin{tabular}[c]{@{}c@{}}Correlação\\\\Marginal\\end{tabular} & ",
      "\\begin{tabular}[c]{@{}c@{}}Correlação\\\\Parcial\\end{tabular} & ",
      "CAR Score & I Score & ",
      "\\begin{tabular}[c]{@{}c@{}}Importância\\\\por Permutação\\end{tabular} & ",
      "\\begin{tabular}[c]{@{}c@{}}Importância\\\\por Impureza\\end{tabular} \\\\"
    ),
    "  \\hline",
    paste0("  ", linhas_tab2, " \\\\"),
    "   \\hline",
    "\\end{tabular}",
    "\\end{table}"
  )
  writeLines(tex_tab2, file.path(dir_saida, paste0("tabela", base_num, "_posicao_causais.tex")))

  invisible(list(scores = scores, ranks_causais = ranks_causais))
}

processar_base(simular_base_dados_1, 1L, "tab 33", "tab 34")
processar_base(simular_base_dados_2, 2L, "tab 35", "tab 36")

message("[artigo] Concluido. Arquivos gravados em: ", dir_saida)
message("[artigo] Figuras gravadas em: ", dir_figuras)
