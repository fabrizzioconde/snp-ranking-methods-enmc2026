# =========================================================
# simulacao_base_dados_1.R
# Simula a "Base de Dados 1" do artigo: Simulacao 1 de Oliveira (2015)
# e cenario "simulacao_01" em ranks_simulacao_somente_inteiros/
# (simulacao_01_bard.R).
#
# Oito SNPs causais independentes (SNP1 a SNP8), sem interacoes epistaticas,
# com modelos geneticos codificados por list.ia:
#   +1 aditivo, +2 dominante, +3 recessivo,
#   -1 aditivo reverso, -2 dominante reverso, -3 recessivo reverso.
#
# Parametros: beta0 = 640; beta = rep(1, 8); maf = c(0.1, 0.4); seed = 123.
# Fenotipo continuo (err.fun = rnorm), como no pipeline de ranks.
# =========================================================

simular_base_dados_1 <- function(seed = 123) {
  set.seed(seed)

  num_individuos <- 1000
  num_snp <- 100

  list.snp <- list(1, 2, 3, 4, 5, 6, 7, 8)
  list.ia  <- list(1, 2, 3, -1, -2, -3, 1, 2)

  simulacao <- scrime::simulateSNPglm(
    n.obs   = num_individuos,
    n.snp   = num_snp,
    list.ia  = list.ia,
    list.snp = list.snp,
    beta0 = 640,
    beta  = c(1, 1, 1, 1, 1, 1, 1, 1),
    maf   = c(0.1, 0.4),
    err.fun = rnorm,
    rand = seed
  )

  genotipo <- as.data.frame(simulacao$x)
  colnames(genotipo) <- paste0("SNP", seq_len(ncol(genotipo)))
  fenotipo <- simulacao$y

  list(
    X = genotipo,
    y = fenotipo,
    list_snp_causal = list.snp,
    snps_causais = paste0("SNP", unique(unlist(list.snp)))
  )
}
