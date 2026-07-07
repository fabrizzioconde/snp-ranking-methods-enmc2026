# =========================================================
# simulacao_base_dados_2.R
# Simula a "Base de Dados 2" descrita no artigo (main.tex, Secao 3.2,
# Equacao \ref{eq:base2}): tres efeitos aditivos isolados (SNP1, SNP2, SNP3),
# uma interacao em par (SNP4 & SNP5) e uma interacao em trio (SNP6, SNP7, SNP8).
#
#   beta0 = 640; beta1 = 2; beta2 = 1,3; beta3 = 0,9; beta4 = 2; beta5 = 3
#   L1 = (SNP1 = 1); L2 = (SNP2 = 2); L3 = (SNP3 = 3)
#   L4 = (SNP4 != 1) & (SNP5 = 2)   [par]
#   L5 = (SNP6 = 1) & (SNP7 = 2) & (SNP8 = 3)  [trio]
#
# Reproduz a Simulacao 4 de Oliveira (2015)
# (Oliveira2015): tres aditivos + interacao em par + interacao em trio.
# Equivalente ao cenario "simulacao_04" do projeto original.
# =========================================================

simular_base_dados_2 <- function(seed = 123) {
  set.seed(seed)

  num_individuos <- 1000
  num_snp <- 100

  list.snp <- list(1, 2, 3, c(4, 5), c(6, 7, 8))
  list.ia  <- list(1, 2, 3, c(-1, 2), c(1, 2, 3))

  simulacao <- scrime::simulateSNPglm(
    n.obs   = num_individuos,
    n.snp   = num_snp,
    list.ia  = list.ia,
    list.snp = list.snp,
    beta0 = 640,
    beta  = c(2, 1.3, 0.9, 2, 3),
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
