# Signature-count robustness: baseline confounding sim, DiffDriver ESTIMATES its BMR
# from the simulated silent count matrix via fastTopics with k signatures (true data
# has 6). Args: k N. Measures DiffDriver-estimated FPR vs the fitted k.
setwd("/dartfs-hpc/rc/home/m/f0052zm/reprotest")
suppressMessages({ library(diffdriver); library(data.table); library(Matrix) })
source("simulate_functions.R"); source("simulate_silent_estimate.R")
load("bmr_sim_inputs.rda"); load("hmm.rda"); load("codeSignature.rda"); load("hotspot1sig.rda")
Lsyn96 <- readRDS("Lsyn.rds")$Lsyn96; TARGET <- mean(readRDS("luad_char.rds")$persample_syn)
b <- bmr_sim
OUT <- "/dartfs/rc/lab/S/Szhao/jiezhou/diffdata/randombmrsingle/tsgsimple_LUAD_variations"; dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
a <- commandArgs(trailingOnly=TRUE); K <- as.integer(a[1]); N <- as.integer(a[2]); Nite <- 100L
R <- data.frame(diff_est=1, nmut=0)[rep(1,Nite),]
set.seed(1)
for (i in 1:Nite) {
  sim <- simulate_confounding(binary=TRUE, sganno=b$sganno, sgmatrix=b$sgmatrix, Nsample=N,
    para=c(0.5,0.5), beta_gc=b$beta_gc, signatures=b$signature, loadings=b$ll, rho=sqrt(0.78),
    hot=1, hmm=hmm, sc=359000)
  cohort_seed <- .Random.seed
  mut <- do.call(rbind, sim$mutlist); e <- sim$pheno; fe <- sim$efsize$diffFe; R$nmut[i] <- sum(mut)
  if (sum(mut)==0) { .Random.seed <- cohort_seed; next }
  Y <- build_silent_matrix(sim$bmrfold$bmr, Lsyn96, target_silent=TARGET)
  est <- diffdriver_estimated_mr(Y, sim$bmrmtxlist, sim$bmrfold$bmr, k=K)
  if (!is.null(est)) R$diff_est[i] <- tryCatch(ddmodel_binary(mut, e, do.call(rbind, est), fe)$pvalue, error=function(x) NA)
  .Random.seed <- cohort_seed
  cat(sprintf("k=%d N=%d i=%d nmut=%d diff_est=%.3f\n", K,N,i,R$nmut[i],R$diff_est[i])); flush.console()
}
saveRDS(R, file.path(OUT, sprintf("sigk_k%d_N%d.rds", K, N)))
cat(sprintf("DONE k=%d N=%d diff_est FP(p<.01)=%.3f meanMut=%.0f\n", K, N, mean(R$diff_est<0.01,na.rm=TRUE), mean(R$nmut)))
