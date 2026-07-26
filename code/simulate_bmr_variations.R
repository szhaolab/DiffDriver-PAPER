# BMR/confounding-sim VARIATIONS worker: one (variation, N) per call, 7 methods per iter.
# Methods: Linear/Fisher/Binomial/LogisticR/MannWhitney (count-based, run_confounding_benchmarks),
# DiffDriver + Coselens: both ESTIMATE the background from a simulated genome-wide silent
# (synonymous) count matrix (build_silent_matrix -> fastTopics / N_syn per site), exactly as in
# real analyses -- no method is handed the true per-group rate.
setwd("/dartfs-hpc/rc/home/m/f0052zm/reprotest")
suppressMessages({ library(diffdriver); library(data.table); library(Matrix); library(coselens)
                   library(brglm) })
source("simulate_functions.R"); source("run_coselens.R"); source("dndscv_truerate.R")
source("simulate_silent_estimate.R")
load("bmr_sim_inputs.rda"); load("hmm.rda"); load("codeSignature.rda"); load("hotspot1sig.rda")
cm <- readRDS("context96_map.rds")
Lsyn96 <- readRDS("Lsyn.rds")$Lsyn96; TARGET <- mean(readRDS("luad_char.rds")$persample_syn)
b <- bmr_sim
OUT <- "/dartfs/rc/lab/S/Szhao/jiezhou/diffdata/randombmrsingle/tsgsimple_LUAD_variations_est"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

# ---- variation definitions ----
VARS <- list(
  baseline    = list(binary=TRUE,  tmb_cv=0,   rho=sqrt(0.78)),
  tmb_cv0.5   = list(binary=TRUE,  tmb_cv=0.5, rho=sqrt(0.78)),
  tmb_cv0.9   = list(binary=TRUE,  tmb_cv=0.9, rho=sqrt(0.78)),
  R2_0.3      = list(binary=TRUE,  tmb_cv=0,   rho=sqrt(0.30)),
  R2_0.5      = list(binary=TRUE,  tmb_cv=0,   rho=sqrt(0.50)),
  R2_0.95     = list(binary=TRUE,  tmb_cv=0,   rho=sqrt(0.95)),
  continuous  = list(binary=FALSE, tmb_cv=0,   rho=sqrt(0.78)))

a <- commandArgs(trailingOnly=TRUE); VAR <- a[1]; N <- as.integer(a[2])
V <- VARS[[VAR]]; Nite <- if (V$binary) 100L else 50L

cos_p <- function(mr1, mr2, maf) {
  r <- tryCatch(run_coselens(maf$group1, maf$group2, subset.genes.by="ERBB3",
        true_rate=list(mutrates1=mr1, mutrates2=mr2, theta=1e6), gene_list="ERBB3", cv=NULL),
        error=function(e) NULL)
  if (is.null(r)) return(1)
  row <- r$substitutions[r$substitutions$gene_name=="ERBB3",]
  if (nrow(row)==1 && is.finite(row$pval)) row$pval else 1
}
R <- data.frame(Linear=1,Fisher=1,Binomial=1,LogisticR=1,MannWhitney=1,Coselens=1,DiffDriver=1,nmut=0)[rep(1,Nite),]
set.seed(1)
for (i in 1:Nite) {
  sim <- simulate_confounding(binary=V$binary, sganno=b$sganno, sgmatrix=b$sgmatrix, Nsample=N,
    para=c(0.5,0.5), beta_gc=b$beta_gc, signatures=b$signature, loadings=b$ll, rho=V$rho, hot=1,
    hmm=hmm, sc=359000, tmb_cv=V$tmb_cv)
  cohort_seed <- .Random.seed
  mut <- do.call(rbind, sim$mutlist); e <- sim$pheno; e_bis <- ifelse(e > mean(e), 1, 0)
  cnt <- as.numeric(Matrix::colSums(mut)); R$nmut[i] <- sum(mut)
  if (sum(mut)==0) { .Random.seed <- cohort_seed; next }
  # count-based (do not use a background rate)
  R$Linear[i]    <- tryCatch(mlr(mut,e)$pvalue, error=function(x) NA)
  R$Fisher[i]    <- tryCatch(genefisher(mut,e_bis)$pvalue, error=function(x) NA)
  R$Binomial[i]  <- tryCatch(genebinom(mut,e_bis)$pvalue, error=function(x) NA)
  R$LogisticR[i] <- tryCatch(genelr(mut,e_bis,covariates=sim$covariate)$pvalue, error=function(x) NA)
  R$MannWhitney[i] <- suppressWarnings(wilcox.test(cnt[e_bis==1], cnt[e_bis==0])$p.value)
  # ---- estimate the background from a simulated silent count matrix (shared by both methods) ----
  Y <- build_silent_matrix(sim$bmrfold$bmr, Lsyn96, target_silent=TARGET)
  fe <- sim$efsize$diffFe
  # DiffDriver (ESTIMATED BMR via fastTopics, k=6)
  R$DiffDriver[i] <- tryCatch({
    est <- diffdriver_estimated_mr(Y, sim$bmrmtxlist, sim$bmrfold$bmr, k=6)
    if (is.null(est)) NA
    else if (V$binary) ddmodel_binary(mut, e, do.call(rbind, est), fe)$pvalue
    else {  # continuous ddmodel needs a per-site label grouping equal-rate/effect sites
      seg <- rep(seq_along(sim$bmrmtxlist),
                 vapply(sim$bmrmtxlist, function(x) if (is.null(x) || length(x)==0) 0L else nrow(x), 0L))
      label <- factor(interaction(seg, fe, drop=TRUE))
      ddmodel(mut, e, do.call(rbind, est), fe, label)$pvalue
    }
  }, error=function(x) NA)
  # Coselens (ESTIMATED per-group rate = N_syn / Lsyn96); group by e>median (equal halves)
  grp <- e > median(e); g1 <- which(grp); g2 <- which(!grp)
  maf <- sim_to_maf(sim$mutlist, sim$annodata, as.integer(grp))
  if (!is.null(maf$group1) && !is.null(maf$group2) && nrow(maf$group1)>0 && nrow(maf$group2)>0) {
    mr1 <- coselens_estimated_rate(Y, g1, Lsyn96, cm)
    mr2 <- coselens_estimated_rate(Y, g2, Lsyn96, cm)
    R$Coselens[i] <- cos_p(mr1, mr2, maf)
  }
  .Random.seed <- cohort_seed
  cat(sprintf("%s N=%d i=%d nmut=%d Lin=%.3f Cos=%.3f DD=%.3f\n", VAR,N,i,R$nmut[i],R$Linear[i],R$Coselens[i],R$DiffDriver[i])); flush.console()
}
saveRDS(R, file.path(OUT, sprintf("var_%s_N%d.rds", VAR, N)))
fp <- function(p) round(mean(p<0.01, na.rm=TRUE),3)
cat(sprintf("DONE %s N=%d FP(p<.01): Lin=%.2f Fis=%.2f Bin=%.2f LR=%.2f MW=%.2f Cos=%.2f DD=%.2f meanMut=%.0f\n",
    VAR,N,fp(R$Linear),fp(R$Fisher),fp(R$Binomial),fp(R$LogisticR),fp(R$MannWhitney),fp(R$Coselens),fp(R$DiffDriver),mean(R$nmut)))
