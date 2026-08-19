# beta_f sensitivity: for one gene, run the power simulation and test each cohort
# with BOTH effect-size vectors -- diffFeFix (averaged/imported DriverMAPS beta_f,
# = the paper's "DiffDriver") and diffFe (ground-truth, context-specific beta_f).
# Paired on identical cohorts (set.seed(10) per size, matching simulate_all_genes.R),
# so the averaged column reproduces the stored diffFix result and the truth column
# is the upper bound achievable by perfectly re-estimating beta_f per context.
# Args: GENE (1..200). Runs both binary-0.8 settings (TSG hot=0, OG hot=1).
setwd("/dartfs-hpc/rc/home/m/f0052zm/reprotest")
suppressMessages({ library(diffdriver); library(data.table) })
source("simulate_functions.R")
load("ri200.Rd"); load("fanno200.Rd"); load("BMR.rda"); load("hmm.rda"); load("parmASHmean.rda")
btsg <- local({ load("beta_tsg.Rd"); beta_gc })
bog  <- local({ load("beta_og.Rd");  beta_gc })
mkfix <- function(bg, idx){ f <- unlist(parmASHmean[[idx]]); names(f)[6] <- "(Intercept)"; f[names(bg)] }

a <- commandArgs(trailingOnly=TRUE); GENE <- as.integer(a[1])
SIZES <- c(200,400,600,800,1000,1200,1500); Nite <- 100L
SETTINGS <- list(
  TSG_binary_0.8 = list(beta_gc=btsg, beta_gcFix=mkfix(btsg,1), hot=0L),
  OG_binary_0.8  = list(beta_gc=bog,  beta_gcFix=mkfix(bog, 2), hot=1L))
OUT <- "/dartfs/rc/lab/S/Szhao/jiezhou/diffdata/betaf_sensitivity"; dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

signal <- GENE <= 50L
res <- list()
for (sname in names(SETTINGS)) {
  S <- SETTINGS[[sname]]
  para <- if (signal) c(0.8,0.2) else c(0.5,0.5)
  per_size <- list()
  for (N in SIZES) {
    avg <- truth <- rep(NA_real_, Nite); nm <- rep(0L, Nite)
    set.seed(10)
    for (it in 1:Nite) {
      s <- simulate_selection(binary=TRUE, sganno=ri200[[GENE]], sgmatrix=fanno200[[GENE]],
        bmrpars=log(BMR), betaf0=0, Nsample=N, beta_gc=S$beta_gc, beta_gcFix=S$beta_gcFix,
        para=para, hot=S$hot, hmm=hmm)
      mut <- do.call(rbind, s$mutlist); mr <- do.call(rbind, s$bmrmtxlist); e <- s$pheno
      nm[it] <- sum(mut); if (nm[it]==0) next
      avg[it]   <- tryCatch(ddmodel_binary_simple(mut,e,mr,s$efsize$diffFeFix)$pvalue, error=function(x) NA)
      truth[it] <- tryCatch(ddmodel_binary_simple(mut,e,mr,s$efsize$diffFe)$pvalue,    error=function(x) NA)
    }
    per_size[[as.character(N)]] <- data.frame(avg=avg, truth=truth, nmut=nm)
  }
  res[[sname]] <- per_size
  cat(sprintf("gene=%d %s signal=%s done\n", GENE, sname, signal)); flush.console()
}
saveRDS(list(gene=GENE, signal=signal, res=res),
        file.path(OUT, sprintf("betafsens_g%03d.rds", GENE)))
cat(sprintf("DONE gene=%d\n", GENE))
