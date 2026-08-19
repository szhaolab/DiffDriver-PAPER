DIR <- "/private/tmp/claude-501/-Users-siming-Dartmouth-College-Dropbox-Siming-Zhao-Project-diffDriver-DiffDriver-PAPER/46eb8cab-f9ce-420f-8359-609ff7701606/scratchpad/betafsens_rds"
DIR <- "betafsens_rds"   # per-gene betafsens_gNNN.rds from simulate_betaf_sensitivity.R
genes <- lapply(1:NG, function(g) readRDS(file.path(DIR, sprintf("betafsens_g%03d.rds", g))))
settings <- names(genes[[1]]$res)

metrics <- function(mat){         # mat: 200 x Niter p-values
  ni <- ncol(mat); sens <- fpr <- numeric(ni)
  for (j in 1:ni){
    q <- p.adjust(mat[,j], method="BH"); hits <- which(q < FDR)
    if (length(hits)>0){ sens[j] <- sum(hits<=NT)/NT; fpr[j] <- sum(hits>NT)/length(hits) }
  }
  c(power=mean(sens), fdr=mean(fpr))
}
out <- data.frame()
for (sn in settings){
  for (N in SIZES){
    ni <- nrow(genes[[1]]$res[[sn]][[as.character(N)]])
    A <- T <- matrix(NA, NG, ni)
    for (g in 1:NG){ d <- genes[[g]]$res[[sn]][[as.character(N)]]; A[g,] <- d$avg; T[g,] <- d$truth }
    ma <- metrics(A); mt <- metrics(T)
    out <- rbind(out, data.frame(setting=sn, N=N,
      power_avg=round(ma["power"],3), power_truth=round(mt["power"],3), dpow=round(mt["power"]-ma["power"],3),
      fdr_avg=round(ma["fdr"],3), fdr_truth=round(mt["fdr"],3)))
  }
}
rownames(out) <- NULL
print(out, row.names=FALSE)
cat("\n== mean power gain (truth - averaged) by setting ==\n")
for (sn in settings) cat(" ", sn, ":", round(mean(out$dpow[out$setting==sn]),3), "\n")
saveRDS(out, file.path(dirname(DIR), "betaf_sensitivity_metrics.rds"))
