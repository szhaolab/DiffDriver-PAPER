suppressMessages(library(data.table))
outbase <- "/dartfs-hpc/rc/home/m/f0052zm/szhao_lab_share/diffDriver/data_run_pancancer299"
man <- fread(file.path(outbase,"manifest299.txt"))
long <- list()
for (i in 1:nrow(man)){
  cid <- man$cohort[i]; tum <- man$tumors[i]; ntum <- man$ntumor[i]; cdir <- file.path(outbase, cid)
  for (cx in c("C1","C2","C3","C4")){
    f <- file.path(cdir, paste0(cid,"_",cx,"_resdd.txt"))
    if (!file.exists(f)) next
    rd <- tryCatch(read.table(f, header=TRUE), error=function(e) NULL)
    if (is.null(rd) || nrow(rd)==0) next
    rd$gene <- rownames(rd)
    long[[paste(cid,cx)]] <- data.table(gene=rd$gene, phenotype=cx, cohort=cid, n_tumor_types=ntum, tumors=tum,
      dd.p=rd$dd.p, alpha_e=if("alpha_e" %in% names(rd)) rd$alpha_e else NA,
      mut.E1=rd$mut.E1, mut.E0=rd$mut.E0, E1=rd$E1, E0=rd$E0)
  }
}
res <- rbindlist(long)
res[, dd.fdr := p.adjust(dd.p, method="fdr"), by=phenotype]
res[, direction := ifelse(is.na(alpha_e),"NA", ifelse(alpha_e>0,"up_in_subtype","up_in_rest"))]
## cohort-level mutation frequency (constant across phenotypes for a gene)
res[, n_samples := E1 + E0]
res[, n_mutated := mut.E1 + mut.E0]
res[, mut_freq := round(n_mutated / n_samples, 4)]
setcolorder(res, c("gene","phenotype","dd.p","dd.fdr","direction","n_tumor_types","tumors","n_samples","n_mutated","mut_freq","mut.E1","mut.E0","E1","E0","alpha_e"))
setorder(res, phenotype, dd.p)
fwrite(res, file.path(outbase,"pancancer299_immune_aggregated.txt"), sep="\t")
## per-gene cohort summary (one row per gene)
gsum <- unique(res[, .(gene, n_tumor_types, tumors, n_samples, n_mutated, mut_freq)])
setorder(gsum, -mut_freq)
fwrite(gsum, file.path(outbase,"pancancer299_gene_cohorts.txt"), sep="\t")
cat("aggregated rows:", nrow(res), " genes:", length(unique(res$gene)), "\n")
cat("gene cohort summary rows:", nrow(gsum), "\n")
cat("mut_freq range:", min(gsum$mut_freq), "-", max(gsum$mut_freq), "\n")
