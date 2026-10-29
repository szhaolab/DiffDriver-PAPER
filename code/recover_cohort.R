suppressMessages({library(diffdriver); library(data.table); library(Matrix)})
ns <- function(x) get(x, envir=asNamespace("diffdriver"))
matrixlistToBMR <- ns("matrixlistToBMR"); prep_positional_data <- ns("prep_positional_data")
prep_pheno_mut_data <- ns("prep_pheno_mut_data"); mut2hotspot <- ns("mut2hotspot")
mlr <- ns("mlr"); mlr.v2 <- ns("mlr.v2"); genefisher <- ns("genefisher")
genebinom <- ns("genebinom"); genelr <- ns("genelr")
aheader <- ns("aheader"); acoltype <- ns("acoltype"); OGs <- ns("OGs")
OGpars <- ns("OGpars"); TSGpars <- ns("TSGpars"); hmm <- ns("hmm")

args <- commandArgs(trailingOnly=TRUE)
outbase <- "/dartfs-hpc/rc/home/m/f0052zm/szhao_lab_share/diffDriver/data_run_pancancer299"
man <- fread(file.path(outbase,"manifest_recover.txt"))
tid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
row <- man[tid]; cid <- row$cohort; kk <- row$k
tset <- strsplit(row$tumors, ",")[[1]]
cdir <- file.path(outbase, cid); dir.create(cdir, showWarnings=FALSE)
anno_dir <- "/dartfs/rc/lab/S/Szhao/library/diffdriver_anno/annodir96"; totalnttype <- 96
afileinfo <- list(file=file.path(anno_dir, paste0("anno",totalnttype,"_nttypeXXX_annodata.txt")),
                  header=aheader, coltype=acoltype, totalntype=totalnttype)
gene0 <- readLines(file.path(cdir,"genes.txt"))
companion <- if ("TP53" %in% gene0) "PIK3CA" else "TP53"
gene <- unique(c(gene0, companion))   # pad to avoid single-gene prep_positional_data bug
## shared master files, filter to this cohort's tumor set
allsamp <- fread(file.path(outbase,"master_pheno.txt"))            # SampleID(12), Study, C1..C4, Age, Sex
samp <- allsamp[Study %in% tset]
setnames(samp, "Study", "Cohort")
allmut <- fread(file.path(outbase,"pancan_all_muts.txt"))          # Chromosome,Position,Ref,Alt,SampleID(12),Tumor
mut <- as.data.frame(allmut[Tumor %in% tset & SampleID %in% samp$SampleID, .(Chromosome,Position,Ref,Alt,SampleID)])
## covariates: Age always; drop Sex if single-sex cohort; drop Cohort if single tumor type
covcols <- "Age"
if (length(unique(samp$Sex)) > 1)    covcols <- c(covcols, "Sex")
if (length(unique(samp$Cohort)) > 1) covcols <- c(covcols, "Cohort")
cov0  <- as.data.frame(samp[, c("SampleID", covcols), with=FALSE])
cat(">>> covariates:", paste(covcols, collapse=", "), "\n")
set.seed(0)
cat(">>>", cid, "k=", kk, "nsamp=", nrow(samp), "genes=", length(gene), "START", format(Sys.time()), "\n"); flush.console()

## ---- BMR once ----
bmrfile <- file.path(cdir, paste0(cid,"_BMRres.Rdata"))
if (file.exists(bmrfile)) { cat(">>> loading saved BMR\n"); load(bmrfile) } else {
  BMRres <- matrixlistToBMR(afileinfo, mut=mut, "signature", k=kk, file.path(cdir, cid)) }
BMRreg <- BMRres[[1]]; bmrsig <- BMRres[[2]]
rdata <- prep_positional_data(gene, afileinfo, BMRreg=BMRreg, output_prefix=cid, output_dir=cdir)
fanno <- rdata$fanno; ri <- rdata$ri
nsyndt <- data.table("SampleID"=names(BMRreg[["ysample"]]), "Nsyn"=BMRreg[["ysample"]])
## establish sample set/order via prep (use C1 to define phename slot); covariates same across phenos
phe1 <- as.data.frame(samp[, .(SampleID, C1)])
cdata <- prep_pheno_mut_data(mut, phe1, covariates=cov0, add_dt=nsyndt)
mutd <- cdata$mut; ci <- cdata$ci; pheno <- cdata$pheno; cov <- cdata$cov
hotspots <- mut2hotspot(mutd)
## signature bmr matrix (phenotype-independent)
bmr <- data.table()
for (i in 1:nrow(ci)){ id <- ci[i,"SampleID"][[1]]
  glambda <- bmrsig$genesp[ri$genename, lambda]; glambda[which(is.na(glambda))] <- 1
  mui <- bmrsig$ysample[id]*bmrsig$sigmtx[id,ri$nttypecode]*exp(as.matrix(fanno[,c("expr","repl","hic")])%*%bmrsig$BMvbeta)*glambda/bmrsig$at[ri$nttypecode]
  bmr <- cbind(bmr, log(mui)) }
if (any(is.na(bmr))) stop("bmr missing")
bmrallg <- split(bmr, ri$genename); riallg <- split(ri, ri$genename); fannoallg <- split(fanno, ri$genename)
## map sample -> each phenotype value, aligned to ci order
phmat <- as.data.frame(samp)[match(ci$SampleID, samp$SampleID), c("C1","C2","C3","C4")]

## ---- loop phenotypes, reuse everything above ----
for (cx in c("C1","C2","C3","C4")){
  if (file.exists(file.path(cdir, paste0(cid,"_",cx,"_resdd.txt")))) { cat(">>> skip existing", cx, "\n"); next }
  e_all <- phmat[[cx]]
  res <- list()
  for (g in names(bmrallg)){
    r <- tryCatch({
      rig <- riallg[[g]]; rig$ridx <- 1:dim(rig)[1]
      muti <- na.omit(ci[rig[mutd, on=c("chrom"="Chromosome","start"="Position","ref"="Ref","alt"="Alt")], on="SampleID"])
      mutmtx <- Matrix::sparseMatrix(i=muti$ridx, j=muti$cidx, dims=c(max(rig$ridx), max(ci$cidx)))
      hotg <- na.omit(rig[hotspots, on=c("chrom"="Chromosome","start"="Position")]); hotmat <- rep(0,nrow(rig)); hotmat[hotg$ridx] <- 1
      if (sum(mutmtx)==0) next
      bmrmtx <- bmrallg[[g]]; ganno <- fannoallg[[g]]
      if (g %in% OGs[,1]){ betaf <- OGpars[names(OGpars)!="beta_f0"]; betaf0 <- OGpars["beta_f0"]
        fe <- as.matrix(ganno[,names(betaf),with=F])%*%betaf + hotmat*hmm[9] + betaf0
      } else { betaf <- TSGpars[names(TSGpars)!="beta_f0"]; betaf0 <- TSGpars["beta_f0"]
        fe <- as.matrix(ganno[,names(betaf),with=F])%*%betaf + betaf0 }
      label <- factor(1:nrow(bmrmtx))
      resg <- list()
      resg[["dd"]] <- ddmodel(mutmtx, e_all, bmrmtx, fe[,1], label=label, cov=cov)
      resg[["mlr"]] <- mlr(mutmtx, e_all); resg[["mlr.v2"]] <- mlr.v2(mutmtx, e_all, pheno$Nsyn)
      eb <- ifelse(e_all>mean(e_all),1,0)
      resg[["fisher"]] <- genefisher(mutmtx, eb); resg[["binom"]] <- genebinom(mutmtx, eb); resg[["lr"]] <- genelr(mutmtx, eb)
      resg
    }, error=function(err){ cat("  [skip]", cx, g, ":", conditionMessage(err), "\n"); NULL })
    if (!is.null(r)) res[[g]] <- r
  }
  if (length(res)==0){ cat("  no genes for", cx, "\n"); next }
  meth <- c("dd","mlr","mlr.v2","fisher","binom","lr")
  resdf <- data.frame(lapply(meth, function(x) unlist(lapply(lapply(res,'[[',x),'[[','pvalue'))))
  colnames(resdf) <- paste0(meth,".p")
  resdf[,paste0(meth,".fdr")] <- apply(resdf,2,p.adjust,method="fdr")
  resdf[,c("mut.E1","mut.E0","E1","E0")] <- do.call(rbind, lapply(lapply(res,'[[',"fisher"),'[[',"count"))
  resdf$alpha_e <- sapply(rownames(resdf), function(g) tryCatch(res[[g]][["dd"]][["res.alt"]]$alpha[2], error=function(x) NA))
  resdf <- resdf[rownames(resdf) %in% gene0, , drop=FALSE]
  write.table(resdf, file.path(cdir, paste0(cid,"_",cx,"_resdd.txt")), quote=F)
  cat(">>>", cx, "done:", nrow(resdf), "genes", format(Sys.time()), "\n"); flush.console()
}
cat(">>> COHORT DONE", format(Sys.time()), "\n")
