#' Generate hotspot indicator sequence via HMM
#'
#' @param hmm Numeric vector of 10 HMM parameters (transition probs and coefficients)
#' @param sgdata List of annotation data.frames (one per nucleotide type)
#' @return A matrix with columns: start (genomic position), seqt (hotspot indicator 0/1)
hotspotseq <- function(hmm, sgdata){
  pos <- sort(unique(do.call(rbind, sgdata)$start))
  seqt <- sample(x=c(0,1), size=1, prob=c(hmm[1], hmm[2]))
  for (i in 2:length(pos)) {
    a <- ifelse(seqt[i-1]==0,
                sample(x=c(0,1), size=1, prob=c(hmm[3], hmm[4])),
                sample(x=c(0,1), size=1, prob=c(hmm[5], hmm[6])))
    seqt <- c(seqt, a)
  }
  cbind(start=pos, seqt=seqt)
}

#' Simulate mutations under a phenotype-dependent selection model (1 functional variable)
#'
#' Ported from szhaolab/diffdriver scripts/simulate_functions.R::simulate_1funcvi
#' (kept under the name simulate_selection for the pipeline).
#' Samples are split into a positively-selected group (size Nsample.ps) and a
#' neutral group (Nsample.neu). The neutral background mutation probability per
#' site of nucleotide type t is exp(bmrpars[t]) * exp(betaf0) and is CONSTANT
#' across samples (no per-sample BMR fold; hence no rho/tau). Selected samples
#' additionally carry the functional fold from beta_gc. `phenotype` is reordered
#' so that selected samples come first, matching the mutation-matrix columns.
#'
#' @param binary Logical; TRUE for binary (case/control), FALSE for continuous phenotype
#' @param sganno List of 9 annotation data.frames (one per nucleotide type)
#' @param sgmatrix List of 9 functional annotation matrices
#' @param bmrpars Log-scale background mutation rates (length 9)
#' @param betaf0 Baseline functional effect size (log-scale)
#' @param Nsample Total number of samples
#' @param beta_gc Named vector of functional effect coefficients (variable effect)
#' @param beta_gcFix Named vector of functional effect coefficients (fixed effect;
#'   defaults to beta_gc). Used to produce the diffFeFix effect-size vector.
#' @param para For binary: c(selection_prob_case, selection_prob_control).
#'   For continuous: c(mean, sd, intercept, slope) where intercept and slope
#'   are logistic regression coefficients for selection probability given phenotype.
#' @param hot Hotspot flag; 0 = no hotspot effect
#' @param hmm HMM parameters for hotspot sequence generation
#' @return A list containing mutlist, pheno, foldlist, annodata, bmrpars,
#'   bmrmtxlist, para, efsize (betaf0, beta_gc, avFe, diffFe, diffFeFix), nsample
simulate_selection <- function(binary=F, sganno, sgmatrix, bmrpars, betaf0=2, Nsample,
                               beta_gc, beta_gcFix=beta_gc, para, hot=0, hmm){
  if (binary==T){ # generate binary phenotype
    Nsamplec <- round(Nsample/2)
    Nsamplen <- Nsample-Nsamplec
    phenotype <- c(rep(1,Nsamplec),rep(0,Nsamplen))
    ss=ifelse(phenotype==1,sample(c(0,1),size=Nsamplec,replace=T,prob = c(1-para[1],para[1])),sample(c(0,1),size=Nsamplen,replace=T,prob = c(1-para[2],para[2])))
  }else{
    phenotype=rnorm(Nsample,mean = para[1],sd=para[2])
    pp=exp(para[3]+para[4]*phenotype)/(1+exp(para[3]+para[4]*phenotype))
    ss=ifelse(runif(Nsample)-pp<0,1,0)
  }
  index=which(ss==1)
  phenotype=c(phenotype[index],phenotype[-index])
  Nsample.ps=sum(ss)
  Nsample.neu <- Nsample - Nsample.ps
  hotseq=hotspotseq(hmm,sganno)
  if (hot==0){
    hotseq[,2]=0
  }
  mutlist <- list()
  countlist <- list()
  annodata <- list()
  bmrmtxlist <- list()
  betagc=c(beta_gc,hmm[9])
  betagcFix=c(beta_gcFix,hmm[9])
  if (all( names(beta_gc) != names(beta_gcFix)) ){stop("beta_gc does not match beta_gcFix!")}
  mutRate <- list()
  foldlist <- list()
  foldlistFix <- list()
  modelmatrix <- list()

  for (t in 1:length(sganno)) {
    hotseqt=merge(sganno[[t]],hotseq,by="start")$seqt
    selename=names(beta_gc)
    ssgdata=cbind(sgmatrix[[t]][,..selename],hotseqt)
    hotindex=which(hotseqt==1)
    pp.neu=rep(exp(bmrpars[t])*exp(betaf0),nrow(ssgdata))
    fold=exp(as.matrix(ssgdata)%*%betagc)
    fold=(Nsample/Nsample.ps)*fold-Nsample.neu/Nsample.ps
    fold[hotindex]=exp(hmm[9])
    if (any(fold<0)){warning("Fold is negative!")}
    fold[which(fold<0)]=1
    pp.ps=ifelse(pp.neu*fold<1,pp.neu*fold,1)
    foldlist[[t]]=data.table::data.table(fold=fold)

    foldFix=exp(as.matrix(ssgdata)%*%betagcFix)
    foldFix[hotindex]=exp(hmm[9])
    foldFix=(Nsample/Nsample.ps)*foldFix-Nsample.neu/Nsample.ps
    if (any(foldFix<0)){warning("Foldfix is negative!")}
    foldFix[which(foldFix<0)]=1
    foldlistFix[[t]]=data.table::data.table(fold=foldFix)

    mutps=replicate(Nsample.ps,rbinom(length(pp.ps),size=1,pp.ps))
    mutneu=replicate(Nsample.neu,rbinom(length(pp.neu),size=1,pp.neu))
    if (!is.matrix(mutps)){
      mutlist[[t]]=as(cbind(t(mutps),t(mutneu)),"sparseMatrix")
    }else{
      mutlist[[t]]=as(cbind(mutps,mutneu),"sparseMatrix")
    }
    countlist[[t]] <- c(sum(mutlist[[t]]))
    bmrmtxlist[[t]] <- matrix(bmrpars[t]+betaf0, ncol = ncol(mutlist[[t]]), nrow = nrow(mutlist[[t]]))
    modelmatrix[[t]] <- ssgdata
  }

  fold <- do.call(rbind,foldlist)
  foldFix <- do.call(rbind,foldlistFix)
  avFe <- rep(log(mean(fold[[1]])*Nsample.ps/Nsample + Nsample.neu/Nsample),nrow(fold))
  diffFe <-  log(fold[[1]]*Nsample.ps/Nsample + Nsample.neu/Nsample)
  diffFeFix <-  log(foldFix[[1]]*Nsample.ps/Nsample + Nsample.neu/Nsample)

  simdata <- list("mutlist"= mutlist, "pheno" = phenotype, "foldlist"=fold,
                  "annodata" = modelmatrix, "bmrpars" = bmrpars, "bmrmtxlist" = bmrmtxlist,
                  "para"=para, "efsize" = list("betaf0" = betaf0, "beta_gc" = betagc,
                  "avFe" = avFe, "diffFe" = diffFe, "diffFeFix" = diffFeFix),
                  "nsample"=c(Nsample.ps,Nsample.neu))
  return(simdata)
}

# -----------------------------------------------------------------------------
# Power-comparison drivers (independent case), ported from the diffdriver repo's
# scripts/simulate_functions.R but calling simulate_selection() (above) instead
# of simulate_1funcvi(). Each runs Niter simulations for one gene/setting and
# returns per-iteration p-values. They depend on the diffdriver package for the
# model/test functions: ddmodel(), ddmodel_binary_simple(), mlr(), genefisher(),
# genebinom(), genelr() -- load it (library(diffdriver) or devtools::load_all)
# before calling these.
# -----------------------------------------------------------------------------

#' DiffDriver power (fixed effect sizes). Uses the diffFeFix effect vector and
#' ddmodel / ddmodel_binary_simple. Saved by callers as `simuresdiffFix`.
power_comparediffiFix <- function(binary, Niter, sganno, sgmatrix, Nsample, para,
                                  bmrpars, betaf0, beta_gc, beta_gcFix, hot = 0, hmm) {
  m1.pvalue <- rep(1, Niter)
  a <- c(); b <- c(); d <- c(); f <- c()
  for (iter in 1:Niter) {
    print(iter)
    simdata <- simulate_selection(binary = binary, sganno = sganno, sgmatrix = sgmatrix,
      bmrpars = bmrpars, betaf0 = betaf0, Nsample = Nsample, beta_gc = beta_gc,
      beta_gcFix = beta_gcFix, para = para, hot = hot, hmm = hmm)
    mut <- do.call(rbind, simdata$mutlist)
    bmrmtx <- do.call(rbind, simdata$bmrmtxlist)
    ssgdata <- do.call(rbind, simdata$annodata)
    indexmtx <- cbind(bmrmtx[, 1], ssgdata)
    label <- factor(interaction(indexmtx))
    e <- simdata$pheno
    ef <- simdata$efsize
    fe <- ef$diffFeFix
    mr <- bmrmtx
    if (sum(mut) == 0) { next }
    if (binary == F) {
      res.m1 <- ddmodel(mut, e, mr, fe, label)
    } else {
      res.m1 <- ddmodel_binary_simple(mut, e, mr, fe)
    }
    m1.pvalue[iter] <- res.m1$pvalue
    parameters <- c(ef$beta_gc, ef$avbetaf1, ef$avbetaf2, ef$betaf1f2, ef$avbetaf1f2)
    a <- rbind(a, parameters)
    null <- c(res.m1$res.null$beta0, res.m1$res.null$alpha, res.m1$res.null$loglikelihood)
    alt  <- c(res.m1$res.alt$beta0,  res.m1$res.alt$alpha,  res.m1$res.alt$loglikelihood)
    d <- rbind(d, null)
    f <- rbind(f, alt)
    b <- c(b, sum(mut))
  }
  list(parameters = a, null = d, alt = f, m1.pvalue = m1.pvalue, "#mut" = b)
}

#' Benchmark-methods power: linear regression (mlr), Fisher, binomial, and
#' logistic tests. Saved by callers as `simuresother`.
power_compareotheri <- function(binary, Niter, sganno, sgmatrix, Nsample, para,
                                bmrpars, betaf0, beta_gc, hot = 0, hmm) {
  m1.pvalue <- m2.pvalue <- m3.pvalue <- m4.pvalue <- rep(1, Niter)
  a <- c(); b <- c()
  for (iter in 1:Niter) {
    print(iter)
    simdata <- simulate_selection(binary = binary, sganno = sganno, sgmatrix = sgmatrix,
      bmrpars = bmrpars, betaf0 = betaf0, Nsample = Nsample, beta_gc = beta_gc,
      para = para, hot = hot, hmm = hmm)
    mut <- do.call(rbind, simdata$mutlist)
    e <- simdata$pheno
    e_bisect <- ifelse(e > mean(e), 1, 0)
    ef <- simdata$efsize
    if (sum(mut) == 0) { next }
    res.m1 <- mlr(mut, e)
    res.m2 <- genefisher(mut, e_bisect)
    res.m3 <- genebinom(mut, e_bisect)
    res.m4 <- genelr(mut, e_bisect)
    m1.pvalue[iter] <- res.m1$pvalue
    m2.pvalue[iter] <- res.m2$pvalue
    m3.pvalue[iter] <- res.m3$pvalue
    m4.pvalue[iter] <- res.m4$pvalue
    parameters <- c(ef$beta_gc, ef$avbetaf1, ef$avbetaf2, ef$betaf1f2, ef$avbetaf1f2)
    a <- rbind(a, parameters)
    b <- c(b, sum(mut))
  }
  list(parameters = a, m1.pvalue = m1.pvalue, m2.pvalue = m2.pvalue,
       m3.pvalue = m3.pvalue, m4.pvalue = m4.pvalue, "#mut" = b)
}
