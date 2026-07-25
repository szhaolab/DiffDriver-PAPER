# =============================================================================
# Simulation functions for the DiffDriver paper.
#
# Two simulation scenarios live here, each with a cohort simulator plus
# "benchmark runner" functions that repeat the simulation Niter times and
# return per-iteration p-values:
#
#   1. SELECTION  (analysis/simulation-power.Rmd) -- the phenotype drives
#      positive SELECTION on the gene (true differential signal). Used to
#      measure statistical POWER.
#        simulate_hotspot_seq()      hotspot indicator sequence via HMM
#        simulate_selection()        simulate one cohort under selection
#        run_selection_benchmarks()  Linear / Fisher / Binomial / LogisticR
#        run_selection_diffdriver()  DiffDriver (fixed effect)
#
#   2. CONFOUNDING (analysis/Simulation-BMR.Rmd) -- the phenotype correlates
#      with a mutational SIGNATURE (background rate), NOT with selection, so
#      there is no true signal. Any "hit" is a FALSE POSITIVE. Used to show
#      that count-based benchmarks are fooled while DiffDriver is robust.
#        build_confounded_bmr()             signature-confounded background rate
#        simulate_confounding()             simulate one cohort under confounding
#        run_confounding_benchmarks()       Linear / Fisher / Binomial /
#                                           LogisticR / Mann-Whitney
#        run_confounding_diffdriver()       DiffDriver (variable effect)
#        run_confounding_diffdriver_fixed() DiffDriver (fixed effect)
#
# The runners depend on the diffdriver package for the model/test functions
# (ddmodel, ddmodel_binary, ddmodel_binary_simple, mlr, genefisher, genebinom,
# genelr) -- load it (library(diffdriver) or devtools::load_all) first. The
# confounding scenario additionally needs the globals `codeSignature` (96
# trinucleotide-context -> number map) and `hotspot1sig` (list of 96 per-segment
# 0/1 hotspot indicators).
# =============================================================================


# =============================================================================
# 1. SELECTION scenario
# =============================================================================

#' Generate a hotspot indicator sequence via HMM
#'
#' @param hmm Numeric vector of 10 HMM parameters (transition probs and coefficients)
#' @param sgdata List of annotation data.frames (one per nucleotide type)
#' @return A matrix with columns: start (genomic position), seqt (hotspot indicator 0/1)
simulate_hotspot_seq <- function(hmm, sgdata){
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
#' Ported from szhaolab/diffdriver scripts/simulate_functions.R::simulate_1funcvi.
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
  hotseq=simulate_hotspot_seq(hmm,sganno)
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

#' Benchmark methods (Linear / Fisher / Binomial / LogisticR) under selection
#'
#' Runs simulate_selection() over Niter iterations and applies the four
#' count-based benchmark tests. Saved by callers as `simuresother`.
#'
#' @return list(parameters, m1.pvalue [Linear], m2.pvalue [Fisher],
#'   m3.pvalue [Binomial], m4.pvalue [LogisticR], "#mut")
run_selection_benchmarks <- function(binary, Niter, sganno, sgmatrix, Nsample, para,
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
    res.m1 <- mlr(mut, e)              # Linear regression
    res.m2 <- genefisher(mut, e_bisect)  # Fisher's exact test
    res.m3 <- genebinom(mut, e_bisect)   # Binomial test
    res.m4 <- genelr(mut, e_bisect)      # Logistic regression
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

#' DiffDriver under selection (fixed effect sizes)
#'
#' Runs simulate_selection() over Niter iterations and applies DiffDriver using
#' the diffFeFix effect vector. Saved by callers as `simuresdiffFix`.
#'
#' @return list(parameters, null, alt, m1.pvalue, "#mut")
run_selection_diffdriver <- function(binary, Niter, sganno, sgmatrix, Nsample, para,
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


# =============================================================================
# 2. CONFOUNDING scenario
#
# The phenotype is correlated with a mutational signature (default signature 3,
# a T>G / smoking-like signature in LUAD) rather than with selection. Mutations
# are simulated over the 96 trinucleotide-context segments of one gene, with two
# hotspots placed at the context with the highest loading in the confounding
# signature (LUAD signature 3 = nttype 34, C[T>G]G) so the hotspot-dominated
# mutations track the phenotype.
#
# NOTE: the fold guard in simulate_confounding() is a clamp (fold kept >= 0.5
# before the 2*fold-1 transform) rather than the original stop(); the original
# scripts/simulate_functions.R is lost, so results are close but not
# bit-identical to the archived runs.
# =============================================================================

#' Build a 96 x Nsample background-mutation-rate matrix confounded with phenotype
#'
#' The per-sample signature loadings are bootstrapped from the REAL loadings
#' (e.g. ll$LUAD, 6 signatures x samples), which are a composition (constant
#' column sums). The columns are then permuted so that the confounding
#' signature's loading tracks the phenotype `e` at correlation `rho`
#' (R^2 = rho^2). Because the summed loadings stay constant per sample, a
#' colSums covariate cannot adjust the confounding away.
#'
#' @param e Phenotype vector (length = number of samples)
#' @param signatures 96 x nsig matrix of signature profiles (columns = signatures)
#' @param loadings nsig x nsample matrix of real per-sample signature loadings
#' @param rho Target correlation between phenotype and the confounding loading
#' @param sc Scaling constant (larger => fewer simulated mutations)
#' @param sig_confound Index of the confounding signature (default 3)
#' @param tmb_cv Per-sample total-mutation-burden variation. 0 (default) keeps the
#'   constant per-sample budget of the LUAD loadings; tmb_cv > 0 multiplies each
#'   sample's whole loading vector by an independent mean-1 log-normal factor with
#'   coefficient of variation tmb_cv, so samples get different total mutation
#'   counts while the (confounded) signature *composition* is preserved.
#' @return list(bmr = 96 x nsample rate matrix, loadings = permuted loadings, signatures)
build_confounded_bmr <- function(e, signatures, loadings, rho, sc = 1, sig_confound = 3, tmb_cv = 0) {
  if (nrow(signatures) != 96) stop("Signature length should be 96!")
  nn <- length(e)
  cc <- as.matrix(loadings)[, sample(ncol(loadings), nn, replace = TRUE)]  # nsig x nn, colSums constant
  complement <- function(y, rho, x) {
    if (missing(x)) x <- runif(length(y), 0, 1)
    y.perp <- residuals(lm(x ~ y))
    rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
  }
  target <- complement(e, rho)                    # correlated with e at rho
  perm <- integer(nn); perm[order(target)] <- order(cc[sig_confound, ])
  cc <- cc[, perm]                                # confounding signature loading now tracks e
  if (tmb_cv > 0) {                               # per-sample total-burden multiplier (mean 1)
    sdlog <- sqrt(log(1 + tmb_cv^2))
    cc <- sweep(cc, 2, rlnorm(nn, meanlog = -sdlog^2 / 2, sdlog = sdlog), "*")
  }
  yy <- as.matrix(signatures) %*% cc
  bmr <- yy / sc
  list(bmr = bmr, loadings = cc, signatures = signatures)
}

#' Simulate mutations with the background rate confounded by a signature
#'
#' Ported from szhaolab/diffdriver scripts/simulate_functions.R::simulate_1funcv96.
#' The neutral per-site probability of nucleotide type t is bmrfold$bmr[t, ]
#' (from build_confounded_bmr), so it varies per sample and is correlated with
#' the phenotype. Selection is drawn from `para`; for the confounding benchmark
#' para = c(0.5, 0.5) gives NO differential selection.
#'
#' @param binary Logical; TRUE = binary case/control phenotype, FALSE = continuous
#'   phenotype ~ N(0,1) (both confounded with the signature loading at `rho`)
#' @param tmb_cv Per-sample total-mutation-burden variation (passed to
#'   build_confounded_bmr); 0 = constant burden
#' @param sganno List of 96 annotation data.frames (one per nucleotide type)
#' @param sgmatrix List of 96 functional annotation matrices
#' @param Nsample Total number of samples
#' @param beta_gc Named vector of functional effect coefficients
#' @param beta_gcFix Named vector of fixed functional effect coefficients (defaults to beta_gc)
#' @param para c(selection_prob_case, selection_prob_control)
#' @param hot Hotspot flag; 1 = apply hotspot effect from hotspot1sig
#' @param hmm HMM parameters (hmm[9] is the log hotspot fold)
#' @param signatures 96 x 7 signature data.frame (type + k1..k6)
#' @param loadings Real per-sample signature loadings (nsig x nsample)
#' @param rho Phenotype-signature correlation (R^2 = rho^2)
#' @param sc Background-rate scaling constant
#' @return list with mutlist, pheno, foldlist, covariate, bmrfold, annodata,
#'   bmrmtxlist, para, efsize, nsample
simulate_confounding <- function(binary = TRUE, sganno, sgmatrix, Nsample, beta_gc,
                                 beta_gcFix = beta_gc, para, hot = 0, hmm,
                                 signatures, loadings, rho, sc, tmb_cv = 0) {

  # --- phenotype + (null) selection ---------------------------------------
  if (isTRUE(binary)) {
    Nsamplec  <- round(Nsample / 2)        # samples with phenotype E = 1
    Nsamplen  <- Nsample - Nsamplec        # samples with phenotype E = 0
    phenotype <- c(rep(1, Nsamplec), rep(0, Nsamplen))
    ss <- ifelse(phenotype == 1,
                 sample(c(0, 1), size = Nsamplec, replace = TRUE, prob = c(1 - para[1], para[1])),
                 sample(c(0, 1), size = Nsamplen, replace = TRUE, prob = c(1 - para[2], para[2])))
  } else {
    # continuous phenotype ~ N(0,1); selection is independent of phenotype
    # (para[1] is the per-sample selection probability; para = c(0.5,0.5) => none)
    phenotype <- rnorm(Nsample)
    ss <- rbinom(Nsample, size = 1, prob = para[1])
  }
  selection   <- rbind(ss, 1 - ss)
  Nsample.ps  <- sum(ss)
  Nsample.neu <- Nsample - Nsample.ps

  # --- phenotype-confounded background rate over the 96 contexts ----------
  sig1 <- merge(signatures, codeSignature)
  sig2 <- sig1[order(sig1$number), c(2, 3, 4, 5, 6, 7)]   # k1..k6 in context order
  bmrfold <- build_confounded_bmr(e = phenotype, signatures = sig2,
                                  loadings = loadings, rho = rho, sc = sc, tmb_cv = tmb_cv)

  # --- simulate mutations, one matrix per nucleotide-type segment ---------
  betagc    <- c(beta_gc,    hmm[9])     # hmm[9] = log hotspot fold
  betagcFix <- c(beta_gcFix, hmm[9])
  mutlist     <- list()
  bmrmtxlist  <- list()
  foldlist    <- list()
  foldlistFix <- list()

  for (t in seq_along(sganno)) {
    if (nrow(sganno[[t]]) == 0) next
    hotseqt  <- hotspot1sig[[t]] * hot
    selename <- names(beta_gc)
    ssgdata  <- cbind(sgmatrix[[t]][, ..selename], hotseqt)
    hotindex <- which(hotseqt == 1)

    # neutral per-site rate: same across sites of type t, varies by sample
    pp.neu <- matrix(1, nrow = nrow(ssgdata), ncol = 1) %x%
              matrix(bmrfold$bmr[t, ], nrow = 1)

    # differential fold from the functional annotations (variable + fixed effect)
    fold <- as.vector(exp(as.matrix(ssgdata) %*% betagc))
    fold[hotindex] <- exp(hmm[9])
    fold[2 * fold < 1] <- 0.5            # clamp (orig. stop): keep 2*fold-1 >= 0
    fold <- 2 * fold - 1

    foldFix <- as.vector(exp(as.matrix(ssgdata) %*% betagcFix))
    foldFix[hotindex] <- exp(hmm[9])
    foldFix <- 2 * foldFix - 1

    foldlist[[t]]    <- data.table::data.table(fold = fold)
    foldlistFix[[t]] <- data.table::data.table(foldFix = foldFix)

    # per-sample mutation probability, then Bernoulli draws
    F        <- cbind(fold, 1) %*% selection
    pp.total <- ifelse(F * pp.neu < 1, F * pp.neu, 0.99)
    if (nrow(pp.total) > 1) {
      mutlist[[t]] <- as(apply(pp.total, 2, rbinom, n = nrow(pp.total), size = 1), "sparseMatrix")
    } else {
      mutlist[[t]] <- as(t(rbinom(length(pp.total), size = 1, pp.total)), "sparseMatrix")
    }
    bmrmtxlist[[t]] <- log(pp.neu)
  }

  # --- effect-size vectors for DiffDriver ---------------------------------
  fold    <- do.call(rbind, foldlist)
  foldFix <- do.call(rbind, foldlistFix)
  diffFe    <- log(fold[[1]]    * Nsample.ps / Nsample + Nsample.neu / Nsample)
  diffFeFix <- log(foldFix[[1]] * Nsample.ps / Nsample + Nsample.neu / Nsample)
  avFe      <- rep(mean(diffFe), nrow(fold))
  covariate <- apply(bmrfold$loadings, 2, sum)   # constant across samples (composition)

  list(mutlist = mutlist, pheno = phenotype, foldlist = fold, covariate = covariate,
       bmrfold = bmrfold, annodata = sganno, bmrmtxlist = bmrmtxlist, para = para,
       efsize = list(avFe = avFe, diffFe = diffFe, diffFeFix = diffFeFix),
       nsample = c(Nsample.ps, Nsample.neu))
}

#' Benchmark methods under signature confounding
#'
#' Runs simulate_confounding() over Niter iterations and applies the count-based
#' benchmarks: Linear (mlr), Fisher, Binomial, LogisticR (with the signature
#' covariate and without), and the Mann-Whitney rank-sum test on per-sample
#' mutation counts. With para = c(0.5, 0.5) there is no true signal, so all
#' rejections are false positives.
#'
#' @return list(m1.pvalue [Linear], m2.pvalue [Fisher], m3.pvalue [Binomial],
#'   m4.pvalue [LogisticR+cov], m5.pvalue [LogisticR no-cov],
#'   mannwhitney.pvalue, "#mut")
run_confounding_benchmarks <- function(binary, Niter, sganno, sgmatrix, Nsample,
                                       para, rho, signatures, loadings, beta_gc,
                                       hot = 0, hmm, sc) {
  m1.pvalue <- m2.pvalue <- m3.pvalue <- m4.pvalue <- m5.pvalue <-
    mannwhitney.pvalue <- rep(1, Niter)
  nmut <- c()
  for (iter in 1:Niter) {
    print(iter)
    sim <- simulate_confounding(binary = binary, sganno = sganno, sgmatrix = sgmatrix,
      Nsample = Nsample, beta_gc = beta_gc, para = para, signatures = signatures,
      loadings = loadings, rho = rho, hot = hot, hmm = hmm, sc = sc)
    mut       <- do.call(rbind, sim$mutlist)
    covariate <- sim$covariate
    e         <- sim$pheno
    e_bisect  <- ifelse(e > mean(e), 1, 0)
    if (sum(mut) == 0) next
    m1.pvalue[iter] <- mlr(mut, e)$pvalue                                         # Linear
    m2.pvalue[iter] <- genefisher(mut, e_bisect)$pvalue                           # Fisher
    m3.pvalue[iter] <- genebinom(mut, e_bisect)$pvalue                            # Binomial
    m4.pvalue[iter] <- genelr(mut, e_bisect, covariates = covariate)$pvalue       # LogisticR (+cov)
    m5.pvalue[iter] <- genelr(mut, e_bisect, covariates = rep(1, length(covariate)))$pvalue  # LogisticR (no cov)
    cnt <- as.numeric(Matrix::colSums(mut))                                       # Mann-Whitney
    mannwhitney.pvalue[iter] <- suppressWarnings(wilcox.test(cnt[e_bisect == 1], cnt[e_bisect == 0])$p.value)
    nmut <- c(nmut, sum(mut))
  }
  list(m1.pvalue = m1.pvalue, m2.pvalue = m2.pvalue, m3.pvalue = m3.pvalue,
       m4.pvalue = m4.pvalue, m5.pvalue = m5.pvalue,
       mannwhitney.pvalue = mannwhitney.pvalue, "#mut" = nmut)
}

#' DiffDriver under signature confounding (variable functional effect)
#'
#' Runs simulate_confounding() over Niter iterations and applies DiffDriver with
#' the diffFe effect vector. Saved by callers as `simuresdiff`.
#'
#' @return list(m1.pvalue, estimates, real, "#mut")
run_confounding_diffdriver <- function(binary, Niter, sganno, sgmatrix, Nsample,
                                       para, signatures, loadings, rho, beta_gc,
                                       hot = 0, hmm, sc) {
  m1.pvalue <- rep(1, Niter)
  nmut <- c(); estimates <- c()
  for (iter in 1:Niter) {
    print(iter)
    sim <- simulate_confounding(binary = binary, sganno = sganno, sgmatrix = sgmatrix,
      Nsample = Nsample, beta_gc = beta_gc, para = para, signatures = signatures,
      loadings = loadings, rho = rho, hot = hot, hmm = hmm, sc = sc)
    mut <- do.call(rbind, sim$mutlist)
    mr  <- do.call(rbind, sim$bmrmtxlist)
    e   <- sim$pheno
    fe  <- sim$efsize$diffFe
    if (sum(mut) == 0) next
    res <- if (binary == F) ddmodel(mut, e, mr, fe) else ddmodel_binary(mut, e, mr, fe)
    m1.pvalue[iter] <- res$pvalue
    nmut <- c(nmut, sum(mut))
    estimates <- rbind(estimates, res$all)
  }
  list(m1.pvalue = m1.pvalue, estimates = estimates, real = para, "#mut" = nmut)
}

#' DiffDriver under signature confounding (fixed functional effect)
#'
#' As run_confounding_diffdriver() but uses the diffFeFix effect vector (fixed
#' effect size). Saved by callers as `simuresdiffFix`.
#'
#' @return list(m1.pvalue, estimates, real, "#mut")
run_confounding_diffdriver_fixed <- function(binary, Niter, sganno, sgmatrix, Nsample,
                                             para, signatures, loadings, rho, beta_gc,
                                             beta_gcFix, hot = 0, hmm, sc) {
  m1.pvalue <- rep(1, Niter)
  nmut <- c(); estimates <- c()
  for (iter in 1:Niter) {
    print(iter)
    sim <- simulate_confounding(binary = binary, sganno = sganno, sgmatrix = sgmatrix,
      Nsample = Nsample, beta_gc = beta_gc, beta_gcFix = beta_gcFix, para = para,
      signatures = signatures, loadings = loadings, rho = rho, hot = hot, hmm = hmm, sc = sc)
    mut <- do.call(rbind, sim$mutlist)
    mr  <- do.call(rbind, sim$bmrmtxlist)
    e   <- sim$pheno
    fe  <- sim$efsize$diffFeFix
    if (sum(mut) == 0) next
    res <- if (binary == F) ddmodel(mut, e, mr, fe) else ddmodel_binary(mut, e, mr, fe)
    m1.pvalue[iter] <- res$pvalue
    nmut <- c(nmut, sum(mut))
    estimates <- rbind(estimates, res$all)
  }
  list(m1.pvalue = m1.pvalue, estimates = estimates, real = para, "#mut" = nmut)
}
