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

#' Simulate mutations under a phenotype-dependent selection model
#'
#' @param binary Logical; TRUE for binary (case/control), FALSE for continuous phenotype
#' @param sganno List of 9 annotation data.frames (one per nucleotide type)
#' @param sgmatrix List of 9 functional annotation matrices
#' @param bmrpars Log-scale background mutation rates (length 9)
#' @param betaf0 Baseline functional effect size (log-scale)
#' @param Nsample Total number of samples
#' @param beta_gc Named vector of functional effect coefficients
#' @param para For binary: c(selection_prob_case, selection_prob_control).
#'   For continuous: c(mean, sd, intercept, slope) where intercept and slope
#'   are logistic regression coefficients for selection probability given phenotype.
#' @param rho Correlation between phenotype and BMR fold change
#' @param tau BMR fold-change scaling factor
#' @param hot Hotspot flag; 0 = no hotspot effect
#' @param hmm HMM parameters for hotspot sequence generation
#' @return A list containing mutlist, pheno, foldlist, bmrfold, annodata,
#'   bmrpars, bmrmtxlist, para, efsize, nsample
simulate_selection <- function(binary=F, sganno, sgmatrix, bmrpars, betaf0=2, Nsample, beta_gc, para, rho, tau, hot=0, hmm){
  if (binary==T){ # generate binary phenotype
    Nsamplec <- round(Nsample/2)
    Nsamplen <- Nsample-Nsamplec
    phenotype <- c(rep(1,Nsamplec),rep(0,Nsamplen))
    ss=ifelse(phenotype==1,sample(c(0,1),size=Nsamplec,replace=T,prob = c(1-para[1],para[1])),sample(c(0,1),size=Nsamplen,replace=T,prob = c(1-para[2],para[2])))
    selection=rbind(ss,1-ss)
    Nsample.ps=sum(ss)
    Nsample.neu=Nsample-Nsample.ps

    complement <- function(y, rho, x) {
      if (missing(x)) x <- runif(length(y), min=0, max=1)
      y.perp <- residuals(lm(x ~ y))
      rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
    }
    b=complement(phenotype,rho)
    b.new=b-min(b)+0.1
    bmrfold= (b.new/mean(b.new))*tau
  }else{
    phenotype=rnorm(Nsample,mean = para[1],sd=para[2])
    pp=exp(para[3]+para[4]*phenotype)/(1+exp(para[3]+para[4]*phenotype))
    bpp=exp(bpara[1]+bpara[2]*phenotype)/(1+exp(bpara[1]+bpara[2]*phenotype))
    ss=ifelse(runif(Nsample)-pp<0,1,0)
    bb=ifelse(runif(Nsample)-bpp<0,1,0)
    Nsamplec <- sum(bb)
    Nsamplen <- Nsample-Nsamplec
    pseudophenotype=c(rep(1,Nsamplec),rep(0,Nsamplen))
    index=which(ss==1)
    phenotype=c(phenotype[index],phenotype[-index])
    Nsample.ps1=sum(ss[1:Nsamplec])
    Nsample.ps0=sum(ss[(Nsamplec+1):Nsample])
    Nsample.ps=Nsample.ps0+Nsample.ps1
    Nsample.neu1 <- Nsamplec - Nsample.ps1
    Nsample.neu0 <- Nsamplen - Nsample.ps0
    Nsample.neu = Nsample.neu0+ Nsample.neu1
    pseudophenotype=c(pseudophenotype[index],pseudophenotype[-index])
  }
  hotseq=hotspotseq(hmm,sganno)
  if (hot==0){
    hotseq[,2]=0
  }

  mutlist <- list()
  countlist <- list()
  bmrmtxlist <- list()
  betagc=c(beta_gc,hmm[9])
  foldlist <- list()
  for (t in 1:length(sganno)) {
    hotseqt=merge(sganno[[t]],hotseq,by="start")$seqt
    selename=names(beta_gc)
    ssgdata=cbind(sgmatrix[[t]][,..selename],hotseqt)
    hotindex=which(hotseqt==1)
    pp.neu=matrix(rep(1,nrow(ssgdata)),ncol=1)%x%matrix(exp(bmrpars[t])*exp(betaf0)*bmrfold,nrow=1)
    fold=as.vector(exp(as.matrix(ssgdata)%*%betagc))
    fold[hotindex]=exp(hmm[9])
    fold=Nsample/Nsample.ps*fold-Nsample.ps/Nsample.neu
    if (any(fold<0)){warning("Inappropriate parameter settings!")}
    index=which(fold<=0)
    fold[index]=1
    F=cbind(fold,1)%*%selection
    pp.total=ifelse(F*pp.neu<1,F*pp.neu,0.99)
    foldlist[[t]]=data.table::data.table(fold=fold)

    if (nrow(pp.total)>1) {
      mutlist[[t]]=as(apply(pp.total,2,rbinom,n=nrow(pp.total),size=1),"sparseMatrix")
    }else{
      mutlist[[t]]=as(t(rbinom(length(pp.total),size=1,pp.total)),"sparseMatrix")
    }

    countlist[[t]] <- c(sum(mutlist[[t]]))
    bmrmtxlist[[t]] <- log(pp.neu)
  }
  fold <- do.call(rbind,foldlist)
  diffFe <- log(fold[[1]]*Nsample.ps/Nsample + Nsample.neu/Nsample)

  simdata <- list("mutlist"= mutlist, "pheno" = phenotype, "foldlist"=fold, "bmrfold"=bmrfold, "annodata" = sganno, "bmrpars" = bmrpars, "bmrmtxlist" = bmrmtxlist, "para"=para, "efsize" = list("betaf0" = betaf0, "beta_gc" = betagc, "diffFe" = diffFe), "nsample"=c(Nsample.ps,Nsample.neu))
  return(simdata)
}
