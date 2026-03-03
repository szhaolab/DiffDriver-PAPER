simulate_selection <- function(binary=F,sganno,sgmatrix, bmrpars, betaf0=2, Nsample, beta_gc, beta_gcFix=beta_gc, para,rho,tau,hot=0, hmm){
  if (binary==T){ # generate binary phenotype

    Nsamplec <- round(Nsample/2) # number of samples with phenotype E=1 (the rest will be 0)
    Nsamplen <- Nsample-Nsamplec
    phenotype <- c(rep(1,Nsamplec),rep(0,Nsamplen))
    ss=ifelse(phenotype==1,sample(c(0,1),size=Nsamplec,replace=T,prob = c(1-para[1],para[1])),sample(c(0,1),size=Nsamplen,replace=T,prob = c(1-para[2],para[2])))
    selection=rbind(ss,1-ss)
    Nsample.ps=sum(ss)
    Nsample.neu=Nsample-Nsample.ps
    #index=which(ss==1)
    #u=rnorm(n=length(phenotype),sd=sqrt((1/rho^2-1)*var(phenotype)))+phenotype
    #u1=u/sd(u)
    #u1=u1-min(u1)+1
    #bmrfold=u1*tau
    complement <- function(y, rho, x) {
      if (missing(x)) x <-runif(length(y), min=0, max =1) # Optional: supply a default if `x` is not given
      y.perp <- residuals(lm(x ~ y))
      rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
    }
    b=complement(phenotype,rho)
    b.new=b-min(b)+0.1
    bmrfold= (b.new/mean(b.new))*tau
    #bmrfold=rep(1,length(b.new))*exp(tau)
  }else{
    phenotype=rnorm(Nsample,mean = para[1],sd=para[2])
    pp=exp(para[3]+para[4]*phenotype)/(1+exp(para[3]+para[4]*phenotype))
    bpp=exp(bpara[1]+bpara[2]*phenotype)/(1+exp(bpara[1]+bpara[2]*phenotype))
    ss=ifelse(runif(Nsample)-pp<0,1,0) # generate positive samples
    bb=ifelse(runif(Nsample)-bpp<0,1,0) # generate pseudo phenotype samples
    Nsamplec <- sum(bb)
    Nsamplen <- Nsample-Nsamplec
    pseudophenotype=c(rep(1,Nsamplec),rep(0,Nsamplen))
    index=which(ss==1)
    phenotype=c(phenotype[index],phenotype[-index])# the positive samples are placed in the front.
    Nsample.ps1=sum(ss[1:Nsamplec]) # number of positive samples
    Nsample.ps0=sum(ss[(Nsamplec+1):Nsample]) # number of positive samples
    Nsample.ps=Nsample.ps0+Nsample.ps1
    Nsample.neu1 <- Nsamplec - Nsample.ps1 # number of neutral samples
    Nsample.neu0 <- Nsamplen - Nsample.ps0 # number of neutral samples
    Nsample.neu = Nsample.neu0+ Nsample.neu1
    pseudophenotype=c(pseudophenotype[index],pseudophenotype[-index])# the positive samples are placed in the front.
  }
  hotseq=hotspotseq(hmm,sganno) # The first column 'start' is the positions, and the second column 'seqt',is the hotspot indicator.
  if (hot==0){
    hotseq[,2]=0
  }

  mutlist <- list() # a list of nine mutation matrices
  countlist <- list()
  annodata <- list() # a list of nine annotation data frames
  bmrmtxlist <- list() # a list of nine background mutation matrices
  betagc=c(beta_gc,hmm[9])
  betagcFix=c(beta_gcFix,hmm[9])
  mutRate <- list()
  foldlist <- list()
  foldlistFix <- list()
  for (t in 1:length(sganno)) {
    hotseqt=merge(sganno[[t]],hotseq,by="start")$seqt
    selename=names(beta_gc)
    ssgdata=cbind(sgmatrix[[t]][,..selename],hotseqt)
    hotindex=which(hotseqt==1)
    pp.neu=matrix(rep(1,nrow(ssgdata)),ncol=1)%x%matrix(exp(bmrpars[t])*exp(betaf0)*bmrfold,nrow=1)
    fold=as.vector(exp(as.matrix(ssgdata)%*%betagc))
    fold[hotindex]=exp(hmm[9])
    fold==Nsample/Nsample.ps*fold-Nsample.ps/Nsample.neu
    if (any(fold<0)){warning("Inappropriate parameter settings!")}
    index=which(fold<=0)
    fold[index]=1
    F=cbind(fold,1)%*%selection
    pp.total=ifelse(F*pp.neu<1,F*pp.neu,0.99)
    foldlist[[t]]=data.table::data.table(fold=fold)

    foldFix=as.vector(exp(as.matrix(ssgdata)%*%betagcFix))
    foldFix[hotindex]=exp(hmm[9])
    foldFix=Nsample/Nsample.ps*foldFix-Nsample.ps/Nsample.neu
    index=which(foldFix<=0)
    foldFix[index]=1
    if (any(foldFix<0)){warning("Inappropriate parameter settings!")}
    foldlistFix[[t]]=data.table::data.table(fold=foldFix)


    if (nrow(pp.total)>1) {
      mutlist[[t]]=as(apply(pp.total,2,rbinom,n=nrow(pp.total),size=1),"sparseMatrix")
    }else{
      mutlist[[t]]=as(t(rbinom(length(pp.total),size=1,pp.total)),"sparseMatrix")
    }

    # if (!is.matrix(mutps)){
    # mutlist[[t]]=as(cbind(t(mutps),t(mutneu)),"sparseMatrix")
    # }else{
    # mutlist[[t]]=as(cbind(mutps,mutneu),"sparseMatrix")
    # }
    countlist[[t]] <- c(sum(mutlist[[t]]))
    bmrmtxlist[[t]] <- log(pp.neu)
  }
  # The forllowings are the ture parameters (???)
  fold <- do.call(rbind,foldlist)
  foldFix <- do.call(rbind,foldlistFix)
  avFe <- rep(log(mean(fold[[1]])*Nsample.ps/Nsample + Nsample.neu/Nsample),nrow(fold))
  diffFe <-  log(fold[[1]]*Nsample.ps/Nsample + Nsample.neu/Nsample)
  diffFeFix <-  log(foldFix[[1]]*Nsample.ps/Nsample + Nsample.neu/Nsample)

  simdata <- list("mutlist"= mutlist, "pheno" = phenotype,"foldlist"=fold,"bmrfold"=bmrfold, "annodata" = sganno, "bmrpars" = bmrpars, "bmrmtxlist" = bmrmtxlist, "para"=para, "efsize" = list( "betaf0" = betaf0,  "beta_gc" = betagc, "avFe" = avFe, "diffFe" = diffFe,"diffFeFix"=diffFeFix),"nsample"=c(Nsample.ps,Nsample.neu))
  return(simdata)
}
