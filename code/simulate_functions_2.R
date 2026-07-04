# Independent Case

simulate_1funcvi <- function(binary=F,sganno,sgmatrix, bmrpars, betaf0=2, Nsample, beta_gc,beta_gcFix=beta_gc, para,hot=0, hmm){
	if (binary==T){ # generate binary phenotype
		Nsamplec <- round(Nsample/2) # number of samples with phenotype E=1 (the rest will be 0)
		Nsamplen <- Nsample-Nsamplec
		phenotype <- c(rep(1,Nsamplec),rep(0,Nsamplen))
		ss=ifelse(phenotype==1,sample(c(0,1),size=Nsamplec,replace=T,prob = c(1-para[1],para[1])),sample(c(0,1),size=Nsamplen,replace=T,prob = c(1-para[2],para[2])))
		}else{
		phenotype=rnorm(Nsample,mean = para[1],sd=para[2])
		pp=exp(para[3]+para[4]*phenotype)/(1+exp(para[3]+para[4]*phenotype))
		ss=ifelse(runif(Nsample)-pp<0,1,0) # generate positive samples
		}
		index=which(ss==1)
		phenotype=c(phenotype[index],phenotype[-index])# the positive samples are placed in the front.
		Nsample.ps=sum(ss) # number of positive samples
		Nsample.neu <- Nsample - Nsample.ps # number of neutral samples
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
		foldlist[[t]]=data.table(fold=fold)
		
		foldFix=exp(as.matrix(ssgdata)%*%betagcFix)
		foldFix[hotindex]=exp(hmm[9])
		foldFix=(Nsample/Nsample.ps)*foldFix-Nsample.neu/Nsample.ps
		if (any(foldFix<0)){warning("Foldfix is negative!")}
		foldFix[which(foldFix<0)]=1
		foldlistFix[[t]]=data.table(fold=foldFix)
		

		mutps=replicate(Nsample.ps,rbinom(length(pp.ps),size=1,pp.ps))
		mutneu=replicate(Nsample.neu,rbinom(length(pp.neu),size=1,pp.neu))
		if (!is.matrix(mutps)){
		mutlist[[t]]=as(cbind(t(mutps),t(mutneu)),"sparseMatrix")
		}else{
		mutlist[[t]]=as(cbind(mutps,mutneu),"sparseMatrix")
		}
		countlist[[t]] <- c(sum(mutlist[[t]]))
		bmrmtxlist[[t]] <- matrix(bmrpars[t]+betaf0, ncol = ncol(mutlist[[t]]), nrow = nrow(mutlist[[t]])) # background mutation matrix for nytype=t
		modelmatrix[[t]] <- ssgdata
		
}

# The forllowings are the ture parameters (???)
	fold <- do.call(rbind,foldlist)
	foldFix <- do.call(rbind,foldlistFix)
	avFe <- rep(log(mean(fold[[1]])*Nsample.ps/Nsample + Nsample.neu/Nsample),nrow(fold))
	diffFe <-  log(fold[[1]]*Nsample.ps/Nsample + Nsample.neu/Nsample)
	diffFeFix <-  log(foldFix[[1]]*Nsample.ps/Nsample + Nsample.neu/Nsample)


	simdata <- list("mutlist"= mutlist, "pheno" = phenotype,"foldlist"=fold, "annodata" = modelmatrix, "bmrpars" = bmrpars, "bmrmtxlist" = bmrmtxlist, "para"=para, "efsize" = list( "betaf0" = betaf0,  "beta_gc" = betagc, "avFe" = avFe, "diffFe" = diffFe,"diffFeFix"=diffFeFix),"nsample"=c(Nsample.ps,Nsample.neu))
	return(simdata)
}





power_compareotheri <- function(binary, Niter, sganno,sgmatrix, Nsample,para,bmrpars,betaf0,beta_gc,hot=0,hmm){

  m1.pvalue <- m2.pvalue <- m3.pvalue <- m4.pvalue <- rep(1,Niter)
	a=c()
	b=c()

	for (iter in 1:Niter) {
	  print(iter)
		simdata <- simulate_1funcvi(binary=binary,sganno=sganno,sgmatrix=sgmatrix, bmrpars=bmrpars, betaf0=betaf0, Nsample=Nsample, beta_gc=beta_gc, para=para,hot=hot,hmm=hmm)
		ssgdata=simdata$annodata
		mut <- do.call(rbind, simdata$mutlist)
		bmrmtx <- do.call(rbind, simdata$bmrmtxlist)
		e <- simdata$pheno
		e_bisect=ifelse(e>mean(e),1,0)
		ef <- simdata$efsize
		print(sum(mut))
		if (sum(mut) ==0) {next}
		res.m1 <- mlr(mut,e)
		res.m2 <- genefisher(mut,e_bisect)
		res.m3 <- genebinom(mut,e_bisect)
		res.m4 <- genelr(mut,e_bisect)
		m1.pvalue[iter] <-  res.m1$pvalue
		m2.pvalue[iter] <-  res.m2$pvalue
		m3.pvalue[iter] <-  res.m3$pvalue
		m4.pvalue[iter] <-  res.m4$pvalue
		parameters=c(ef$beta_gc,ef$avbetaf1,ef$avbetaf2,ef$betaf1f2,ef$avbetaf1f2)
		a=rbind(a,parameters)
		nummut=sum(mut)
		b=c(b,nummut)
		}
		return(list("parameters"=a, "m1.pvalue" =m1.pvalue, "m2.pvalue" =m2.pvalue,
               "m3.pvalue" =m3.pvalue,"m4.pvalue" =m4.pvalue,"#mut"=b))
		}




power_comparebasei <- function(binary, Niter, sganno,sgmatrix, Nsample,para,bmrpars,betaf0,beta_gc,hot=0,hmm){
  m1.pvalue <-  rep(1,Niter)
	a=c()
	b=c()
	d=c()
	f=c()
	for (iter in 1:Niter) {
		print(iter)
		simdata <- simulate_1funcvi(binary=binary,sganno=sganno,sgmatrix=sgmatrix, bmrpars=bmrpars, betaf0=betaf0, Nsample=Nsample, beta_gc=beta_gc, para=para,hot=hot,hmm=hmm)
		
		ssgdata=do.call(rbind,simdata$annodata)
		mut <- do.call(rbind, simdata$mutlist)
		bmrmtx <- do.call(rbind, simdata$bmrmtxlist)
		e <- simdata$pheno
		ef <- simdata$efsize
		fe <- ef$avFe
		mr <- bmrmtx
		indexmtx=cbind(bmrmtx[,1],ssgdata)
		label=factor(interaction(indexmtx))	
		if (sum(mut) ==0) {
		next()
}
		if (binary==F){
		res.m1 <- ddmodel(mut,e, mr, fe, label=label)
		}else{
		res.m1 <- ddmodel_binary_simple(mut,e,mr,fe)
		}
		m1.pvalue[iter] <-  res.m1$pvalue
		parameters=c(ef$beta_gc,ef$avbetaf1,ef$avbetaf2,ef$betaf1f2,ef$avbetaf1f2)
		null=c(res.m1$res.null$beta0,res.m1$res.null$alpha,res.m1$res.null$loglikelihood)
		alt=c(res.m1$res.alt$beta0,res.m1$res.alt$alpha,res.m1$res.alt$loglikelihood)
		d=rbind(d,null)
		f=rbind(f,alt)
		a=rbind(a,parameters)
		nummut=sum(mut)
		b=c(b,nummut)
  }
	return(list("parameters"=a,"null"=d,"alt"=f, "m1.pvalue" =m1.pvalue,"#mut"=b))
	#return(list("mut"=mut,"e"=e,"mr"=mr,"fe"=fe,"output"=res.m1)) 
}




power_comparediffiFix <-function(binary, Niter, sganno,sgmatrix, Nsample,para,bmrpars,betaf0,beta_gc,beta_gcFix,hot=0,hmm){
  m1.pvalue <- rep(1,Niter)
	a=c()
	b=c()
	d=c()
	f=c()
	for (iter in 1:Niter) {
		print(iter)
		simdata <- simulate_1funcvi(binary=binary,sganno=sganno,sgmatrix=sgmatrix, bmrpars=bmrpars, betaf0=betaf0, Nsample=Nsample, beta_gc=beta_gc,beta_gcFix=beta_gcFix, para=para,hot=hot,hmm=hmm)
		mut <- do.call(rbind, simdata$mutlist)
		bmrmtx <- do.call(rbind, simdata$bmrmtxlist)
		ssgdata=do.call(rbind,simdata$annodata)
		indexmtx=cbind(bmrmtx[,1],ssgdata)
		label=factor(interaction(indexmtx))	
		e <- simdata$pheno
		ef <- simdata$efsize
		fe <- ef$diffFeFix
		mr <- bmrmtx
		if (sum(mut) ==0) {next}
		if (binary == F){
		res.m1 <- ddmodel(mut,e, mr, fe,label)
		}else{
		res.m1 <- ddmodel_binary_simple(mut,e,mr,fe)
		}
		m1.pvalue[iter] <-  res.m1$pvalue
		parameters=c(ef$beta_gc,ef$avbetaf1,ef$avbetaf2,ef$betaf1f2,ef$avbetaf1f2)
		a=rbind(a,parameters)
		null=c(res.m1$res.null$beta0,res.m1$res.null$alpha,res.m1$res.null$loglikelihood)
		alt=c(res.m1$res.alt$beta0,res.m1$res.alt$alpha,res.m1$res.alt$loglikelihood)
		d=rbind(d,null)
		f=rbind(f,alt)
		nummut=sum(mut)
		b=c(b,nummut)
  }
   return(list("parameters"=a,"null"=d,"alt"=f, "m1.pvalue" =m1.pvalue,"#mut"=b))
	#return(res.m1)  
}



power_comparediffi <-function(binary, Niter, sganno,sgmatrix, Nsample,para,bmrpars,betaf0,beta_gc,hot=0,hmm){
  m1.pvalue <- rep(1,Niter)
	a=c()
	b=c()
	d=c()
	f=c()
	for (iter in 1:Niter) {
		print(iter)
		simdata <- simulate_1funcvi(binary=binary,sgann=sganno,sgmatrix=sgmatrix, bmrpars=bmrpars, betaf0=betaf0, Nsample=Nsample, beta_gc=beta_gc, para=para,hot=hot,hmm=hmm)
		mut <- do.call(rbind, simdata$mutlist)
		bmrmtx <- do.call(rbind, simdata$bmrmtxlist)
		ssgdata=do.call(rbind,simdata$annodata)
		indexmtx=cbind(bmrmtx[,1],ssgdata)
		label=factor(interaction(indexmtx))	
		e <- simdata$pheno
		ef <- simdata$efsize
		fe <- ef$diffFe
		mr <- bmrmtx
		if (sum(mut) ==0) {next}
		if (binary == F){
		res.m1 <- ddmodel(mut,e, mr, fe,label)
		}else{
		res.m1 <- ddmodel_binary_simple(mut,e,mr,fe)
		}
		m1.pvalue[iter] <-  res.m1$pvalue
		parameters=c(ef$beta_gc,ef$avbetaf1,ef$avbetaf2,ef$betaf1f2,ef$avbetaf1f2)
		a=rbind(a,parameters)
		null=c(res.m1$res.null$beta0,res.m1$res.null$alpha,res.m1$res.null$loglikelihood)
		alt=c(res.m1$res.alt$beta0,res.m1$res.alt$alpha,res.m1$res.alt$loglikelihood)
		d=rbind(d,null)
		f=rbind(f,alt)
		nummut=sum(mut)
		b=c(b,nummut)
  }
   return(list("parameters"=a,"null"=d,"alt"=f, "m1.pvalue" =m1.pvalue,"#mut"=b))
	#return(res.m1)  
}



power_comparediffiold <-function(binary, Niter, sganno,sgmatrix, Nsample,para,bmrpars,betaf0,beta_gc,hot=0,hmm){
  m1.pvalue <- rep(1,Niter)
	a=c()
	b=c()
	d=c()
	f=c()
	for (iter in 1:Niter) {
		print(iter)
		simdata <- simulate_1funcvi(binary=binary,sganno,sgmatrix, bmrpars, betaf0, Nsample, beta_gc, para,hot,hmm)
		mut <- do.call(rbind, simdata$mutlist)
		bmrmtx <- do.call(rbind, simdata$bmrmtxlist)
		ssgdata=do.call(rbind,simdata$annodata)
		indexmtx=cbind(bmrmtx[,1],ssgdata)
		label=factor(interaction(indexmtx))	
		e <- simdata$pheno
		ef <- simdata$efsize
		fe <- ef$diffFe
		mr <- bmrmtx
		if (sum(mut) ==0) {next}
		if (binary == F){
		res.m1 <- ddmodelold(mut,e, mr, fe,label)
		}else{
		res.m1 <- ddmodel_binary_simple(mut,e,mr,fe)
		}
		m1.pvalue[iter] <-  res.m1$pvalue
		parameters=c(ef$beta_gc,ef$avbetaf1,ef$avbetaf2,ef$betaf1f2,ef$avbetaf1f2)
		a=rbind(a,parameters)
		null=c(res.m1$res.null$beta0,res.m1$res.null$alpha,res.m1$res.null$loglikelihood)
		alt=c(res.m1$res.alt$beta0,res.m1$res.alt$alpha,res.m1$res.alt$loglikelihood)
		d=rbind(d,null)
		f=rbind(f,alt)
		nummut=sum(mut)
		b=c(b,nummut)
  }
   return(list("parameters"=a,"null"=d,"alt"=f, "m1.pvalue" =m1.pvalue,"#mut"=b))
	#return(res.m1)  
}





# Correlation 9 Case

simulate_1funcv9 <- function(binary=F,sganno,sgmatrix, bmrpars, betaf0=2, Nsample, beta_gc,beta_gcFix=beta_gc, para,rho,tau,hot=0, hmm){
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





power_compareother9 <- function(binary, Niter, sganno,sgmatrix, Nsample,para,rho,tau,bmrpars,betaf0,beta_gc,hot=0,hmm){

  m1.pvalue <- m2.pvalue <- m3.pvalue <- m4.pvalue <- rep(1,Niter)
	a=c()
	b=c()

	for (iter in 1:Niter) {
print(iter)
		simdata <- simulate_1funcv9(binary=binary,sganno=sganno,sgmatrix=sgmatrix, bmrpars=bmrpars, betaf0=betaf0, Nsample=Nsample, beta_gc=beta_gc, para=para, rho=rho,tau=tau,hot=hot,hmm=hmm)
		ssgdata=simdata$annodata
		mut <- do.call(rbind, simdata$mutlist)
		bmrmtx <- do.call(rbind, simdata$bmrmtxlist)
		bmrfold <- simdata$bmrfold
		e <- simdata$pheno
		e_bisect=ifelse(e>mean(e),1,0)
		ef <- simdata$efsize
		if (sum(mut) ==0) {next}
		res.m1 <- mlr(mut=mut,e=e,covariates= bmrfold)
		res.m2 <- genefisher(mut=mut,e=e_bisect)
		res.m3 <- genebinom(mut=mut,e=e_bisect)
		res.m4 <- genelr(mut=mut,e=e_bisect,covariates= bmrfold)
		m1.pvalue[iter] <-  res.m1$pvalue
		m2.pvalue[iter] <-  res.m2$pvalue
		m3.pvalue[iter] <-  res.m3$pvalue
		m4.pvalue[iter] <-  res.m4$pvalue
		parameters=c(ef$beta_gc,ef$avbetaf1,ef$avbetaf2,ef$betaf1f2,ef$avbetaf1f2)
		a=rbind(a,parameters)
		nummut=sum(mut)
		b=c(b,nummut)
		}
		return(list("m1.pvalue" =m1.pvalue, "m2.pvalue" =m2.pvalue,
               "m3.pvalue" =m3.pvalue,"m4.pvalue" =m4.pvalue,"#mut"=b))
		}




power_comparebase9 <- function(binary, Niter, sganno,sgmatrix,Nsample,para,rho,tau=1,bmrpars,betaf0,beta_gc,hot=0,hmm){
  m1.pvalue <- rep(1,Niter)
	a=c()
	b=c()

	for (iter in 1:Niter) {
		simdata <- simulate_1funcv9(binary=binary,sganno,sgmatrix, bmrpars, betaf0, Nsample, beta_gc, para,rho,tau,hot,hmm)
		ssgdata=simdata$annodata
		mut <- do.call(rbind, simdata$mutlist)
		bmrmtx <- do.call(rbind, simdata$bmrmtxlist)
		e <- simdata$pheno
		ef <- simdata$efsize
		fe <- ef$avFe
		mr <- bmrmtx
		if (sum(mut) ==0) {next}

		if (binary == F){
		res.m1 <- ddmodel(mut,e,mr, fe)
		}else{
		res.m1 <- ddmodel_binary(mut,e,mr,fe)
		}
		m1.pvalue[iter] <-  res.m1$pvalue
		parameters=c(ef$beta_gc,ef$avbetaf1,ef$avbetaf2,ef$betaf1f2,ef$avbetaf1f2)
		a=rbind(a,parameters)
		nummut=sum(mut)
		b=c(b,nummut)
  }
   return(list("parameters"=a, "m1.pvalue" =m1.pvalue,"#mut"=b))
  }





power_comparediff9 <- function(binary, Niter, sganno,sgmatrix,Nsample,para,rho,tau=1,bmrpars,betaf0,beta_gc,hot=0,hmm){
  m1.pvalue <- rep(1,Niter)
	a=c()
	b=c()

	for (iter in 1:Niter) {
print(iter)
		simdata <- simulate_1funcv9(binary=binary,sganno=sganno,sgmatrix=sgmatrix, bmrpars=bmrpars, betaf0=betaf0, Nsample=Nsample, beta_gc=beta_gc, para=para,rho=rho,tau=tau,hot=hot,hmm=hmm)
		ssgdata=simdata$annodata
		mut <- do.call(rbind, simdata$mutlist)
		bmrmtx <- do.call(rbind, simdata$bmrmtxlist)
		e <- simdata$pheno
		ef <- simdata$efsize
		fe <- ef$diffFe
		mr <- bmrmtx
		if (sum(mut) ==0) {next}

		if (binary == F){
		res.m1 <- ddmodel(mut,e,mr, fe)
		}else{
		res.m1 <- ddmodel_binary(mut,e,mr,fe)
		}
		m1.pvalue[iter] <-  res.m1$pvalue
		parameters=c(ef$beta_gc,ef$avbetaf1,ef$avbetaf2,ef$betaf1f2,ef$avbetaf1f2)
		a=rbind(a,parameters)
		nummut=sum(mut)
		b=c(b,nummut)
  }
   return(list("parameters"=a, "m1.pvalue" =m1.pvalue,"#mut"=b))
  }




power_comparediff9Fix <- function(binary, Niter, sganno,sgmatrix,Nsample,para,rho,tau=1,bmrpars,betaf0,beta_gc,beta_gcFix=beta_gc,hot=0,hmm){
  m1.pvalue <- rep(1,Niter)
	a=c()
	b=c()

	for (iter in 1:Niter) {
		print(iter)
		simdata <- simulate_1funcv9(binary=binary,sganno=sganno,sgmatrix=sgmatrix, bmrpars=bmrpars, betaf0=betaf0, Nsample=Nsample, beta_gc=beta_gc,beta_gcFix=beta_gcFix, para=para,rho=rho,tau=tau,hot=hot,hmm=hmm)
		ssgdata=simdata$annodata
		mut <- do.call(rbind, simdata$mutlist)
		bmrmtx <- do.call(rbind, simdata$bmrmtxlist)
		e <- simdata$pheno
		ef <- simdata$efsize
		fe <- ef$diffFeFix
		mr <- bmrmtx
		if (sum(mut) ==0) {next}

		if (binary == F){
		res.m1 <- ddmodel(mut,e,mr, fe)
		}else{
		res.m1 <- ddmodel_binary(mut,e,mr,fe)
		}
		m1.pvalue[iter] <-  res.m1$pvalue
		parameters=c(ef$beta_gc,ef$avbetaf1,ef$avbetaf2,ef$betaf1f2,ef$avbetaf1f2)
		a=rbind(a,parameters)
		nummut=sum(mut)
		b=c(b,nummut)
  }
   return(list("parameters"=a, "m1.pvalue" =m1.pvalue,"#mut"=b))
  }



# Corrleation 96 Case

simulate_1funcv96 <- function(binary=F,sganno,sgmatrix, Nsample, beta_gc,beta_gcFix=beta_gc, para,hot=0, hmm=hmm,signatures,rho,sc){
	if (binary==T){ # generate binary phenotype
		Nsamplec <- round(Nsample/2) # number of samples with phenotype E=1 (the rest will be 0)
		Nsamplen <- Nsample-Nsamplec
		phenotype <- c(rep(1,Nsamplec),rep(0,Nsamplen))
		ss=ifelse(phenotype==1,sample(c(0,1),size=Nsamplec,replace=T,prob = c(1-para[1],para[1])),sample(c(0,1),size=Nsamplen,replace=T,prob = c(1-para[2],para[2])))
		selection=rbind(ss,1-ss)
		Nsample.ps=sum(ss)
   	Nsample.neu=Nsample-Nsample.ps
   	sig1=merge(signatures,codeSignature)
   	sig2=sig1[order(sig1$number),c(3,2,4,5)]
		bmrfold= bmrSignature(e=phenotype,signatures=sig2,rho=rho,sc=sc)# generate the 96 times n bmr matrix
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
		# hotseq=hotspotseq(hmm,sganno) # The first column 'start' is the positions, and the second column 'seqt',is the hotspot indicator.
		# if (hot==0){
		# 	hotseq[,2]=0
		# }

	mutlist <- list() # a list of nine mutation matrices
	countlist <- list()
	annodata <- list() # a list of nine annotation data frames
	bmrmtxlist <- list() # a list of nine background mutation matrices
	betagc=c(beta_gc,hmm[9])
	betagcFix=c(beta_gcFix,hmm[9])
	mutRate <- list()
	foldlist <- list()
	foldlistFix <- list()
	size=c()
	for (t in 1:length(sganno)) {
	  	size=c(size,nrow(sganno[[t]]))
	  	if (nrow(sganno[[t]])==0) {next}
	  	hotseqt= hotspot1sig[[t]]*hot
	  	selename=names(beta_gc)
		ssgdata=cbind(sgmatrix[[t]][,..selename],hotseqt)
		hotindex=which(hotseqt==1)
		pp.neu=matrix(rep(1,nrow(ssgdata)),ncol=1)%x%matrix(bmrfold$bmr[t,],nrow=1)
		fold=as.vector(exp(as.matrix(ssgdata)%*%betagc))
		fold[hotindex]=exp(hmm[9])
		if (any(2*fold<1)){stop("Error:inappropriate parameter settings!")}
		fold=2*fold-1

		foldFix=as.vector(exp(as.matrix(ssgdata)%*%betagcFix))
		foldFix[hotindex]=exp(hmm[9])
		#if (any(2*foldFix<1)){stop("Error:inappropriate parameter settings!")}
		foldFix=2*foldFix-1






		F=cbind(fold,1)%*%selection
		pp.total=ifelse(F*pp.neu<1,F*pp.neu,0.99)
		foldlist[[t]]=data.table::data.table(fold=fold)
		foldlistFix[[t]]=data.table::data.table(foldFix=foldFix)
		if (nrow(pp.total)>1) {
		mutlist[[t]]=as(apply(pp.total,2,rbinom,n=nrow(pp.total),size=1),"sparseMatrix")
		}else{
		  mutlist[[t]]=as(t(rbinom(length(pp.total),size=1,pp.total)),"sparseMatrix")
		  }
		countlist[[t]] <- c(sum(mutlist[[t]]))
		bmrmtxlist[[t]] <- log(pp.neu)
}
# The forllowings are the ture parameters (???)

	fold <- do.call(rbind,foldlist)
	foldFix <- do.call(rbind,foldlistFix)
	#avFe <- rep(log(mean(fold[[1]])*Nsample.ps/Nsample + Nsample.neu/Nsample),nrow(fold))
	#avFe <- rep(log(mean(fold[[1]])*Nsample.ps/Nsample + Nsample.neu/Nsample),nrow(fold))
	diffFe <-  log(fold[[1]]*Nsample.ps/Nsample + Nsample.neu/Nsample)
 	diffFeFix <-  log(foldFix[[1]]*Nsample.ps/Nsample + Nsample.neu/Nsample)
 
	avFe <-rep(mean(diffFe),nrow(fold))
	#covariate=apply(bmrfold$bmr, 2, "%*%",size)
	covariate=apply(bmrfold$loadings, 2, sum)
	simdata <- list("mutlist"= mutlist, "pheno" = phenotype,"foldlist"=fold,"covariate"=covariate,"bmrfold"=bmrfold, "annodata" = sganno, "bmrmtxlist" = bmrmtxlist, "para"=para, "efsize" = list(  "avFe" = avFe, "diffFe" = diffFe,"diffFeFix"=diffFeFix),"nsample"=c(Nsample.ps,Nsample.neu))
	return(simdata)
}





power_compareother96 <- function(binary, Niter, sganno,sgmatrix, Nsample,para,rho,signatures,beta_gc,hot=0,hmm,sc){

m1.pvalue <- m2.pvalue <- m3.pvalue <- m4.pvalue<- m5.pvalue <- rep(1,Niter)
	a=c()
	b=c()

	for (iter in 1:Niter) {
	  print(iter)
		simdata <- simulate_1funcv96(binary=binary,sganno=sganno,sgmatrix=sgmatrix, signatures=signatures,Nsample=Nsample, beta_gc=beta_gc, para=para, rho=rho,hot=hot,hmm=hmm,sc=sc)
		ssgdata=simdata$annodata
		mut <- do.call(rbind, simdata$mutlist)
		bmrmtx <- do.call(rbind, simdata$bmrmtxlist)
		covariate <- simdata$covariate
		e <- simdata$pheno
		e_bisect=ifelse(e>mean(e),1,0)
		ef <- simdata$efsize
		if (sum(mut) ==0) {next}
		res.m1 <- mlr(mut,e)
		res.m2 <- genefisher(mut,e_bisect)
		res.m3 <- genebinom(mut,e_bisect)
		res.m4 <- genelr(mut,e_bisect,covariates= covariate)
		res.m5 <- genelr(mut,e_bisect,covariates=rep(1,length(covariate)))
		m1.pvalue[iter] <-  res.m1$pvalue
		m2.pvalue[iter] <-  res.m2$pvalue
		m3.pvalue[iter] <-  res.m3$pvalue
		m4.pvalue[iter] <-  res.m4$pvalue
		m5.pvalue[iter] <-  res.m5$pvalue
		parameters=c(ef$beta_gc,ef$avbetaf1,ef$avbetaf2,ef$betaf1f2,ef$avbetaf1f2)
		a=rbind(a,parameters)
		nummut=sum(mut)
		b=c(b,nummut)
		}
		return(list("m1.pvalue" =m1.pvalue, "m2.pvalue" =m2.pvalue,
               "m3.pvalue" =m3.pvalue,"m4.pvalue" =m4.pvalue,"m5.pvalue"=m5.pvalue,"#mut"=b))
		}




power_comparebase96 <- function(binary, Niter, sganno,sgmatrix,Nsample,para,signatures,rho,beta_gc,hot=0,hmm,sc){
  m1.pvalue <- rep(1,Niter)
	a=c()
	b=c()
	d=c()
	for (iter in 1:Niter) {
	  print(iter)
		simdata <- simulate_1funcv96(binary=binary,sganno=sganno,sgmatrix=sgmatrix, Nsample=Nsample, beta_gc=beta_gc, para=para,signatures=signatures, rho,hot=hot,hmm=hmm,sc=sc)
		ssgdata=simdata$annodata
		mut <- do.call(rbind, simdata$mutlist)
		bmrmtx <- do.call(rbind, simdata$bmrmtxlist)
		e <- simdata$pheno
		ef <- simdata$efsize
		fe <- ef$avFe
		mr <- bmrmtx
		if (sum(mut) ==0) {next}
		if (binary == F){
		res.m1 <- ddmodel(mut,e, mr, fe)
		}else{
		res.m1 <- ddmodel_binary(mut,e,mr,fe)
		}
		m1.pvalue[iter] <-  res.m1$pvalue
		#m1.pvalue <-  res.m1$pvalue
		all=res.m1$all
		#parameters=c(ef$beta_gc,ef$avbetaf1,ef$avbetaf2,ef$betaf1f2,ef$avbetaf1f2)
		#a=rbind(a,parameters)
		nummut=sum(mut)
		b=c(b,nummut)
		d=rbind(d,all)
		#input=list("mut"=mut,"e"=e,"bmr"=mr,"fe"=fe)
}
   return(list("m1.pvalue" =m1.pvalue,"estimates"=d,"real"=para,"#mut"=b))
  }




power_comparediff96 <- function(binary, Niter, sganno,sgmatrix,Nsample,para,signatures,rho,beta_gc,hot=0,hmm,sc){
  m1.pvalue <- rep(1,Niter)
	a=c()
	b=c()
	d=c()
	for (iter in 1:Niter) {
	  print(iter)
		simdata <- simulate_1funcv96(binary=binary,sganno=sganno,sgmatrix=sgmatrix, Nsample=Nsample, beta_gc=beta_gc, para=para,signatures=signatures, rho=rho,hot=hot,hmm=hmm,sc=sc)
		ssgdata=simdata$annodata
		mut <- do.call(rbind, simdata$mutlist)
		bmrmtx <- do.call(rbind, simdata$bmrmtxlist)
		e <- simdata$pheno
		ef <- simdata$efsize
		fe <- ef$diffFe
		mr <- bmrmtx
		if (sum(mut) ==0) {next}
		if (binary == F){
		res.m1 <- ddmodel(mut,e, mr, fe)
		}else{
		res.m1 <- ddmodel_binary(mut,e,mr,fe)
		}
		m1.pvalue[iter] <-  res.m1$pvalue
		#m1.pvalue <-  res.m1$pvalue
		#parameters=c(ef$beta_gc,ef$avbetaf1,ef$avbetaf2,ef$betaf1f2,ef$avbetaf1f2)
		#a=rbind(a,parameters)
		nummut=sum(mut)
		b=c(b,nummut)
		all=res.m1$all
		d=rbind(d,all) 
		#input=list("mut"=mut,"e"=e,"bmr"=mr,"fe"=fe)
 
}
   return(list("m1.pvalue" =m1.pvalue,"estimates"=d,"real"=para,"#mut"=b))
  }



power_comparediff96Fix <- function(binary, Niter, sganno,sgmatrix,Nsample,para,signatures,rho,beta_gc,beta_gcFix,hot=0,hmm,sc){
  m1.pvalue <- rep(1,Niter)
	a=c()
	b=c()
	d=c()
	for (iter in 1:Niter) {
	  print(iter)
		simdata <- simulate_1funcv96(binary=binary,sganno=sganno,sgmatrix=sgmatrix, Nsample=Nsample, beta_gc=beta_gc,beta_gcFix=beta_gcFix, para=para,signatures=signatures, rho=rho,hot=hot,hmm=hmm,sc=sc)
		ssgdata=simdata$annodata
		mut <- do.call(rbind, simdata$mutlist)
		bmrmtx <- do.call(rbind, simdata$bmrmtxlist)
		e <- simdata$pheno
		ef <- simdata$efsize
		fe <- ef$diffFeFix
		mr <- bmrmtx
		if (sum(mut) ==0) {next}
		if (binary == F){
		res.m1 <- ddmodel(mut,e, mr, fe)
		}else{
		res.m1 <- ddmodel_binary(mut,e,mr,fe)
		}
		m1.pvalue[iter] <-  res.m1$pvalue
		#m1.pvalue <-  res.m1$pvalue
		#parameters=c(ef$beta_gc,ef$avbetaf1,ef$avbetaf2,ef$betaf1f2,ef$avbetaf1f2)
		#a=rbind(a,parameters)
		nummut=sum(mut)
		b=c(b,nummut)
		all=res.m1$all
		d=rbind(d,all) 
		#input=list("mut"=mut,"e"=e,"bmr"=mr,"fe"=fe)
 
}
   return(list("m1.pvalue" =m1.pvalue,"estimates"=d,"real"=para,"#mut"=b))
  }



