#'
#' @param Argument binary indicate whether the phenetype is binary or not. 
#' @param List sgdata provide the annotation data 
#' @param Vector bmrpars provides nine parameters for each of the ny type. 
#' @param Scalar betaf0 is the gene effect. 
#' @param Scalar Nsample specifies the sample size
#' @param Two-compoent vector beta_gc specify the effect of functional type 7 and 8. 
#' @param The vector para are the parameters for generating under-selection genes from phenotype.   
#' @param The scalar hot indicate whether there are hotspot effect, in which hot=0 stands for no hotspots while hot=1 hotspots.  
#' @param The vector hmm provides probabilities for generating the hotspot sequence. 
#'
#' @return A list including the mutation data, phenotype daga, annotation data and true parameters.  
#' @export
#'
#' @examples
simulate_1funcv <- function(binary=F,sgdata, bmrpars, betaf0=2, Nsample, beta_gc, para,hot=0, hmm){
	if (binary==T){ # generate binary phenotype 
		Nsamplec <- round(Nsample/2) # number of samples with phenotype E=1 (the rest will be 0)
		Nsamplen <- Nsample-Nsamplec
		phenotype <- c(rep(1,Nsamplec),rep(0,Nsamplen))
		ss=ifelse(phenotype==1,sample(c(0,1),size=Nsamplec,replace=T,prob = c(1-para[1],para[1])),sample(c(0,1),size=Nsamplen,replace=T,prob = c(1-para[2],para[2]))) # generate positive samples
		}else{ # generate conitnuous phenotype
		phenotype=rnorm(Nsample,mean = para[1],sd=para[2])
		pp=exp(para[3]+para[4]*phenotype)/(1+exp(para[3]+para[4]*phenotype))
		ss=ifelse(runif(Nsample)-pp<0,1,0) # generate positive samples
		}
		index=which(ss==1)
		phenotype=c(phenotype[index],phenotype[-index])# the positive samples are placed in the front.
		Nsample.ps=sum(ss) # number of positive samples
		Nsample.neu <- Nsample - Nsample.ps # number of neutral samples 
		hotseq=hotspotseq(hmm,sgdata) # The first column'start' is the positions, and the second column 'seqt',is the hotspot indicator.  
		if (hot==0){
			hotseq[,2]=0
			}
	mutlist <- list() # a list of nine mutation matrices
	countlist <- list()
	annodata <- list() # a list of nine annotation data frames
	bmrmtxlist <- list() # a list of nine background mutation matrices
	hotsize <- list() # a list of nine hotspot sequences. 
	for (t in 1:length(sgdata)) {
		ssgdata=merge(sgdata[[t]],hotseq,by="start")
		tnpos1 <- dim(ssgdata[functypecode==7])[1]
		tnpos2 <- dim(ssgdata[functypecode==8])[1]
		k1=ssgdata[functypecode==7 & seqt==1,.N]
		k3=ssgdata[functypecode==8 & seqt==1,.N]
		k2=tnpos1-k1 
		k4=tnpos2-k3
		pp2=exp(bmrpars[t])*exp(betaf0)*exp(beta_gc[1]) # mutation rate for functypecode=7
		pp3=exp(bmrpars[t])*exp(betaf0)*exp(beta_gc[1] + beta_gc[2])  # mutation rate for functypecode=8
		pp1=exp(bmrpars[t])*exp(betaf0) # background mutation rate
		hotpp= min(pp1*exp(hmm[9]),1) # mutation rate for hotspot positions
		annodata[[t]] <- rbind(ssgdata[functypecode==7], ssgdata[functypecode==8]) # order the annodata based on functional type
		if (k1>0){ # generate mutation for functional type=7 & positive samples
			mut11= rsparsematrix(k1,Nsample.ps,nnz = rbinom(1, Nsample.ps * k1, hotpp),rand.x=NULL)
			mut12= rsparsematrix(k2,Nsample.ps,nnz = rbinom(1, Nsample.ps * k2, pp2),rand.x=NULL)
			mut1=rbind(mut11,mut12)
			}else{
				mut1= rsparsematrix(k2,Nsample.ps,nnz = rbinom(1, Nsample.ps * k2, pp2),rand.x=NULL)
			}
		if (k3>0){# generate mutation for functional type=8 & positive samples
			mut21= rsparsematrix(k3,Nsample.ps,nnz = rbinom(1, Nsample.ps * k3, hotpp),rand.x=NULL)
			mut22= rsparsematrix(k4,Nsample.ps,nnz = rbinom(1, Nsample.ps * k4, pp3),rand.x=NULL)
			mut2=rbind(mut21,mut22)
			}else{
				mut2=  rsparsematrix(k4,Nsample.ps,nnz = rbinom(1, Nsample.ps * k4, pp3),rand.x=NULL)
				}
		mut.ps=rbind(mut1,mut2) # mutation matrix for positive samples 
		mut.neu= rsparsematrix(tnpos1+tnpos2,Nsample.neu,nnz = rbinom(1, Nsample.neu * (tnpos1+tnpos2), pp1),rand.x=NULL) # mutation for neutral samples
		mutlist[[t]] <- cbind(mut.ps,mut.neu) # mutation matrix for all samples
		countlist[[t]] <- c(tnpos1, tnpos2,sum(mut.ps),sum(mut.neu))
		bmrmtxlist[[t]] <- matrix(bmrpars[t], ncol = ncol(mutlist[[t]]), nrow = nrow(mutlist[[t]])) # background mutation matrix for nytype=t 
		hotsize[[t]]=hmm[9]*c(rep(1,k1),rep(0,k2),rep(1,k3),rep(0,k4)) # hot spot sequence for nytype=t
}

# The forllowings are the ture parameters
	avbetaf1 <- log(exp(beta_gc[1]) * Nsample.ps/Nsample + Nsample.neu/Nsample)
	avbetaf2 <- log(exp(beta_gc[1] + beta_gc[2]) * Nsample.ps/Nsample + Nsample.neu/Nsample)
	pos1pos2ratio <- colSums(do.call(rbind, countlist))[1]/colSums(do.call(rbind, countlist))[2]
	avbetaf1f2 <- log((pos1pos2ratio*exp(avbetaf1) + exp(avbetaf2))/(pos1pos2ratio+1))
	betaf1f2 <- log((pos1pos2ratio*exp(beta_gc[1]) + exp(beta_gc[2]))/(pos1pos2ratio+1))
	simdata <- list("mutlist"= mutlist, "hotsize"=hotsize, "pheno" = phenotype, "annodata" = annodata, "bmrpars" = bmrpars, "bmrmtxlist" = bmrmtxlist, "para"=para, "efsize" = list( "betaf0" = betaf0,  "beta_gc" = beta_gc, "avbetaf1" = avbetaf1, "avbetaf2" = avbetaf2, "avbetaf1f2" = avbetaf1f2,"betaf1f2"=betaf1f2),"nsample"=c(Nsample.ps,Nsample.neu))
	return(simdata)
}


