load("data/pppower_bmr_variations.Rd")
SIZES <- bmrvar$sizes
MCOL <- c(Linear="#00AFBB", Fisher="#2b8cbe", Binomial="#aa00ff", LogisticR="#0001f6",
          MannWhitney="#1a9850", Coselens="#ff7f00", DiffDriver="#FF3D2E")
MPCH <- c(Linear=19, Fisher=17, Binomial=15, LogisticR=18, MannWhitney=3, Coselens=8, DiffDriver=17)
MLAB <- c("Linear regression","Fisher","Binomial","Logistic reg.","Mann–Whitney","coselens","DiffDriver")
al_of <- function(d) 1/nrow(d)
fpN   <- function(vlist, m) sapply(as.character(SIZES), function(n){d<-vlist[[n]]; mean(d[[m]]<al_of(d),na.rm=TRUE)})
CA<-0.82; CM<-0.9; CL<-1.2               # cex: axis / title / panel-letter

sub_fpN <- function(vlist, ttl, letter=NULL, xt=c(1,4,7)){
  x<-seq_along(SIZES); al<-al_of(vlist[["1500"]])
  plot(NA,xlim=range(x),ylim=c(0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  axis(1,at=x[xt],labels=SIZES[xt],cex.axis=CA,mgp=c(2,0.45,0),tcl=-0.25)
  axis(2,at=seq(0,1,.25),cex.axis=CA,las=1,mgp=c(2,0.55,0),tcl=-0.25)
  grid(col="grey88"); abline(h=al,lty=3,col="grey50")
  for(m in names(MCOL)){y<-fpN(vlist,m);lines(x,y,col=MCOL[m],lwd=1.3);points(x,y,col=MCOL[m],pch=MPCH[m],cex=0.85)}
  mtext(ttl,side=3,line=0.25,cex=CM)
  if(!is.null(letter)) mtext(letter,side=3,line=0.9,adj=-0.30,font=2,cex=CL)
}
draw <- function(){
  m <- rbind(c(1,2,3,4),      # a: R2 sweep (0.3,0.5,0.78,0.95)
             c(5,6,7,10),     # b: burden (tmb 0,0.5,0.9) + legend cell(10)
             c(8,8,9,9))      # c: continuous (span 2) ; d: k (span 2)
  layout(m, heights=c(1,1,1))
  par(mar=c(2.3,2.5,1.9,0.6), oma=c(1.9,2.4,0.3,0.3), cex=1)   # cex=1: no auto-shrink

  ## (a) confounding strength — the four R^2 panels
  sub_fpN(bmrvar$partA$R2_0.3, "R² = 0.3", "a")
  sub_fpN(bmrvar$partA$R2_0.5, "R² = 0.5")
  sub_fpN(bmrvar$partA$baseline, "R² = 0.78")
  sub_fpN(bmrvar$partA$R2_0.95, "R² = 0.95")
  ## (b) mutation-burden heterogeneity — the three tmb_cv panels (alpha=1/Niter differs)
  sub_fpN(bmrvar$partA$baseline,  "tmb_cv = 0  (α=.01)", "b")
  sub_fpN(bmrvar$partA$tmb_cv0.5, "tmb_cv = 0.5  (α=.0067)")
  sub_fpN(bmrvar$partA$tmb_cv0.9, "tmb_cv = 0.9  (α=.004)")
  ## (c) continuous phenotype
  sub_fpN(bmrvar$partA$continuous, "Continuous phenotype, R² = 0.78  (α=.02)", "c")
  ## (d) number of fitted signatures
  ks<-names(bmrvar$partB); kcols<-c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e")
  al<-al_of(bmrvar$partB[[1]][["1500"]]); x<-seq_along(SIZES)
  plot(NA,xlim=range(x),ylim=c(0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  axis(1,at=x,labels=SIZES,cex.axis=CA,mgp=c(2,0.45,0),tcl=-0.25)
  axis(2,at=seq(0,1,.25),cex.axis=CA,las=1,mgp=c(2,0.55,0),tcl=-0.25)
  grid(col="grey88"); abline(h=al,lty=3,col="grey50")
  for(i in seq_along(ks)){y<-sapply(as.character(SIZES),function(n){d<-bmrvar$partB[[ks[i]]][[n]];mean(d$diff_est<al_of(d),na.rm=TRUE)})
    lines(x,y,col=kcols[i],lwd=1.5);points(x,y,col=kcols[i],pch=19,cex=0.9)}
  mtext("Number of fitted signatures, DiffDriver (truth = 6)",side=3,line=0.25,cex=CM)
  mtext("d",side=3,line=0.9,adj=-0.14,font=2,cex=CL)
  legend("topright",paste("k =",sub("k","",ks)),col=kcols,pch=19,lwd=1.5,bty="n",cex=0.78,ncol=2)
  ## shared 7-method legend (cell 10)
  par(mar=c(0,0,0,0)); plot.new()
  legend("center", MLAB, col=MCOL, pch=MPCH, lwd=1.4, lty=1, bty="n", cex=0.8, title="Method", seg.len=1.2)
  ## shared axis titles
  mtext("Number of samples", side=1, outer=TRUE, line=0.6, cex=CM)
  mtext("False positive rate (p < α = 1/Niter)", side=2, outer=TRUE, line=0.9, cex=CM)
}
OUT <- getwd()   
pdf(file.path(OUT,"SuppFig3_confounding_robustness.pdf"), width=7.2, height=7.4, pointsize=10); draw(); dev.off()
png(file.path(OUT,"SuppFig3_confounding_robustness.png"), width=7.2, height=7.4, units="in", res=200, pointsize=10); draw(); dev.off()
cat("done\n")
