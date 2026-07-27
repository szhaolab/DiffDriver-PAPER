load("data/pppower_bmr_variations.Rd")
SIZES <- bmrvar$sizes
MCOL <- c(Linear="#00AFBB", Fisher="#aa00ff", Binomial="#0001f6", LogisticR="#fdae61",
          MannWhitney="#a65628", Coselens="#1a9850", DiffDriver="#FF3D2E")
MPCH <- c(Linear=17, Fisher=15, Binomial=18, LogisticR=8, MannWhitney=3, Coselens=4, DiffDriver=19)
MLAB <- c("Linear regression","Fisher","Binomial","Logistic reg.","Mann–Whitney","Coselens","DiffDriver")
al_of <- function(d) 1/nrow(d)
fpN   <- function(vlist, m) sapply(as.character(SIZES), function(n){d<-vlist[[n]]; mean(d[[m]]<al_of(d),na.rm=TRUE)})
CA<-0.8; CM<-0.88; CLAB<-0.8; CL<-1.15    # cex: axis ticks / title / axis-label / panel-letter

sub_fpN <- function(vlist, ttl, letter=NULL, xt=c(1,4,7)){
  x<-seq_along(SIZES); al<-al_of(vlist[["1500"]])
  plot(NA,xlim=range(x),ylim=c(0,1),xaxt="n",yaxt="n",
       xlab="Number of samples", ylab="FPR", cex.lab=CLAB)
  axis(1,at=x[xt],labels=SIZES[xt],cex.axis=CA,mgp=c(2,0.45,0),tcl=-0.25)
  axis(2,at=seq(0,1,.25),cex.axis=CA,las=1,mgp=c(2,0.55,0),tcl=-0.25)
  grid(col="grey88"); abline(h=al,lty=3,col="grey50")
  for(m in names(MCOL)){y<-fpN(vlist,m);lines(x,y,col=MCOL[m],lwd=1.3);points(x,y,col=MCOL[m],pch=MPCH[m],cex=0.85)}
  mtext(ttl,side=3,line=0.25,cex=CM)
  if(!is.null(letter)) mtext(letter,side=3,line=0.9,adj=-0.42,font=2,cex=CL,xpd=NA)
}
draw <- function(){
  m <- rbind(c(1,2,3,4), c(5,6,7,10), c(8,8,9,9))
  layout(m, heights=c(1,1,1))
  par(mar=c(3.1,3.5,1.9,0.7), oma=c(0.3,0.3,0.3,0.3), cex=1, mgp=c(1.9,0.5,0))
  ## (a) confounding strength — four R^2 panels
  sub_fpN(bmrvar$partA$R2_0.3, "R² = 0.3", "a")
  sub_fpN(bmrvar$partA$R2_0.5, "R² = 0.5")
  sub_fpN(bmrvar$partA$baseline, "R² = 0.78")
  sub_fpN(bmrvar$partA$R2_0.95, "R² = 0.95")
  ## (b) burden — three tmb_cv panels
  sub_fpN(bmrvar$partA$baseline,  "tmb_cv = 0", "b")
  sub_fpN(bmrvar$partA$tmb_cv0.5, "tmb_cv = 0.5")
  sub_fpN(bmrvar$partA$tmb_cv0.9, "tmb_cv = 0.9")
  ## (c) continuous
  sub_fpN(bmrvar$partA$continuous, "Continuous context, R² = 0.78", "c")
  ## (d) number of fitted signatures
  ks<-names(bmrvar$partB); kcols<-c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e")
  al<-al_of(bmrvar$partB[[1]][["1500"]]); x<-seq_along(SIZES)
  plot(NA,xlim=range(x),ylim=c(0,1),xaxt="n",yaxt="n",
       xlab="Number of samples", ylab="DiffDriver FPR", cex.lab=CLAB)
  axis(1,at=x,labels=SIZES,cex.axis=CA,mgp=c(2,0.45,0),tcl=-0.25)
  axis(2,at=seq(0,1,.25),cex.axis=CA,las=1,mgp=c(2,0.55,0),tcl=-0.25)
  grid(col="grey88"); abline(h=al,lty=3,col="grey50")
  for(i in seq_along(ks)){y<-sapply(as.character(SIZES),function(n){d<-bmrvar$partB[[ks[i]]][[n]];mean(d$diff_est<al_of(d),na.rm=TRUE)})
    lines(x,y,col=kcols[i],lwd=1.5);points(x,y,col=kcols[i],pch=19,cex=0.9)}
  mtext("Number of fitted signatures (truth = 6)",side=3,line=0.25,cex=CM)
  mtext("d",side=3,line=0.9,adj=-0.16,font=2,cex=CL,xpd=NA)
  legend("topright",paste("k =",sub("k","",ks)),col=kcols,pch=19,lwd=1.5,bty="n",cex=0.78,ncol=2)
  ## shared 7-method legend (cell 10)
  par(mar=c(0,0,0,0)); plot.new()
  legend("center", MLAB, col=MCOL, pch=MPCH, lwd=1.4, lty=1, bty="n", cex=0.82, title="Method", seg.len=1.3)
}
OUT <- "/private/tmp/claude-501/-Users-siming-Dartmouth-College-Dropbox-Siming-Zhao-Project-diffDriver-DiffDriver-PAPER/46eb8cab-f9ce-420f-8359-609ff7701606/scratchpad"
pdf(file.path(OUT,"SuppFig4_confounding_robustness.pdf"), width=7.5, height=7.6, pointsize=10); draw(); dev.off()
png(file.path(OUT,"SuppFig4_confounding_robustness.png"), width=7.5, height=7.6, units="in", res=200, pointsize=10); draw(); dev.off()
cat("done\n")
