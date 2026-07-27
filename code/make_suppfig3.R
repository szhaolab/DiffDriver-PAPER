load("data/pppower_bmr_variations.Rd")
SIZES <- bmrvar$sizes
MCOL <- c(Linear="#00AFBB", Fisher="#2b8cbe", Binomial="#aa00ff", LogisticR="#0001f6",
          MannWhitney="#1a9850", Coselens="#ff7f00", DiffDriver="#FF3D2E")
MPCH <- c(Linear=19, Fisher=17, Binomial=15, LogisticR=18, MannWhitney=3, Coselens=8, DiffDriver=17)
MLAB <- c("Linear regression","Fisher","Binomial","Logistic reg.","Mann–Whitney","coselens","DiffDriver")
al_of  <- function(d) 1/nrow(d)
fpN    <- function(vlist, m) sapply(as.character(SIZES), function(n){d<-vlist[[n]]; mean(d[[m]]<al_of(d),na.rm=TRUE)})
panel_letter <- function(L) mtext(L, side=3, line=0.6, adj=-0.18, font=2, cex=1.15)

draw <- function(){
  par(mfrow=c(2,2), mar=c(3.5,3.7,2.1,0.6), mgp=c(2.15,0.65,0), tcl=-0.3)

  ## (a) confounding strength: FPR vs R^2 at N=1500 (100 iters -> alpha=0.01)
  r2v <- c(0.3,0.5,0.78,0.95)
  r2d <- list(bmrvar$partA$R2_0.3, bmrvar$partA$R2_0.5, bmrvar$partA$baseline, bmrvar$partA$R2_0.95)
  al  <- al_of(r2d[[1]][["1500"]])
  plot(NA, xlim=range(r2v), ylim=c(0,1), xlab="Phenotype–signature R²",
       ylab=sprintf("False positive rate (p < %.4g)",al), las=1, cex.axis=0.9, cex.lab=1.0)
  grid(); abline(h=al, lty=3, col="grey55")
  for (m in names(MCOL)) { y<-sapply(r2d, function(v) mean(v[["1500"]][[m]]<al_of(v[["1500"]]),na.rm=TRUE))
    lines(r2v,y,col=MCOL[m],lwd=1.5); points(r2v,y,col=MCOL[m],pch=MPCH[m],cex=1.1) }
  title("Confounding strength (N = 1,500)", cex.main=1.0, font.main=1)
  legend("topleft", MLAB, col=MCOL, pch=MPCH, lwd=1.5, lty=1, bty="n", cex=0.8, ncol=1, seg.len=1.3)
  panel_letter("a")

  ## (b) mutation-burden heterogeneity: FPR vs N at tmb_cv=0.9 (250 iters -> 0.004)
  v<-bmrvar$partA$tmb_cv0.9; al<-al_of(v[["1500"]]); x<-seq_along(SIZES)
  plot(NA, xlim=range(x), ylim=c(0,1), xaxt="n", xlab="Number of samples",
       ylab=sprintf("False positive rate (p < %.4g)",al), las=1, cex.axis=0.9, cex.lab=1.0)
  axis(1,x,SIZES,cex.axis=0.9); grid(); abline(h=al,lty=3,col="grey55")
  for (m in names(MCOL)){ y<-fpN(v,m); lines(x,y,col=MCOL[m],lwd=1.5); points(x,y,col=MCOL[m],pch=MPCH[m],cex=1.1) }
  title(expression("Mutation-burden heterogeneity (tmb"[cv]*" = 0.9)"), cex.main=1.0, font.main=1)
  panel_letter("b")

  ## (c) context type: continuous phenotype, FPR vs N (50 iters -> 0.02)
  v<-bmrvar$partA$continuous; al<-al_of(v[["1500"]]); x<-seq_along(SIZES)
  plot(NA, xlim=range(x), ylim=c(0,1), xaxt="n", xlab="Number of samples",
       ylab=sprintf("False positive rate (p < %.4g)",al), las=1, cex.axis=0.9, cex.lab=1.0)
  axis(1,x,SIZES,cex.axis=0.9); grid(); abline(h=al,lty=3,col="grey55")
  for (m in names(MCOL)){ y<-fpN(v,m); lines(x,y,col=MCOL[m],lwd=1.5); points(x,y,col=MCOL[m],pch=MPCH[m],cex=1.1) }
  title("Continuous phenotype (R² = 0.78)", cex.main=1.0, font.main=1)
  panel_letter("c")

  ## (d) number of fitted signatures: DiffDriver FPR vs N, k in {4,5,6,8,10}
  ks<-names(bmrvar$partB); kcols<-c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e")
  al<-al_of(bmrvar$partB[[1]][["1500"]]); x<-seq_along(SIZES)
  plot(NA, xlim=range(x), ylim=c(0,1), xaxt="n", xlab="Number of samples",
       ylab=sprintf("DiffDriver FPR (p < %.4g)",al), las=1, cex.axis=0.9, cex.lab=1.0)
  axis(1,x,SIZES,cex.axis=0.9); grid(); abline(h=al,lty=3,col="grey55")
  for (i in seq_along(ks)){ y<-sapply(as.character(SIZES),function(n){d<-bmrvar$partB[[ks[i]]][[n]]; mean(d$diff_est<al_of(d),na.rm=TRUE)})
    lines(x,y,col=kcols[i],lwd=1.6); points(x,y,col=kcols[i],pch=19,cex=1.1) }
  title("Number of fitted signatures (truth = 6)", cex.main=1.0, font.main=1)
  legend("topleft", paste("k =",sub("k","",ks)), col=kcols, pch=19, lwd=1.6, bty="n", cex=0.8)
  panel_letter("d")
}

OUT <- getwd()   
pdf(file.path(OUT,"SuppFig3_confounding_robustness.pdf"), width=7.2, height=6.2, pointsize=12); draw(); dev.off()
png(file.path(OUT,"SuppFig3_confounding_robustness.png"), width=7.2, height=6.2, units="in", res=200, pointsize=12); draw(); dev.off()
cat("wrote SuppFig3 pdf+png\n")
