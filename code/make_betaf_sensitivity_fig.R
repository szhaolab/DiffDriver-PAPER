SP <- "/private/tmp/claude-501/-Users-siming-Dartmouth-College-Dropbox-Siming-Zhao-Project-diffDriver-DiffDriver-PAPER/46eb8cab-f9ce-420f-8359-609ff7701606/scratchpad"
SP <- "."   # dir with betaf_sensitivity_metrics.rds; writes betaf_sensitivity.{pdf,png}
SIZES <- c(200,400,600,800,1000,1200,1500)
CAVG<-"#FF3D2E"; CTRU<-"#1f78b4"           # averaged (red) vs ground-truth (blue)
titles <- c(TSG_binary_0.8="Tumor-suppressor genes (TSG)", OG_binary_0.8="Oncogenes (OG)")
draw <- function(){
  par(mfrow=c(2,2), mar=c(3.6,4.0,2.2,0.8), mgp=c(2.3,0.7,0), cex=1)
  for (metric in c("power","fdr")){
    for (sn in c("TSG_binary_0.8","OG_binary_0.8")){
      d <- out[out$setting==sn,]; x <- seq_along(SIZES)
      ya <- d[[paste0(metric,"_avg")]]; yt <- d[[paste0(metric,"_truth")]]
      ymax <- if (metric=="power") max(0.35, yt)*1.05 else 1
      plot(NA,xlim=range(x),ylim=c(0,ymax),xaxt="n",las=1,
           xlab="Number of samples",
           ylab=if(metric=="power")"Power (FDR = 0.1)" else "Observed FDR",
           main=titles[[sn]], cex.main=1.0, cex.lab=0.95, cex.axis=0.9)
      axis(1,x,SIZES,cex.axis=0.7); grid(col="grey90")
      if (metric=="fdr") abline(h=0.1,lty=3,col="grey45")
      lines(x,ya,col=CAVG,lwd=1.8); points(x,ya,col=CAVG,pch=19,cex=1.2)
      lines(x,yt,col=CTRU,lwd=1.8,lty=2); points(x,yt,col=CTRU,pch=17,cex=1.2)
      if (metric=="power" && sn=="TSG_binary_0.8")
        legend("topleft", c(expression(paste("Averaged ", beta^f, " (imported, DiffDriver)")),
                            expression(paste("Ground-truth ", beta^f, " (context-specific)"))),
               col=c(CAVG,CTRU), pch=c(19,17), lty=c(1,2), lwd=1.8, bty="n", cex=0.78)
    }
  }
}
pdf(file.path(SP,"betaf_sensitivity.pdf"), width=7.2, height=6.4, pointsize=11); draw(); dev.off()
png(file.path(SP,"betaf_sensitivity.png"), width=7.2, height=6.4, units="in", res=200, pointsize=11); draw(); dev.off()
cat("done\n")
