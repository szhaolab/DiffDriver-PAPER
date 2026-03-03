# TODO: Generate the hotspot sequence
#' @title Generate the hotspot sequence
#' @param hmm transition probability vector
#' @description Generate hotspot sequence
#' @param sgdata annotation list
#' @return A n by 2  data frame. The first column start provide the position information, the second column seqt indicates whether this position being a hotspot.
#' @noRd
hotspotseq=function(hmm,sgdata){
  pos=sort(unique(do.call(rbind,sgdata)$start)) # distract the positions
  seqt=sample(x=c(0,1),size=1,prob =c(hmm[1],hmm[2])) # generate the inital state
  for (i in 2:length(pos)) { # generate the following states
    a=ifelse(seqt[i-1]==0,sample(x=c(0,1),size=1,prob =c(hmm[3],hmm[4])),
             sample(x=c(0,1),size=1,prob =c(hmm[5],hmm[6])))
    seqt=c(seqt,a)
  }
  hotseq=cbind(start=pos,seqt=seqt)
  return(hotseq)
}
