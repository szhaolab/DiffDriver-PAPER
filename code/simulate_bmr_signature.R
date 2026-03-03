#' @param e phenotype vector
#' @param sigatures A 96 times m matrix M. Each column of it is a signaturee.
#' @param rho Correlation between e and the loadings of first signature.
#' @param s Scale parameter
#' @keywords internal
#' @return matrix bmr
#' @noRd
bmrSignature=function(e,signatures,rho,sc=1,adjustment=F){
  if (adjustment==T) {
    w=(N+1)/sum((N+1))
  }else{
    w=rep(1,96)
  } # weight based on the number of silent mutations
  nn=length(e) ## sample size
  if (nrow(signatures)!=96) {stop("Signature length should be 96!")}
  m=ncol(signatures) ## number of signatures
  cc= matrix(runif(m*nn),nrow = m, ncol = nn) ## loading matrix
  # for (i in 1:nrow(cc)) {
  #   cc[i,]=cc[i,]/mean(cc[i,])
  # }
  complement <- function(y, rho, x) {
    if (missing(x)) x <-runif(length(y), min=0, max =1) # Optional: supply a default if `x` is not given
    y.perp <- residuals(lm(x ~ y))
    rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
  }
  b=complement(e,rho)
  b.new=b-min(b)+0.1
  cc[1,]=b.new/mean(b.new)
  yy=as.matrix(signatures)%*%as.matrix(cc)
  bmr=diag(1/w)%*%yy/sc
  return(list("bmr"=bmr,"loadings"=cc,"signatures"=signatures))
}
