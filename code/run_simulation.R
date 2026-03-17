#' Run a single simulation iteration: simulate mutations, then test
run_one_iteration <- function(binary, sganno, sgmatrix, bmrpars, betaf0, Nsample, beta_gc, para, rho, tau, hot, hmm) {
  simdata <- simulate_selection(
    binary = binary, sganno = sganno, sgmatrix = sgmatrix,
    bmrpars = bmrpars, betaf0 = betaf0, Nsample = Nsample,
    beta_gc = beta_gc, para = para, rho = rho, tau = tau,
    hot = hot, hmm = hmm
  )
  mut <- do.call(rbind, simdata$mutlist)
  mr  <- do.call(rbind, simdata$bmrmtxlist)
  e   <- simdata$pheno
  fe  <- simdata$efsize$diffFe
  nummut <- sum(mut)

  if (nummut == 0) {
    return(list(pvalue = NA, nummut = 0))
  }

  if (binary == FALSE) {
    res <- ddmodel(mut, e, mr, fe)
  } else {
    res <- ddmodel_binary(mut, e, mr, fe)
  }

  list(pvalue = res$pvalue, nummut = nummut)
}

#' Run simulation over multiple iterations
run_simulation <- function(binary, Niter, sganno, sgmatrix, Nsample, para, rho, tau = 1, bmrpars, betaf0, beta_gc, hot = 0, hmm) {
  pvalues <- rep(NA, Niter)
  nummuts <- rep(0, Niter)

  for (iter in 1:Niter) {
    print(iter)
    res <- run_one_iteration(
      binary = binary, sganno = sganno, sgmatrix = sgmatrix,
      bmrpars = bmrpars, betaf0 = betaf0, Nsample = Nsample,
      beta_gc = beta_gc, para = para, rho = rho, tau = tau,
      hot = hot, hmm = hmm
    )
    pvalues[iter] <- res$pvalue
    nummuts[iter] <- res$nummut
  }

  list(pvalues = pvalues, nummuts = nummuts)
}
