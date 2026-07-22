# =============================================================================
# Estimated-background helpers for the signature-confounding simulation.
#
# The confounding simulation (code/simulate_functions.R) can hand DiffDriver and
# coselens the TRUE background rate. These helpers instead let each method
# ESTIMATE the background from a simulated genome-wide SILENT (synonymous)
# mutation count matrix, calibrated to the real LUAD synonymous burden.
#
#   build_silent_matrix()     simulate the 96 x Nsample silent count matrix Y
#   coselens_estimated_rate() per-group 192-context rate from Y (for dndscv_truerate)
#   diffdriver_estimated_mr() per-site log-background from fastTopics(Y)
#
# All three take a single confounding-sim cohort `sim` (from simulate_confounding)
# and the site inventory `Lsyn96` (data/Lsyn.rds; synonymous site counts per
# 96 trinucleotide context, from dndscv's RefCDS). The 96->192 mapping used by
# coselens is data/context96_map.rds.
#
# Because the true bmr already yields ~63 synonymous mutations per sample and
# LUAD's observed burden is ~64, the calibration factor g is ~1.02 -- the silent
# count matrix and the ERBB3 mutations are on the same realistic scale, so no
# rescaling of the ERBB3 sim is needed.
# =============================================================================

#' Simulate the genome-wide silent (synonymous) mutation count matrix
#'
#' Y[c,s] ~ Poisson(Lsyn96[c] * bmr[c,s] * g), where g calibrates the mean
#' per-sample synonymous total to `target_silent` (the LUAD observed burden).
#'
#' @param bmr 96 x Nsample true per-site background probability (sim$bmrfold$bmr)
#' @param Lsyn96 Length-96 synonymous site counts per context (data/Lsyn.rds$Lsyn96)
#' @param target_silent Target mean synonymous mutations per sample (LUAD ~64.2)
#' @return 96 x Nsample integer matrix of simulated silent counts
build_silent_matrix <- function(bmr, Lsyn96, target_silent = 64.2) {
  N <- ncol(bmr)
  g <- target_silent / mean(colSums(Lsyn96 * bmr))
  lambda <- (Lsyn96 %o% rep(1, N)) * bmr * g
  matrix(rpois(length(lambda), lambda), nrow = 96)
}

#' coselens ESTIMATED per-group substitution rate from the silent count matrix
#'
#' The per-context rate is counts/sites (exactly what dndscv's Poisson
#' substitution model fits): rate96[c] = sum over the group's samples of Y[c,]
#' divided by the synonymous site count Lsyn96[c]. This is the true rate up to
#' the (near-1) calibration constant plus Poisson noise. The 96-vector is mapped
#' to the 192 dndscv contexts (each 96-context maps to two) via context96_map,
#' ready to pass to run_coselens(..., true_rate = list(mutrates1 = , mutrates2 = )).
#'
#' @param Y 96 x Nsample silent count matrix (build_silent_matrix)
#' @param group_idx Integer indices of the samples in this group
#' @param Lsyn96 Length-96 synonymous site counts per context
#' @param context96_map data.frame with column `number` (96-context index per
#'   192 dndscv context, in substmodel order); data/context96_map.rds
#' @return Length-192 per-site rate vector in dndscv substmodel order
coselens_estimated_rate <- function(Y, group_idx, Lsyn96, context96_map) {
  rate96 <- rowSums(Y[, group_idx, drop = FALSE]) / Lsyn96
  rate96[context96_map$number]
}

#' DiffDriver ESTIMATED per-site log-background from the silent count matrix
#'
#' Runs the diffdriver signature BMR (fastTopics on the Nsample x 96 count
#' matrix) exactly as diffdriver::matrixlistToBMR does, reconstructing any
#' all-zero contexts (fastTopics drops them) and low-count samples, to get the
#' per-sample x 96 spectrum. The spectrum is scale-matched to the true bmr
#' magnitude (a common constant; absorbed by the model baseline and cancelling
#' in the case/control test) and expanded to a per-site log-rate list matching
#' `bmrmtxlist` (one matrix per nucleotide-type segment).
#'
#' @param Y 96 x Nsample silent count matrix
#' @param bmrmtxlist The sim's per-segment true log-background list (sim$bmrmtxlist);
#'   used only for its per-segment dimensions
#' @param bmr 96 x Nsample true per-site rate (sim$bmrfold$bmr); for scale-matching
#' @param k Number of signatures for fastTopics (default 6)
#' @return rbind-able list of per-segment log-background matrices, or NULL on failure
diffdriver_estimated_mr <- function(Y, bmrmtxlist, bmr, k = 6) {
  N <- ncol(Y)
  ymat <- t(Y); colnames(ymat) <- 1:96; rownames(ymat) <- paste0("S", 1:N)
  keep <- rowSums(ymat) >= 5
  ft <- tryCatch(fastTopics::fit_topic_model(ymat[keep, ], k = k), error = function(x) NULL)
  if (is.null(ft)) return(NULL)
  # reconstruct dropped contexts (as matrixlistToBMR) and dropped samples
  ff <- matrix(min(ft$F) / 2, 96, k); ff[as.numeric(rownames(ft$F)), ] <- ft$F
  ll <- matrix(rep(colMeans(ft$L), N), N, byrow = TRUE)
  rownames(ll) <- rownames(ymat); ll[rownames(ft$L), ] <- ft$L
  sig_s <- ll %*% t(ff)                       # Nsample x 96 spectrum
  sig_s <- sig_s * (mean(bmr) / mean(sig_s))  # scale-match to true bmr magnitude
  lapply(seq_len(96), function(tt) {
    nr <- nrow(bmrmtxlist[[tt]])
    if (is.null(nr) || nr == 0) return(bmrmtxlist[[tt]])
    log(matrix(1, nr, 1) %x% matrix(sig_s[, tt], nrow = 1))
  })
}
