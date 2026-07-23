#!/usr/bin/env Rscript
# =============================================================================
# simulate_all_genes.R
#
# Consolidates the thousands of per-gene simulation scripts (one .R per gene per
# setting, e.g. phynotypepppower_betaf0=0_sample1200_40_0.9.R) into a single
# driver that loops over ALL genes for a given setting. For each gene it runs the
# same call each per-gene script did:
#
#   set.seed(10)
#   power_compare<method>(binary=, Niter=, sganno=ri200[[gene]],
#       sgmatrix=fanno200[[gene]], bmrpars=log(BMR), Nsample=, betaf0=,
#       beta_gc=, [beta_gcFix=,] para=, hot=, hmm=hmm)
#
# Instead of writing one .Rd per gene, ALL results are collected into a single
# list `simures_all` and saved to ONE file (OUT_FILE, default
# pppower_all_<SETTING>_<METHOD>.Rd). Each element is named by the stem the old
# per-gene files used (e.g. "pppower_betaf0=0_sample1200_40"), so results still
# map 1:1 to those file names -- load(one_file) then index simures_all[[stem]].
#
# A "setting" = one of the result folders. Pick it with SETTING + METHOD below;
# sample sizes and betaf0 are looped within. The power_compare* drivers and the
# simulation function (simulate_selection) come from this repo's
# code/simulate_functions.R; the diffdriver package supplies the underlying
# model/test functions (ddmodel, mlr, ...).
#
# Usage:
#   Rscript simulate_all_genes.R                      # uses the config block below
#   Rscript simulate_all_genes.R TSG_binary_0.9 other # override SETTING METHOD
# =============================================================================

## ------------------------- PATHS (edit to your env) -------------------------
DIFFDRIVER_DIR <- "~/SimingLab/jiezhou/diffdriver"            # diffdriver package only (load_all / ddmodel / mlr ...)
DATA_DIR       <- "data"                                      # ri200.Rd, fanno200.Rd, beta_tsg.Rd, beta_og.Rd, BMR.rda, hmm.rda, parmASHmean.rda
OUT_DIR        <- "."                                         # directory for the single output file
OUT_FILE       <- NULL   # single .Rd to save all results to; NULL = OUT_DIR/pppower_all_<SETTING>_<METHOD>.Rd
USE_LOAD_ALL   <- TRUE   # TRUE: devtools::load_all(DIFFDRIVER_DIR); FALSE: library(diffdriver)
CODE_DIR       <- NULL   # folder holding simulate_functions.R (this repo's code/); NULL = same dir as this script

## ------------------------- SETTING / METHOD --------------------------------
# Presets (one per result folder). Override via command-line args if given.
SETTING <- "TSG_binary_0.8"   # TSG_binary_0.8 | TSG_binary_0.9 | TSG_conti |
                              # OG_binary_0.8  | OG_binary_0.9  | OG_conti
METHOD  <- "diffFix"          # "diffFix" (run_selection_diffdriver) or "other" (run_selection_benchmarks)

# Which sample sizes / betaf0 values to run (defaults match the paper's grids;
# NULL -> filled in from the setting's phenotype below).
SAMPLE_SIZES  <- NULL
BETAF0_VALUES <- NULL
N_GENES       <- 200
SIGNAL_GENES  <- 50           # genes 1..SIGNAL_GENES carry differential selection

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) SETTING <- args[1]
if (length(args) >= 2) METHOD  <- args[2]

## ------------------------- setting definition ------------------------------
# geneset drives the beta file, the parmASHmean index used for beta_gcFix, and
# the hotspot flag; phenotype drives Nite and the para vector; para_var is the
# case-selection prob for signal genes (binary only); fname() encodes the
# folder's output-name convention.
presets <- list(
  TSG_binary_0.8 = list(geneset="TSG", phenotype="binary", para_var=0.8,
    fname=function(b,N,g,pv) sprintf("pppower_betaf0=%s_sample%s_%s.Rd", b, N, g)),
  TSG_binary_0.9 = list(geneset="TSG", phenotype="binary", para_var=0.9,
    fname=function(b,N,g,pv) sprintf("pppower_betaf0=%s_sample%s_%s.Rd", b, N, g)),
  TSG_conti      = list(geneset="TSG", phenotype="conti", para_var=NA, conti_slope=2,
    fname=function(b,N,g,pv) sprintf("pppower_betaf0=%s_sample%s_%s.Rd", b, N, g)),
  # As TSG_conti but with the signal-gene selection slope set to 1 (para = c(0,2,0,1)).
  TSG_conti_slope1 = list(geneset="TSG", phenotype="conti", para_var=NA, conti_slope=1,
    fname=function(b,N,g,pv) sprintf("pppower_betaf0=%s_sample%s_%s.Rd", b, N, g)),
  OG_binary_0.8  = list(geneset="OG",  phenotype="binary", para_var=0.8,
    fname=function(b,N,g,pv) sprintf("pppower_betaf0=%s_sample%s_%s_%s.Rd", b, N, pv, g)),
  OG_binary_0.9  = list(geneset="OG",  phenotype="binary", para_var=0.9,
    fname=function(b,N,g,pv) sprintf("pppower_betaf0=%s_sample%s_%s_%s.Rd", b, N, pv, g)),
  OG_conti       = list(geneset="OG",  phenotype="conti", para_var=NA, conti_slope=2,
    fname=function(b,N,g,pv) sprintf("pppower_betaf0=%s_sample%s_%s.Rd", b, N, g))
)
if (!SETTING %in% names(presets)) stop("Unknown SETTING: ", SETTING,
  " (choose one of: ", paste(names(presets), collapse=", "), ")")
S <- presets[[SETTING]]

geneset   <- S$geneset
phenotype <- S$phenotype
para_var  <- S$para_var
ash_idx   <- if (geneset == "TSG") 1L else 2L          # parmASHmean index
hot       <- if (geneset == "TSG") 0L else 1L          # OG runs used hotspots on
Nite      <- if (phenotype == "binary") 100L else 50L
if (is.null(SAMPLE_SIZES)) {
  SAMPLE_SIZES <- if (phenotype == "binary") c(200,400,600,800,1000,1200,1500) else c(200,400,600,800,1000,1200)
}
if (is.null(BETAF0_VALUES)) {
  BETAF0_VALUES <- if (phenotype == "binary") c(0) else c(0,1)
}

## ------------------------- environment / data ------------------------------
suppressMessages({
  library("data.table")
  if (phenotype == "conti" && METHOD == "diffFix") library("brglm")  # conti diffFix needs brglm
})
if (USE_LOAD_ALL) { suppressMessages(library("devtools")); load_all(DIFFDRIVER_DIR) } else library("diffdriver")
# Simulation + power-comparison functions live in this repo's code/ and call
# simulate_selection(); sourced AFTER the package so they take precedence.
if (is.null(CODE_DIR)) {
  .argfile <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  CODE_DIR <- if (length(.argfile)) dirname(sub("^--file=", "", .argfile[1])) else "code"
}
source(file.path(CODE_DIR, "simulate_functions.R"))   # simulate_selection() + power_compare*()

load(file.path(DATA_DIR, "ri200.Rd"))       # -> ri200
load(file.path(DATA_DIR, "fanno200.Rd"))    # -> fanno200
load(file.path(DATA_DIR, if (geneset == "TSG") "beta_tsg.Rd" else "beta_og.Rd"))  # -> beta_gc: TRUE simulation effect sizes
load(file.path(DATA_DIR, "BMR.rda"))        # -> BMR
load(file.path(DATA_DIR, "hmm.rda"))        # -> hmm

# beta_gcFix = DiffDriver's built-in ("fixed") parameters used only by the diffFix
# method when RUNNING DiffDriver -- a property of the fitted model, distinct from
# the simulation coefficients beta_gc above. Derived from parmASHmean
# (parmASHmean[[1]] = TSG, [[2]] = OG), exactly as in the per-gene scripts.
beta_gcFix <- beta_gc
if (METHOD == "diffFix") {
  pf <- file.path(DATA_DIR, "parmASHmean.rda")
  if (!file.exists(pf)) {
    stop("parmASHmean.rda not found in '", DATA_DIR,
         "' -- required for the diffFix method's beta_gcFix.")
  }
  load(pf)                                  # -> parmASHmean
  beta_gcFix <- unlist(parmASHmean[[ash_idx]])
  names(beta_gcFix)[6] <- "(Intercept)"
  beta_gcFix <- beta_gcFix[names(beta_gc)]
}

## ------------------------- run: loop over all genes -------------------------
out_file <- if (is.null(OUT_FILE)) {
  file.path(OUT_DIR, sprintf("pppower_all_%s_%s.Rd", SETTING, METHOD))
} else OUT_FILE
message(sprintf("[%s / %s] genes=1:%d  betaf0=%s  sizes=%s  Nite=%d  hot=%d  -> %s",
                SETTING, METHOD, N_GENES, paste(BETAF0_VALUES, collapse=","),
                paste(SAMPLE_SIZES, collapse=","), Nite, hot, out_file))

# Every gene's result is collected in one list and saved once, at the end, to a
# single file. Each element is one gene's power_compare* output, named by the stem
# the old per-gene files used (pppower_betaf0=..._sample..._<gene>[...]) so results
# still map 1:1 to those file names.
simures_all <- list()

for (betaf0 in BETAF0_VALUES) {
  for (Nsample in SAMPLE_SIZES) {
    for (gene in 1:N_GENES) {
      signal <- gene <= SIGNAL_GENES
      if (phenotype == "binary") {
        para <- if (signal) c(para_var, 1 - para_var) else c(0.5, 0.5)
      } else {
        # continuous: c(mean, sd, intercept, slope); slope for signal genes is
        # the setting's conti_slope (default 2), 0 for null genes.
        slope <- if (is.null(S$conti_slope)) 2 else S$conti_slope
        para <- c(0, 2, 0, if (signal) slope else 0)
      }

      set.seed(10)
      simures <- if (METHOD == "diffFix") {
        run_selection_diffdriver(
          binary = (phenotype == "binary"), Niter = Nite,
          sganno = ri200[[gene]], sgmatrix = fanno200[[gene]], bmrpars = log(BMR),
          Nsample = Nsample, betaf0 = betaf0, beta_gc = beta_gc, beta_gcFix = beta_gcFix,
          para = para, hot = hot, hmm = hmm)
      } else {
        run_selection_benchmarks(
          binary = (phenotype == "binary"), Niter = Nite,
          sganno = ri200[[gene]], sgmatrix = fanno200[[gene]], bmrpars = log(BMR),
          Nsample = Nsample, betaf0 = betaf0, beta_gc = beta_gc,
          para = para, hot = hot, hmm = hmm)
      }

      key <- sub("\\.Rd$", "", S$fname(betaf0, Nsample, gene, para_var))
      simures_all[[key]] <- simures
    }
    message(sprintf("  done betaf0=%s sample=%s (%d genes)", betaf0, Nsample, N_GENES))
  }
}

save(simures_all, file = out_file)
message(sprintf("Saved %d results to %s", length(simures_all), out_file))
