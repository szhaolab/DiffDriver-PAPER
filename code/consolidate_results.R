#!/usr/bin/env Rscript
# =============================================================================
# consolidate_results.R
#
# Reads the existing per-gene result files (one .Rd per gene, in each result
# folder's diffFix/ or other/ subdir) for ONE setting x method and combines them
# into a single file `pppower_all_<SETTING>_<METHOD>.Rd` containing a list
# `simures_all` keyed by the per-gene filename stem -- exactly the format
# produced by simulate_all_genes.R.
#
# Usage: Rscript consolidate_results.R <SETTING> <METHOD> [OUT_DIR] [DATA_DIR]
#   SETTING: TSG_binary_0.8 | TSG_binary_0.9 | TSG_conti |
#            OG_binary_0.8  | OG_binary_0.9  | OG_conti
#   METHOD : diffFix | other
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: consolidate_results.R <SETTING> <METHOD> [OUT_DIR] [DATA_DIR]")
SETTING <- args[1]
METHOD  <- args[2]
OUT_DIR <- if (length(args) >= 3) args[3] else "."
DATA_DIR<- if (length(args) >= 4) args[4] else "/dartfs/rc/lab/S/Szhao/jiezhou/diffdata"

BINARY_SIZES <- c(200, 400, 600, 800, 1000, 1200, 1500)
CONTI_SIZES  <- c(200, 400, 600, 800, 1000, 1200)
N_GENES      <- 200

# base = result folder; sizes / betaf0 / para variant / filename convention must
# match the Rmd's load_and_compute and simulate_all_genes.R's presets.
presets <- list(
  TSG_binary_0.8 = list(base = file.path(DATA_DIR, "gene239/conti/TSG_BRCA"),
    sizes = BINARY_SIZES, betaf0 = 0, pv = 0.8,
    fname = function(b, N, g, pv) sprintf("pppower_betaf0=%s_sample%s_%s.Rd", b, N, g)),
  TSG_binary_0.9 = list(base = file.path(DATA_DIR, "gene239/conti/TSG_BRCA/TSG0.9/TSG_gene"),
    sizes = BINARY_SIZES, betaf0 = 0, pv = 0.9,
    fname = function(b, N, g, pv) sprintf("pppower_betaf0=%s_sample%s_%s.Rd", b, N, g)),
  TSG_conti      = list(base = file.path(DATA_DIR, "gene239/conti/TSG_BRCA/conti/conti2"),
    sizes = CONTI_SIZES, betaf0 = c(0, 1), pv = NA,
    fname = function(b, N, g, pv) sprintf("pppower_betaf0=%s_sample%s_%s.Rd", b, N, g)),
  OG_binary_0.8  = list(base = file.path(DATA_DIR, "binaryEM/OG_gene"),
    sizes = BINARY_SIZES, betaf0 = 0, pv = 0.8,
    fname = function(b, N, g, pv) sprintf("pppower_betaf0=%s_sample%s_%s_%s.Rd", b, N, pv, g)),
  OG_binary_0.9  = list(base = file.path(DATA_DIR, "binaryEM/OG_gene"),
    sizes = BINARY_SIZES, betaf0 = 0, pv = 0.9,
    fname = function(b, N, g, pv) sprintf("pppower_betaf0=%s_sample%s_%s_%s.Rd", b, N, pv, g)),
  OG_conti       = list(base = file.path(DATA_DIR, "gene239/conti/OCG/OG_gene"),
    sizes = CONTI_SIZES, betaf0 = c(0, 1), pv = NA,
    fname = function(b, N, g, pv) sprintf("pppower_betaf0=%s_sample%s_%s.Rd", b, N, g))
)
if (!SETTING %in% names(presets)) stop("Unknown SETTING: ", SETTING)
if (!METHOD %in% c("diffFix", "other")) stop("METHOD must be diffFix or other")
S       <- presets[[SETTING]]
objname <- if (METHOD == "diffFix") "simuresdiffFix" else "simuresother"

message(sprintf("[%s / %s] reading from %s", SETTING, METHOD, file.path(S$base, METHOD)))

simures_all <- list()
n_ok <- 0L; missing <- character(0)
for (b in S$betaf0) {
  for (N in S$sizes) {
    for (g in seq_len(N_GENES)) {
      fn  <- S$fname(b, N, g, S$pv)
      fp  <- file.path(S$base, METHOD, fn)
      key <- sub("\\.Rd$", "", fn)
      if (!file.exists(fp)) { missing <- c(missing, fn); next }
      e <- new.env(); load(fp, envir = e)
      simures_all[[key]] <- get(objname, envir = e)
      n_ok <- n_ok + 1L
    }
  }
}

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
out <- file.path(OUT_DIR, sprintf("pppower_all_%s_%s.Rd", SETTING, METHOD))
save(simures_all, file = out)
message(sprintf("Saved %d results (%d missing) to %s", n_ok, length(missing), out))
if (length(missing)) message("  first missing: ", paste(head(missing, 5), collapse = ", "))
