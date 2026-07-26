# Rebuild data/pppower_bmr_variations.Rd partA from the estimated-background variation
# worker output (simulate_bmr_variations.R -> var_<VAR>_N<N>.rds), keeping partB
# (signature-count robustness, simulate_signature_count.R) as-is.
SRC <- "/Volumes/Szhao/jiezhou/diffdata/randombmrsingle/tsgsimple_LUAD_variations_est"
VARS  <- c("baseline","tmb_cv0.5","tmb_cv0.9","R2_0.3","R2_0.5","R2_0.95","continuous")
SIZES <- c(200,400,600,800,1000,1200,1500)

load("data/pppower_bmr_variations.Rd")   # existing `bmrvar` (partA, partB, sizes)
stopifnot(identical(bmrvar$sizes, SIZES))

partA <- setNames(vector("list", length(VARS)), VARS)
missing <- character(0)
for (v in VARS) {
  lst <- setNames(vector("list", length(SIZES)), as.character(SIZES))
  for (N in SIZES) {
    f <- file.path(SRC, sprintf("var_%s_N%d.rds", v, N))
    if (file.exists(f)) lst[[as.character(N)]] <- readRDS(f)
    else missing <- c(missing, sprintf("%s_N%d", v, N))
  }
  partA[[v]] <- lst
}
if (length(missing)) { cat("MISSING:\n"); print(missing); stop("incomplete results") }

bmrvar$partA <- partA            # replace partA with estimated-rate results; keep partB
save(bmrvar, file = "data/pppower_bmr_variations.Rd")

fp <- function(p) round(mean(p < 0.01, na.rm = TRUE), 2)
for (v in VARS) {
  cat(sprintf("\n%-11s ", v))
  for (N in SIZES) {
    d <- partA[[v]][[as.character(N)]]
    cat(sprintf("N%d[Cos=%.2f DD=%.2f Lin=%.2f] ", N, fp(d$Coselens), fp(d$DiffDriver), fp(d$Linear)))
  }
}
cat("\n")
