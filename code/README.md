# Code

Command-line scripts and shared R code for the main-text analyses.

## Simulation / benchmarking (Fig. 3)
- `run_simulation.R`, `simulate_all_genes.R`, `simulate_functions*.R` — power/FPR simulations.
- `run_coselens.R` — Coselens baseline for the method comparison.
- `consolidate_results.R`, `dndscv_truerate.R` — collate simulation outputs.

## Pan-cancer immune subtypes (Fig. 5a, `pancancer299`)
- `run_cohort_multi.R` — runs DiffDriver on the 299 driver-gene cohorts (Bailey et al.),
  testing each gene for differential selection between each immune subtype (Thorsson
  C1–C4, one-vs-rest) across 8,253 tumors / 27 TCGA cancers, adjusting for age, sex,
  and tumor type (sex/tumor-type dropped for single-level cohorts).
- `recover_cohort.R` — reruns cohorts that failed in the first pass.
- `aggregate_pancancer299.R` — combines per-cohort C1–C4 results into
  `data/pancancer299_immune_aggregated.txt`; BH FDR is computed within each immune
  subtype (`p.adjust(dd.p, "fdr"), by = phenotype`). Yields 18/220 genes (8.2%).
- `job_pancan299.sh`, `job_recover.sh` — cluster submission scripts.

The clinical per-trait analysis (Fig. 4a, 45/185 = 24.3%) is reproduced directly in
`analysis/Clinical.Rmd` from `data/clinical_pertrait_pvalues.tsv` (BH across tumor
types within each trait).
