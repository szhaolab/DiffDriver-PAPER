---
name: coselens-benchmark
description: How coselens is benchmarked against DiffDriver in the simulation, incl. true-BMR injection and the shared-cohort requirement
metadata:
  type: project
---

coselens (github ggruenhagen3/coselens; dN/dS via dndscv) is benchmarked against DiffDriver on the DiffDriver-PAPER simulations.

Code: `code/run_coselens.R` (`sim_to_maf`, `run_coselens`, `build_true_mutrates`, `patch_dndscv_indels`), `code/dndscv_truerate.R` (edited copy of `coselens::dndscv` adding `true_mutrates`/`true_theta`), `data/context_nttype_map.rds` (192 dndscv trinucleotide contexts -> sim's 9 nttypes).

Two non-obvious constraints:
- **Indels**: the sim produces only substitutions; stock coselens crashes (calc_ex_ind / lik_func expect indel columns). `patch_dndscv_indels` + `use_indel_sites=FALSE` fixes it.
- **Shared cohort**: coselens runs genome-wide on ONE cohort, but the sim draws an independent phenotype per gene. A valid benchmark must simulate all genes against a SINGLE phenotype vector (test harness forces this via a per-gene `set.seed`).

Simulation function (June 2026): `simulate_selection` was replaced by `simulate_1funcvi` (ported from szhaolab/diffdriver scripts/simulate_functions.R). Key change: the neutral background rate per site is `exp(bmrpars[t]+betaf0)`, **constant across samples** (no per-sample bmrfold, so the old `rho`/`tau` args are gone). `run_simulation.R` callers updated accordingly. Because BMR is sample-constant, `true_rate$bmrfold_sum1/2` is just the number of samples assigned to each phenotype group (cases / controls), counting zero-mutation samples too — e.g. 600 each for a 1200-sample 50/50 split.

**Why:** DiffDriver is given the true BMR (`mr` arg to `ddmodel`); to compare fairly the user wanted coselens/dndscv to also use the true rate instead of estimating it. **How to apply:** pass `true_rate=list(bmrpars=log(BMR), betaf0=, bmrfold_sum1=sum(pheno==1), bmrfold_sum2=sum(pheno==0))`; this injects exp_* = colSums(L*mutrates_true), skips the NB regression, and forces mrfold=1 via near-Poisson theta (1e6). Verified: injected exp_* equal truth; p-values theta-stable to ~1e-6. Related: [[simulation-pipeline]].
