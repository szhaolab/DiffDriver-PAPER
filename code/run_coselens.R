#' Patch dndscv output to add dummy indel columns when none exist.
#' Coselens assumes indel columns are always present in sel_cv.
patch_dndscv_indels <- function(dndsout) {
  if (!"n_ind" %in% colnames(dndsout$sel_cv)) {
    dndsout$sel_cv$n_ind <- 0
    dndsout$sel_cv$exp_ind <- 1e-10
    dndsout$sel_cv$wind_cv <- 0
    dndsout$sel_cv$pind_cv <- 1
    dndsout$sel_cv$qind_cv <- 1
    dndsout$sel_cv$pglobal_cv <- dndsout$sel_cv$pallsubs_cv
    dndsout$sel_cv$qglobal_cv <- dndsout$sel_cv$qallsubs_cv
  }
  if (is.null(dndsout$nbregind)) {
    dndsout$nbregind <- list(theta = 1)
  }
  dndsout
}

#' Convert simulated mutation matrices to MAF-format data.frame.
#'
#' @param mutlist List of 9 sparse mutation matrices (sites x samples).
#' @param sganno  List of 9 annotation data.tables with chrom/start/ref/alt.
#' @param phenotype Numeric phenotype vector (1=case, 0=control).
#' @return List with group1 (cases) and group2 (controls) data.frames,
#'   each with columns: sampleID, chr, pos, ref, alt.
sim_to_maf <- function(mutlist, sganno, phenotype) {
  all_muts <- list()
  for (t in seq_along(mutlist)) {
    m <- mutlist[[t]]
    if (sum(m) == 0) next
    idx <- Matrix::summary(m)
    anno <- sganno[[t]]
    all_muts[[length(all_muts) + 1]] <- data.frame(
      sampleID = paste0("S", idx$j),
      chr      = sub("^chr", "", anno$chrom[idx$i]),
      pos      = anno$start[idx$i],
      ref      = anno$ref[idx$i],
      alt      = anno$alt[idx$i],
      stringsAsFactors = FALSE
    )
  }
  if (length(all_muts) == 0) return(list(group1 = NULL, group2 = NULL))
  all_muts_df <- do.call(rbind, all_muts)

  case_ids <- paste0("S", which(phenotype == 1))
  ctrl_ids <- paste0("S", which(phenotype == 0))
  list(
    group1 = all_muts_df[all_muts_df$sampleID %in% case_ids, ],
    group2 = all_muts_df[all_muts_df$sampleID %in% ctrl_ids, ]
  )
}

#' Build the true per-context (192) substitution rate from the simulation's BMR.
#'
#' The simulation's neutral per-site mutation probability factorizes as
#' exp(bmrpars[nttype]) * exp(betaf0) * bmrfold[sample], with no position
#' dependence within a nttype. dndscv works in 192 trinucleotide contexts, so
#' we map each context to its nttype (data/context_nttype_map.rds, aligned to
#' the substmodel row order) and return the per-context base rate
#' exp(bmrpars[nttype] + betaf0). Multiply this by a group's summed bmrfold to
#' get that group's total per-context exposure (what dndscv's `mutrates` encodes).
#'
#' @param bmrpars Length-9 log background rates (= log(BMR)), indexed by nttype.
#' @param betaf0 Baseline functional effect (log scale).
#' @param map_path Path to the context->nttype map rds.
#' @return Named numeric vector of length 192 in substmodel context order.
build_true_mutrates <- function(bmrpars, betaf0,
                                map_path = "data/context_nttype_map.rds") {
  stopifnot(length(bmrpars) == 9)
  map <- readRDS(map_path)
  stopifnot(nrow(map) == 192, all(map$nttype %in% 1:9))
  rates <- exp(as.numeric(bmrpars)[map$nttype] + betaf0)
  setNames(rates, map$context)
}

#' Run coselens on simulated data with patched indel handling.
#'
#' @param group1 Data.frame of case mutations (sampleID, chr, pos, ref, alt).
#' @param group2 Data.frame of control mutations.
#' @param subset.genes.by Character vector of gene names to subset results.
#' @param true_rate Optional list to replace dndscv's estimated background with
#'   the simulation's known rate. Two ways to supply it:
#'   (a) SELECTION sim (9 nttypes): `bmrpars` (log(BMR), length 9), `betaf0`
#'       (scalar), and `bmrfold_sum1`/`bmrfold_sum2` (sum of per-sample bmrfold
#'       over the samples in group1/group2). The per-context 192 rate is built
#'       as build_true_mutrates(bmrpars, betaf0) * bmrfold_sum.
#'   (b) CONFOUNDING sim (96 contexts): `mutrates1`/`mutrates2`, already-built
#'       length-192 per-context rate vectors in substmodel order (one per group),
#'       used directly. Build them from the sim's bmrfold$bmr with
#'       data/context96_map.rds: rate96 = rowSums(bmr[, group]); the 192 vector
#'       is rate96[context96_map$number].
#'   Optional `theta` (default 1e6) sets the near-Poisson overdispersion used to
#'   force mrfold=1. When supplied, both groups are run through `dndscv_truerate`
#'   instead of coselens::dndscv; source("code/dndscv_truerate.R") first.
#' @param ... Additional arguments passed to dndscv.
#' @return coselens result list, or NULL on error.
run_coselens <- function(group1, group2, subset.genes.by = NULL,
                         true_rate = NULL, ...) {
  requireNamespace("coselens", quietly = TRUE)

  dnds_args <- list(
    refdb = "hg19", sm = "192r_3w", kc = "cgc81", cv = "hg19",
    max_muts_per_gene_per_sample = Inf,
    max_coding_muts_per_sample = 3000,
    use_indel_sites = FALSE,
    maxcovs = 20, constrain_wnon_wspl = TRUE,
    outp = 3, numcode = 1, mingenecovs = 500,
    outmats = TRUE, outmutrates = TRUE
  )
  dnds_args <- modifyList(dnds_args, list(...))

  if (is.null(true_rate)) {
    message("[coselens] Running dndscv for group 1")
    dnds1 <- do.call(coselens::dndscv, c(list(mutations = group1), dnds_args))
    message("[coselens] Running dndscv for group 2")
    dnds2 <- do.call(coselens::dndscv, c(list(mutations = group2), dnds_args))
  } else {
    if (!exists("dndscv_truerate"))
      stop("true_rate supplied but dndscv_truerate not found; ",
           "source('code/dndscv_truerate.R') first.")
    if (!is.null(true_rate$mutrates1)) {
      # (b) confounding sim: per-group 192-context rates supplied directly
      mr1 <- true_rate$mutrates1
      mr2 <- true_rate$mutrates2
      # covariates are unused when the rate is given; a restricted gene_list
      # would not match the genome-wide cv matrix, so force cv=NULL. The
      # list()-assignment keeps an explicit NULL (modifyList would drop it).
      dnds_args["cv"] <- list(NULL)
    } else {
      # (a) selection sim: build from log(BMR)/betaf0 and per-group bmrfold sum
      base <- build_true_mutrates(true_rate$bmrpars, true_rate$betaf0)
      mr1 <- base * true_rate$bmrfold_sum1
      mr2 <- base * true_rate$bmrfold_sum2
    }
    theta <- if (is.null(true_rate$theta)) 1e6 else true_rate$theta
    message("[coselens] Running dndscv_truerate (true BMR) for group 1")
    dnds1 <- do.call(dndscv_truerate, c(
      list(mutations = group1, true_mutrates = mr1, true_theta = theta), dnds_args))
    message("[coselens] Running dndscv_truerate (true BMR) for group 2")
    dnds2 <- do.call(dndscv_truerate, c(
      list(mutations = group2, true_mutrates = mr2, true_theta = theta), dnds_args))
  }

  dnds1 <- patch_dndscv_indels(dnds1)
  dnds2 <- patch_dndscv_indels(dnds2)

  calc_ex     <- coselens:::calc_ex
  lik_func    <- coselens:::lik_func
  get_effect_size <- coselens::get_effect_size

  group1_ex <- calc_ex(dnds1)
  group2_ex <- calc_ex(dnds2)
  group1_ex$ex_tot <- group1_ex[, 2] + group1_ex[, 3]
  group2_ex$ex_tot <- group2_ex[, 2] + group2_ex[, 3]
  ex_values <- merge(
    group1_ex[, c("gene_name", "ex_mis", "ex_non", "ex_tot")],
    group2_ex[, c("gene_name", "ex_mis", "ex_non", "ex_tot")],
    by = "gene_name", suffixes = c(".group1", ".group2")
  )

  num_patients_1 <- length(unique(dnds1$annotmuts$sampleID))
  num_patients_2 <- length(unique(dnds2$annotmuts$sampleID))
  n_genes <- nrow(dnds1$sel_cv)
  ex_ind_1 <- data.frame(gene_name = dnds1$sel_cv$gene_name,
                          ex_ind = (dnds1$sel_cv$n_ind - dnds1$sel_cv$exp_ind) / num_patients_1)
  ex_ind_2 <- data.frame(gene_name = dnds2$sel_cv$gene_name,
                          ex_ind = (dnds2$sel_cv$n_ind - dnds2$sel_cv$exp_ind) / num_patients_2)
  ex_ind_values <- merge(ex_ind_1, ex_ind_2, by = "gene_name",
                         suffixes = c(".group1", ".group2"))
  ex_values <- merge(ex_values, ex_ind_values, by = "gene_name")

  group1_sel_cv <- dnds1$sel_cv[order(as.numeric(rownames(dnds1$sel_cv))), , drop = FALSE]
  group2_sel_cv <- dnds2$sel_cv[order(as.numeric(rownames(dnds2$sel_cv))), , drop = FALSE]

  h0_sel_cv <- lik_func(dnds1, dnds2)

  h0_cols <- c("gene_name", "llall", "llmis", "lltrunc",
               "llmis_check_A", "lltrunc_check_B", "llind0", "llind1")
  g1_cols <- c("gene_name", "wmis_cv", "wnon_cv", "llall")

  lldf <- merge(h0_sel_cv[, h0_cols], group1_sel_cv[, g1_cols],
                by = "gene_name", suffixes = c(".0", ".group1"))
  lldf <- merge(lldf, group2_sel_cv[, g1_cols],
                by = "gene_name", suffixes = c(".group1", ".group2"))
  nc <- ncol(lldf)
  colnames(lldf)[(nc - 2):nc] <- c("wmis_cv.group2", "wnon_cv.group2", "llall.group2")

  lldf$llall.1 <- lldf$llall.group1 + lldf$llall.group2
  lldf <- lldf[, c(1, ncol(lldf), 2:(ncol(lldf) - 1))]

  if (length(subset.genes.by) > 0) {
    lldf <- lldf[lldf$gene_name %in% subset.genes.by, ]
    ex_values <- ex_values[ex_values$gene_name %in% subset.genes.by, ]
  }

  lldf$pall    <- pchisq(-2 * (lldf$llall.0 - lldf$llall.1), df = 2, lower.tail = FALSE)
  lldf$pmis    <- pchisq(-2 * (lldf$llmis - lldf$llall.1), df = 1, lower.tail = FALSE)
  lldf$ptrunc  <- pchisq(-2 * (lldf$lltrunc - lldf$llall.1), df = 1, lower.tail = FALSE)
  lldf$pind    <- 1 - pchisq(2 * (lldf$llind1 - lldf$llind0), df = 1)
  lldf$pglobal <- 1 - pchisq(-2 * (log(lldf$pall) + log(lldf$pind)), df = 4)

  lldf <- lldf[order(lldf$pglobal), ]
  lldf$qind    <- p.adjust(lldf$pind, method = "BH")
  lldf$qall    <- p.adjust(lldf$pall, method = "BH")
  lldf$qmis    <- p.adjust(lldf$pmis, method = "BH")
  lldf$qtrunc  <- p.adjust(lldf$ptrunc, method = "BH")
  lldf$qglobal <- p.adjust(lldf$pglobal, method = "BH")

  lldf <- merge(lldf, ex_values[, c("gene_name", "ex_tot.group1", "ex_tot.group2",
                                     "ex_ind.group1", "ex_ind.group2",
                                     "ex_mis.group1", "ex_mis.group2",
                                     "ex_non.group1", "ex_non.group2")],
                by = "gene_name")

  single.test.names <- c("psub.group", "pind.group", "pmis.group", "ptrunc.group",
                          "pglobal.group", "qsub.group", "qind.group", "qmis.group",
                          "qtrunc.group", "qglobal.group")
  single_cols <- c("gene_name", "pallsubs_cv", "pind_cv", "pmis_cv", "ptrunc_cv",
                    "pglobal_cv", "qallsubs_cv", "qind_cv", "qmis_cv",
                    "qtrunc_cv", "qglobal_cv")
  lldf <- merge(lldf, group1_sel_cv[, single_cols], by = "gene_name")
  colnames(lldf)[(ncol(lldf) - length(single.test.names) + 1):ncol(lldf)] <-
    paste0(single.test.names, "1")

  lldf <- merge(lldf, group2_sel_cv[, single_cols], by = "gene_name")
  colnames(lldf)[(ncol(lldf) - length(single.test.names) + 1):ncol(lldf)] <-
    paste0(single.test.names, "2")

  full.data.from.lldf <- c("gene_name", "ex_tot.group1", "ex_tot.group2",
                            "ex_ind.group1", "ex_ind.group2",
                            "pall", "pind", "pglobal", "qall", "qind", "qglobal",
                            "ex_mis.group1", "ex_mis.group2",
                            "ex_non.group1", "ex_non.group2",
                            "pmis", "ptrunc", "qmis", "qtrunc",
                            paste0(single.test.names, "1"),
                            paste0(single.test.names, "2"))
  full.data.from.lldf.names <- c("gene_name", "num.driver.sub.group1",
                                  "num.driver.sub.group2", "num.driver.ind.group1",
                                  "num.driver.ind.group2",
                                  "psub", "pind", "pglobal", "qsub", "qind", "qglobal",
                                  "num.driver.mis.group1", "num.driver.mis.group2",
                                  "num.driver.trunc.group1", "num.driver.trunc.group2",
                                  "pmis", "ptrunc", "qmis", "qtrunc",
                                  paste0(single.test.names, "1"),
                                  paste0(single.test.names, "2"))
  full.data <- lldf[, full.data.from.lldf]
  colnames(full.data) <- full.data.from.lldf.names

  out.list <- list()
  out.list[["substitutions"]] <- get_effect_size(full.data, mutation.class = "sub")
  out.list[["indels"]]        <- get_effect_size(full.data, mutation.class = "ind")
  out.list[["missense_sub"]]  <- get_effect_size(full.data, mutation.class = "mis")
  out.list[["truncating_sub"]]<- get_effect_size(full.data, mutation.class = "trunc")
  out.list[["overall_mut"]]   <- get_effect_size(full.data, mutation.class = "global")
  out.list[["dndscv"]]        <- list(dndscv_group1 = dnds1, dndscv_group2 = dnds2)

  message("[coselens] Done.")
  out.list
}
