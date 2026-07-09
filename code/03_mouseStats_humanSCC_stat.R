## Reviewer comment 5: statistics for Figure 1 driver-variant comparisons
##
## 1. Mouse tumor burden:
##    Kruskal-Wallis test across four driver variants, followed by pairwise
##    Wilcoxon rank-sum tests with multiple-testing correction.
##
## 2. Human SSCs:
##    Evaluate feasibility of likelihood-ratio tests (LRTs). SSCs are
##    model-derived parameters from cancereffectsizeR, not replicate-level
##    observations, so ANOVA/Kruskal-Wallis on SSC point estimates is not
##    appropriate.

suppressPackageStartupMessages({
  library(cancereffectsizeR)
  library(data.table)
})

location_data <- "../data/"
rdata_output <- "./R_Data/"
stats_output <- "./SuppleStatistics/"

if (!dir.exists(stats_output)) {
  dir.create(stats_output, recursive = TRUE)
}

driver_variants <- c("KRAS_G12D", "BRAF_V600E", "KRAS_G12C", "EGFR_L858R")
driver_labels <- c(
  "KRAS_G12D" = "KRAS G12D",
  "BRAF_V600E" = "BRAF V600E",
  "KRAS_G12C" = "KRAS G12C",
  "EGFR_L858R" = "EGFR L858R"
)

## Mouse tumor burden ---------------------------------------------------------

mouse_tumor <- fread(paste0(
  location_data,
  "fig3f_sgInertTumorBurden_webPlotDigitizer_final.csv"
))
mouse_tumor[, Genotype := fifelse(
  Genotype == "G12D", "KRAS_G12D",
  fifelse(
    Genotype == "G12C", "KRAS_G12C",
    fifelse(Genotype == "EGFR", "EGFR_L858R",
            fifelse(Genotype == "BRAF", "BRAF_V600E", Genotype))
  )
)]
mouse_tumor <- mouse_tumor[Genotype %in% driver_variants]
mouse_tumor[, Genotype := factor(Genotype, levels = driver_variants)]

mouse_summary <- mouse_tumor[, .(
  n = .N,
  median_TM_M = median(TM_M, na.rm = TRUE),
  mean_TM_M = mean(TM_M, na.rm = TRUE),
  q025_TM_M = as.numeric(quantile(TM_M, 0.025, na.rm = TRUE)),
  q975_TM_M = as.numeric(quantile(TM_M, 0.975, na.rm = TRUE)),
  min_TM_M = min(TM_M, na.rm = TRUE),
  max_TM_M = max(TM_M, na.rm = TRUE)
), by = Genotype]
mouse_summary[, Genotype_label := driver_labels[as.character(Genotype)]]
setcolorder(mouse_summary, c("Genotype", "Genotype_label"))

kw_mouse <- kruskal.test(TM_M ~ Genotype, data = mouse_tumor)
kw_mouse_table <- data.table(
  test = "Kruskal-Wallis rank-sum test",
  statistic = unname(kw_mouse$statistic),
  df = unname(kw_mouse$parameter),
  p_value = kw_mouse$p.value,
  method = kw_mouse$method
)

mouse_pairwise_raw <- rbindlist(lapply(
  combn(driver_variants, 2, simplify = FALSE),
  function(pair) {
    test <- wilcox.test(
      TM_M ~ Genotype,
      data = mouse_tumor[Genotype %in% pair],
      exact = FALSE
    )
    data.table(
      variant_1 = pair[1],
      variant_2 = pair[2],
      statistic = unname(test$statistic),
      p_value = test$p.value
    )
  }
))
mouse_pairwise_raw[, p_adj_holm := p.adjust(p_value, method = "holm")]

pairwise_mouse_table <- copy(mouse_pairwise_raw)
pairwise_mouse_table[, `:=`(
  variant_1_label = driver_labels[as.character(variant_1)],
  variant_2_label = driver_labels[as.character(variant_2)],
  p_adjust_method = "Holm",
  test = "Pairwise Wilcoxon rank-sum test"
)]
setcolorder(
  pairwise_mouse_table,
  c(
    "test", "variant_1", "variant_1_label", "variant_2",
    "variant_2_label", "statistic", "p_value", "p_adj_holm",
    "p_adjust_method"
  )
)

fwrite(mouse_summary, paste0(stats_output, "mouse_tumor_burden_summary.csv"))
fwrite(kw_mouse_table, paste0(stats_output, "mouse_tumor_burden_kw_test.csv"))
fwrite(
  pairwise_mouse_table,
  paste0(stats_output, "mouse_tumor_burden_pairwise_wilcoxon_holm.csv")
)

## Human SSC LRT feasibility --------------------------------------------------

extract_driver_ssc <- function(cesa_file, group_label) {
  cesa <- readRDS(cesa_file)
  selection_results <- slot(cesa, "selection_results")
  recurrent <- selection_results[[1]]
  out <- recurrent[variant_name %in% driver_variants]
  out[, group := group_label]
  setcolorder(out, "group")
  out[]
}

human_ssc_summary <- rbindlist(list(
  extract_driver_ssc(paste0(rdata_output, "cesa_nonsmoking.rds"), "Never-smoker"),
  extract_driver_ssc(paste0(rdata_output, "cesa_smoking.rds"), "Smoker")
), use.names = TRUE, fill = TRUE)
human_ssc_summary[, variant_label := driver_labels[variant_name]]
setcolorder(human_ssc_summary, c("group", "variant_name", "variant_label"))
fwrite(human_ssc_summary, paste0(stats_output, "human_driver_ssc_summary.csv"))

## The functions below reconstruct the one-parameter SSWM likelihood used by
## ces_variant(model = "default") and fit a constrained shared-gamma null model.
## They are kept in this script so that feasibility can be assessed by checking
## whether the reconstructed per-variant alternative log-likelihoods match the
## saved cancereffectsizeR loglikelihood values.

get_model_samples <- function(cesa) {
  samples <- slot(cesa, "samples")
  samples[
    !is.na(sig_analysis_grp) & !is.na(gene_rate_grp),
    Unique_Patient_Identifier
  ]
}

get_variant_likelihood_inputs <- function(cesa, variant_rows, model_samples) {
  samples <- slot(cesa, "samples")[
    Unique_Patient_Identifier %in% model_samples
  ]
  maf <- slot(cesa, "maf")[
    Unique_Patient_Identifier %in% model_samples
  ]
  mutations <- slot(cesa, "mutations")

  tumors_with_variants_by_gene <- unique(
    maf[, .(gene = unlist(genes)), by = "Unique_Patient_Identifier"]
  )[, .(samples = list(Unique_Patient_Identifier)), by = "gene"]
  tumors_with_variants_by_gene <- list2env(
    setNames(
      tumors_with_variants_by_gene$samples,
      tumors_with_variants_by_gene$gene
    ),
    parent = emptyenv()
  )

  snv_aac_of_interest <- mutations$aac_snv_key[
    variant_rows$variant_id,
    on = "aac_id"
  ]
  samples_by_aac <- maf[
    snv_aac_of_interest,
    .(variant_id),
    on = c(variant_id = "snv_id"),
    by = "Unique_Patient_Identifier",
    nomatch = NULL
  ]
  samples_by_aac[
    snv_aac_of_interest,
    aac_id := aac_id,
    on = c(variant_id = "snv_id")
  ]
  samples_by_aac <- samples_by_aac[
    ,
    .(samples = list(unique(Unique_Patient_Identifier))),
    by = "aac_id"
  ]
  samples_by_aac <- list2env(
    setNames(samples_by_aac$samples, samples_by_aac$aac_id),
    parent = emptyenv()
  )

  setkey(samples, "covered_regions")
  genome_wide_cov_samples <- samples[
    "genome",
    Unique_Patient_Identifier,
    nomatch = NULL
  ]

  lapply(seq_len(nrow(variant_rows)), function(i) {
    variant_row <- variant_rows[i]
    coverage_group <- mutations$amino_acid_change[
      variant_row$variant_id,
      covered_in,
      on = "aac_id"
    ][[1]]

    if (length(coverage_group) == 1 && is.na(coverage_group)) {
      covered_samples <- genome_wide_cov_samples
    } else {
      covered_samples <- c(
        samples[coverage_group, Unique_Patient_Identifier, nomatch = NULL],
        genome_wide_cov_samples
      )
    }
    covered_samples <- unique(covered_samples)

    baseline_rates <- baseline_mutation_rates(
      cesa = cesa,
      aac_ids = variant_row$variant_id,
      samples = covered_samples
    )
    rates <- baseline_rates[[variant_row$variant_id]]
    names(rates) <- baseline_rates$Unique_Patient_Identifier

    tumors_with_variant <- intersect(
      samples_by_aac[[variant_row$variant_id]],
      covered_samples
    )
    tumors_with_gene_mutated <- tumors_with_variants_by_gene[[variant_row$gene]]
    tumors_without <- setdiff(covered_samples, tumors_with_gene_mutated)

    list(
      group = variant_row$group,
      variant_name = variant_row$variant_name,
      variant_id = variant_row$variant_id,
      saved_loglikelihood = variant_row$loglikelihood,
      rates_tumors_with = rates[tumors_with_variant],
      rates_tumors_without = rates[tumors_without],
      included_with_variant = length(tumors_with_variant),
      included_total = length(tumors_with_variant) + length(tumors_without)
    )
  })
}

sswm_minus_loglik <- function(gamma, rates_tumors_with, rates_tumors_without) {
  if (!is.finite(gamma) || gamma <= 0) {
    return(Inf)
  }
  loglik <- 0
  if (length(rates_tumors_without) > 0) {
    loglik <- loglik - sum(gamma * rates_tumors_without)
  }
  if (length(rates_tumors_with) > 0) {
    loglik <- loglik + sum(log1p(-exp(-gamma * rates_tumors_with)))
  }
  -loglik
}

fit_one_gamma <- function(likelihood_inputs) {
  opt <- optimize(
    f = function(gamma) {
      sum(vapply(
        likelihood_inputs,
        function(x) {
          sswm_minus_loglik(
            gamma,
            x$rates_tumors_with,
            x$rates_tumors_without
          )
        },
        numeric(1)
      ))
    },
    interval = c(0.001, 1e9)
  )
  list(gamma = opt$minimum, loglikelihood = -opt$objective)
}

fit_each_gamma <- function(likelihood_inputs) {
  rbindlist(lapply(likelihood_inputs, function(x) {
    opt <- optimize(
      f = function(gamma) {
        sswm_minus_loglik(gamma, x$rates_tumors_with, x$rates_tumors_without)
      },
      interval = c(0.001, 1e9)
    )
    data.table(
      group = x$group,
      variant_name = x$variant_name,
      variant_id = x$variant_id,
      fitted_gamma = opt$minimum,
      reconstructed_loglikelihood = -opt$objective,
      saved_loglikelihood = x$saved_loglikelihood,
      loglikelihood_delta = (-opt$objective) - x$saved_loglikelihood,
      reconstructed_included_with_variant = x$included_with_variant,
      reconstructed_included_total = x$included_total
    )
  }))
}

run_omnibus_lrt <- function(cesa_file, group_label) {
  cesa <- readRDS(cesa_file)
  model_samples <- get_model_samples(cesa)
  variant_rows <- extract_driver_ssc(cesa_file, group_label)
  inputs <- get_variant_likelihood_inputs(cesa, variant_rows, model_samples)

  alternative_fit <- fit_each_gamma(inputs)
  null_fit <- fit_one_gamma(inputs)
  alt_loglik <- sum(alternative_fit$reconstructed_loglikelihood)
  lrt_stat <- 2 * (alt_loglik - null_fit$loglikelihood)

  lrt <- data.table(
    group = group_label,
    test = "Omnibus LRT: one shared SSC vs four separate SSCs",
    null_df = 1L,
    alternative_df = length(driver_variants),
    lrt_df = length(driver_variants) - 1L,
    null_gamma = null_fit$gamma,
    null_loglikelihood = null_fit$loglikelihood,
    alternative_loglikelihood = alt_loglik,
    lrt_statistic = lrt_stat,
    p_value = pchisq(lrt_stat, df = length(driver_variants) - 1L, lower.tail = FALSE)
  )

  list(lrt = lrt, alternative_diagnostic = alternative_fit)
}

run_pairwise_lrt <- function(cesa_file, group_label) {
  cesa <- readRDS(cesa_file)
  model_samples <- get_model_samples(cesa)
  variant_rows <- extract_driver_ssc(cesa_file, group_label)
  inputs_all <- get_variant_likelihood_inputs(cesa, variant_rows, model_samples)
  names(inputs_all) <- vapply(inputs_all, `[[`, character(1), "variant_name")

  pair_results <- rbindlist(lapply(combn(driver_variants, 2, simplify = FALSE), function(pair) {
    inputs <- inputs_all[pair]
    alternative_fit <- fit_each_gamma(inputs)
    null_fit <- fit_one_gamma(inputs)
    alt_loglik <- sum(alternative_fit$reconstructed_loglikelihood)
    lrt_stat <- 2 * (alt_loglik - null_fit$loglikelihood)
    data.table(
      group = group_label,
      variant_1 = pair[1],
      variant_2 = pair[2],
      variant_1_label = driver_labels[pair[1]],
      variant_2_label = driver_labels[pair[2]],
      test = "Pairwise LRT: one shared SSC vs two separate SSCs",
      lrt_df = 1L,
      null_gamma = null_fit$gamma,
      null_loglikelihood = null_fit$loglikelihood,
      alternative_loglikelihood = alt_loglik,
      lrt_statistic = lrt_stat,
      p_value = pchisq(lrt_stat, df = 1, lower.tail = FALSE)
    )
  }))
  pair_results[, p_adj_holm := p.adjust(p_value, method = "holm"), by = group]
  pair_results[]
}

## This is intentionally TRUE: it runs quickly for only four driver variants
## per human group and provides a concrete feasibility diagnostic.
run_lrt_prototype <- TRUE

human_lrt_feasibility <- data.table(
  item = c(
    "Appropriate test family",
    "Existing alternative model",
    "Null model required for LRT",
    "Implementation difficulty",
    "Main validation check",
    "Interpretation caution"
  ),
  assessment = c(
    "Model-based likelihood-ratio tests are appropriate; ANOVA/Kruskal-Wallis on SSC point estimates is not.",
    "Available: ces_variant() saved per-variant selection_intensity and loglikelihood for each driver in smoker and never-smoker groups.",
    "Needed: constrained fit with a shared selection_intensity across selected drivers, using the same per-sample baseline rates and same gene hold-out rule.",
    "Moderate: requires reconstructing baseline rates and variant carrier/eligible sample sets from the CESAnalysis object; no package-level LRT wrapper is exposed.",
    "Reconstructed separate-SSC log-likelihoods should closely match the saved ces_variant() loglikelihood values before using the LRT P values.",
    "Pairwise LRTs among KRAS_G12C and KRAS_G12D share a gene and use gene-level hold-outs; report as model-based comparisons, not replicate-level rank tests."
  )
)

if (run_lrt_prototype) {
  lrt_ns <- run_omnibus_lrt(
    paste0(rdata_output, "cesa_nonsmoking.rds"),
    "Never-smoker"
  )
  lrt_sm <- run_omnibus_lrt(
    paste0(rdata_output, "cesa_smoking.rds"),
    "Smoker"
  )
  human_omnibus_lrt <- rbindlist(list(lrt_ns$lrt, lrt_sm$lrt))
  human_lrt_diagnostic <- rbindlist(list(
    lrt_ns$alternative_diagnostic,
    lrt_sm$alternative_diagnostic
  ))
  human_pairwise_lrt <- rbindlist(list(
    run_pairwise_lrt(paste0(rdata_output, "cesa_nonsmoking.rds"), "Never-smoker"),
    run_pairwise_lrt(paste0(rdata_output, "cesa_smoking.rds"), "Smoker")
  ))

  max_abs_loglik_delta <- max(
    abs(human_lrt_diagnostic$loglikelihood_delta),
    na.rm = TRUE
  )
  human_lrt_feasibility <- rbind(
    human_lrt_feasibility,
    data.table(
      item = "Prototype LRT status",
      assessment = paste0(
        "Prototype ran. Maximum absolute difference between reconstructed ",
        "and saved alternative log-likelihoods: ",
        signif(max_abs_loglik_delta, 4),
        ". Inspect human_lrt_reconstruction_diagnostic.csv before reporting."
      )
    )
  )

  fwrite(
    human_omnibus_lrt,
    paste0(stats_output, "human_driver_ssc_omnibus_lrt_prototype.csv")
  )
  fwrite(
    human_pairwise_lrt,
    paste0(stats_output, "human_driver_ssc_pairwise_lrt_prototype.csv")
  )
  fwrite(
    human_lrt_diagnostic,
    paste0(stats_output, "human_lrt_reconstruction_diagnostic.csv")
  )
}

fwrite(
  human_lrt_feasibility,
  paste0(stats_output, "human_driver_ssc_lrt_feasibility.csv")
)

## Unified supplementary table ------------------------------------------------
## The table reports inferential results once. Its legend explains that the
## mouse tests apply to both the Figure 1 mean display and the Supplementary
## Figure 1 median display, which use the same replicate-level observations.

mouse_overall_table <- data.table(
  Group = "Mouse tumor burden",
  Comparison = "Overall among four variants",
  Test = "Kruskal-Wallis",
  `P value` = kw_mouse$p.value,
  `Adjusted P value` = NA_real_,
  Conclusion = ifelse(
    kw_mouse$p.value < 0.05,
    "Significant",
    "Not significant"
  )
)

mouse_pairwise_table_for_supp <- pairwise_mouse_table[, .(
  Group = "Mouse tumor burden",
  Comparison = paste(variant_1_label, "vs", variant_2_label),
  Test = "Wilcoxon rank-sum",
  `P value` = p_value,
  `Adjusted P value` = p_adj_holm,
  Conclusion = ifelse(
    p_adj_holm < 0.05,
    "Significant",
    "Not significant"
  )
)]

supplementary_table <- rbindlist(
  list(mouse_overall_table, mouse_pairwise_table_for_supp),
  use.names = TRUE,
  fill = TRUE
)

if (run_lrt_prototype) {
  human_omnibus_table <- human_omnibus_lrt[, .(
    Group = paste("Human SSC,", group),
    Comparison = "Overall among four driver variants",
    Test = "LRT",
    `P value` = p_value,
    `Adjusted P value` = NA_real_,
    Conclusion = ifelse(
      p_value < 0.05,
      "Significant",
      "Not significant"
    )
  )]

  human_pairwise_table_for_supp <- human_pairwise_lrt[, .(
    Group = paste("Human SSC,", group),
    Comparison = paste(variant_1_label, "vs", variant_2_label),
    Test = "Pairwise LRT",
    `P value` = p_value,
    `Adjusted P value` = p_adj_holm,
    Conclusion = ifelse(
      p_adj_holm < 0.05,
      "Significant",
      "Not significant after Holm adjustment"
    )
  )]

  supplementary_table <- rbindlist(
    list(
      supplementary_table,
      human_omnibus_table,
      human_pairwise_table_for_supp
    ),
    use.names = TRUE,
    fill = TRUE
  )
}

supplementary_table[, group_order := match(
  Group,
  c("Mouse tumor burden", "Human SSC, Never-smoker", "Human SSC, Smoker")
)]
supplementary_table[, comparison_order := fifelse(
  startsWith(Comparison, "Overall among four"),
  0L,
  1L
)]
setorder(supplementary_table, group_order, comparison_order, Comparison)
supplementary_table[, c("group_order", "comparison_order") := NULL]

fwrite(
  supplementary_table,
  paste0(
    stats_output,
    "Supplementary_Table_driver_variant_statistical_comparisons.csv"
  )
)


