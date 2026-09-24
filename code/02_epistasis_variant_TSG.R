library(ggplot2)
library(ggpubr)
library(patchwork)
library(data.table)
library(cancereffectsizeR)
library(ces.refset.hg19)

## epistasis analysis between a variant and a gene; using select_variants() to include more samples #####
get_epistasis_cesa <- function(
  cesa_object,
  variant,
  gene_name,
  specific_samples
) {
  variants_geneName <- cancereffectsizeR::select_variants(
    cesa = cesa_object,
    genes = gene_name
  )
  variants_variant <- cancereffectsizeR::select_variants(
    cesa_object,
    variant_ids = variant
  )
  variants_target = rbind(variants_variant, variants_geneName)
  comp_target <- define_compound_variants(
    cesa = cesa_object,
    variant_table = variants_target,
    by = "gene",
    merge_distance = Inf
  )

  run_name <- paste0(variant, "_vs_", gene_name)
  cesa_object <- ces_epistasis(
    cesa = cesa_object,
    variants = comp_target,
    run_name = run_name,
    samples = specific_samples
  )
  epi <- cesa_object$epistasis[[run_name]]
  set(epi, j = "run_type", value = run_name)
  result <- list(epi, cesa_object, comp_target, variants_target)
  return(result)
}

tsg <- c(
  "PTEN",
  "STK11",
  "Nf1",
  "Setd2",
  "Stag2",
  "Tsc1",
  "Rb1",
  "Fbxw7",
  "Rnf43",
  "Rbm10",
  "Mga",
  "Keap1",
  "Kmt2d",
  "Smad4",
  "Cmtr2",
  "Tp53",
  "Nf2",
  "Apc",
  "Ptprd",
  "Atm",
  "Msh2",
  "Arid2",
  "Kdm6a",
  "Atrx",
  "BRCA2",
  "Palb2"
)
tsg <- toupper(tsg)
tsg <- c(tsg, "CDKN2A.p14arf", "CDKN2A.p16INK4a")

# Pre-allocate epi_final with column names from the first iteration
KRAS_G12D_vs_PTEN <- get_epistasis_cesa(
  cesa_object = cesa_smoking_w_panel,
  variant = oncogene_variants[1],
  gene_name = tsg[1],
  specific_samples = sample_smo_wP_forCesa
)
epi_initial <- KRAS_G12D_vs_PTEN[[1]]
epi_final <- data.table(matrix(NA, nrow = 0, ncol = length(epi_initial))) ## create an empty data table with 0 rows and length(epi_initial) columns
setnames(epi_final, names(epi_initial)) ## then set the column names
epi_final_smoking_wP <- epi_final
epi_final_nonsmoking_wP <- epi_final

# for smoker and Never smokers
for (i in 1:length(oncogene_variants)) {
  for (j in 1:length(tsg)) {
    epi_smoking_wP <- get_epistasis_cesa(
      cesa_object = cesa_smoking_w_panel,
      variant = oncogene_variants[i],
      gene_name = tsg[j],
      specific_samples = sample_smo_wP_forCesa
    )[[1]]
    epi_final_smoking_wP <- rbindlist(
      list(epi_final_smoking_wP, epi_smoking_wP),
      fill = F
    )

    epi_nonsmoking_wP <- get_epistasis_cesa(
      cesa_object = cesa_nonsmoking_w_panel,
      variant = oncogene_variants[i],
      gene_name = tsg[j],
      specific_samples = sample_nonsmo_wP_forCesa
    )[[1]]
    epi_final_nonsmoking_wP <- rbindlist(
      list(epi_final_nonsmoking_wP, epi_nonsmoking_wP),
      fill = F
    )
  }
}

## for all samples, needed for checking sample sizes for each pair: nAB ---
epi_final_pan <- epi_final
for (i in 1:length(oncogene_variants)) {
  for (j in 1:length(tsg)) {
    epi_pan <- get_epistasis_cesa(
      cesa_object = cesa,
      variant = oncogene_variants[i],
      gene_name = tsg[j],
      specific_samples = sample_pan_forCsea
    )[[1]]
    epi_final_pan <- rbindlist(
      list(epi_final_pan, epi_pan),
      fill = F
    )
  }
}
epi_final_pan[, data_type := "Pan"]
epi_final_pan[nAB == "0"]

epi_data <- rbind(
  epi_final_smoking_wP[, data_type := "Smo"],
  epi_final_nonsmoking_wP[, data_type := "nonSmo"]
)
epi_data <- epi_data[,
  c("variant_A_type", "remaining_column") := tstrsplit(
    run_type,
    split = "_vs",
    keep = 1
  )
]
epi_data[, gene_Name := stringr::str_split(run_type, "_", simplify = T)[, 4]]

## remove all rows with "CDKN2A.p16INK4a" in "gene_Name" column and "CDKN2A.p14arf" in "variant_B" column
## because select_variants() function also select all variant information of CDKN2A.p14arf" when the "gene_name" is "CDKN2A.p16INK4a", which cause a duplication.
epi_data <- epi_data[
  !.("CDKN2A.p16INK4a", "CDKN2A.p14arf"),
  on = .(gene_Name, variant_B)
]

## add variant_B_s column, to shorten CDKN2A.p16INK4a and CDKN2A.p14arf to p16INK4a and p14arf
epi_data[, variant_B_s := variant_B]
epi_data[
  stringr::str_detect(variant_B_s, "\\."),
  variant_B_s := stringr::str_split(variant_B_s, "\\.", simplify = T)[, 2]
]

save(
  epi_data,
  file = paste0(rdata_output, "epistasis_variant_TSG_output.Rdata")
)


### for KRAS_G12C ####
epi_KRAS_G12C <- epi_data["KRAS_G12C", on = "variant_A_type"]
epi_KRAS_G12C_S <- epi_data[
  .("KRAS_G12C", "Smo"),
  on = .(variant_A_type, data_type)
]
epi_KRAS_G12C_nonS <- epi_data[
  .("KRAS_G12C", "nonSmo"),
  on = .(variant_A_type, data_type)
]

## sorting gene by ces_B_on_A of TSG under KRSA_G12C background in NS, then smokers
variant_order_KRAS_G12C <- epi_KRAS_G12C_nonS[ces_B_on_A > 0.001000][order(
  -ces_B_on_A
)]$variant_B_s
variant_lf_KRAS_G12C <- epi_KRAS_G12C_nonS[ces_B_on_A == 0.001000]$variant_B_s
variant_lf_order_KRAS_G12C <- epi_KRAS_G12C_S[
  variant_lf_KRAS_G12C,
  on = "variant_B_s"
][ces_B_on_A > 0.001000][order(-ces_B_on_A)]$variant_B_s
variant_lf_2_KRAS_G12C <- epi_KRAS_G12C_S[
  variant_lf_KRAS_G12C,
  on = "variant_B_s"
][ces_B_on_A == 0.001000]$variant_B_s

gene_order_KRAS_G12C <- c(
  variant_order_KRAS_G12C,
  variant_lf_order_KRAS_G12C,
  variant_lf_2_KRAS_G12C
)
gene_order_KRAS_G12C

epi_KRAS_G12C[,
  variant_B_s := factor(variant_B_s, levels = gene_order_KRAS_G12C, ordered = T)
]
epi_KRAS_G12C[,
  data_type := factor(data_type, levels = c("nonSmo", "Smo"), ordered = T)
]

## edit CI bar which has NA values,CI should be zero if it is NA
epi_KRAS_G12C$ci_low_95_ces_B_on_A_forplot <- ifelse(
  is.na(epi_KRAS_G12C$ci_low_95_ces_B_on_A),
  0,
  epi_KRAS_G12C$ci_low_95_ces_B_on_A
)
epi_KRAS_G12C$errorBarColor <- ifelse(
  epi_KRAS_G12C$nAB %in% c("0", "1"),
  "Gray",
  "Black"
)
epi_KRAS_G12C$data_type2 <- as.vector(epi_KRAS_G12C$data_type)
epi_KRAS_G12C$data_type2[epi_KRAS_G12C$nAB == "0"] <- "Gray"


p_G12C <- ggplot(data = epi_KRAS_G12C, aes(x = data_type, y = ces_B_on_A)) +
  geom_errorbar(
    aes(
      ymin = ci_low_95_ces_B_on_A_forplot,
      ymax = ci_high_95_ces_B_on_A,
      color = errorBarColor
    ),
    width = 0,
    na.rm = T
  ) +
  geom_point(
    aes(fill = data_type2, color = errorBarColor),
    shape = 21,
    size = 2
  ) +
  scale_fill_manual(
    values = c("Smo" = "red", "nonSmo" = "blue", "Gray" = "gray")
  ) + ## "All" = "black",
  scale_color_manual(values = c("Gray" = "gray", "Black" = "black")) +
  labs(x = "Tumor suppressor gene") +
  ylab(bquote(bold(atop("Cancer effect", "(KRAS G12C background)")))) +
  scale_y_log10(
    labels = function(x) format(x, big.mark = ",", scientific = F)
  ) +
  scale_x_discrete(labels = c("Never-smoker", "Smoker")) + ## change the x-axis label ## "All",
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
    axis.text.y = element_text(size = 8, face = "bold"),
    axis.title.y = element_text(size = 8, face = "bold"),
    axis.title.x = element_text(size = 8, face = "bold"),
    legend.position = "none"
  ) +
  facet_wrap(~variant_B_s, nrow = 1) +
  theme(
    strip.text.x = element_text(
      angle = 45,
      hjust = 0.5,
      vjust = 0.5,
      size = 8,
      face = "bold"
    )
  )

### for KRAS_G12D ####
epi_KRAS_G12D <- epi_data["KRAS_G12D", on = "variant_A_type"]
epi_KRAS_G12D_S <- epi_data[
  .("KRAS_G12D", "Smo"),
  on = .(variant_A_type, data_type)
]
epi_KRAS_G12D_nonS <- epi_data[
  .("KRAS_G12D", "nonSmo"),
  on = .(variant_A_type, data_type)
]

## use the gene_order_KRAS_G12C
epi_KRAS_G12D[,
  variant_B_s := factor(variant_B_s, levels = gene_order_KRAS_G12C, ordered = T)
]
epi_KRAS_G12D[,
  data_type := factor(data_type, levels = c("nonSmo", "Smo"), ordered = T)
]

## edit CI bar which has NA values, CI should be zero if it is NA
epi_KRAS_G12D$ci_low_95_ces_B_on_A_forplot <- ifelse(
  is.na(epi_KRAS_G12D$ci_low_95_ces_B_on_A),
  0,
  epi_KRAS_G12D$ci_low_95_ces_B_on_A
)
epi_KRAS_G12D$errorBarColor <- ifelse(
  epi_KRAS_G12D$nAB %in% c("0", "1"),
  "Gray",
  "Black"
)
epi_KRAS_G12D$data_type2 <- as.vector(epi_KRAS_G12D$data_type)
epi_KRAS_G12D$data_type2[epi_KRAS_G12D$nAB == "0"] <- "Gray"


p_G12D <- ggplot(data = epi_KRAS_G12D, aes(x = data_type, y = ces_B_on_A)) +
  geom_errorbar(
    aes(
      ymin = ci_low_95_ces_B_on_A_forplot,
      ymax = ci_high_95_ces_B_on_A,
      color = errorBarColor
    ),
    width = 0,
    na.rm = T
  ) +
  geom_point(
    aes(fill = data_type2, color = errorBarColor),
    shape = 21,
    size = 2
  ) +
  scale_fill_manual(
    values = c("Smo" = "red", "nonSmo" = "blue", "Gray" = "gray")
  ) +
  scale_color_manual(values = c("Gray" = "gray", "Black" = "black")) +
  labs(x = "Tumor suppressor gene") +
  ylab(bquote(bold(atop("Cancer effect", "(KRAS G12D background)")))) +
  scale_y_log10(
    labels = function(x) format(x, big.mark = ",", scientific = F)
  ) +
  scale_x_discrete(labels = c("Never-smoker", "Smoker")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
    axis.text.y = element_text(size = 8, face = "bold"),
    axis.title.y = element_text(size = 8, face = "bold"),
    axis.title.x = element_text(size = 8, face = "bold"),
    legend.position = "none"
  ) +
  facet_wrap(~variant_B_s, nrow = 1) +
  theme(
    strip.text.x = element_text(
      angle = 45,
      hjust = 0.5,
      vjust = 0.5,
      size = 8,
      face = "bold"
    )
  ) #+

### for BRAF_V600E ####
epi_BRAF_V600E <- epi_data["BRAF_V600E", on = "variant_A_type"]
epi_BRAF_V600E_S <- epi_data[
  .("BRAF_V600E", "Smo"),
  on = .(variant_A_type, data_type)
]
epi_BRAF_V600E_nonS <- epi_data[
  .("BRAF_V600E", "nonSmo"),
  on = .(variant_A_type, data_type)
]

## use the gene_order_KRAS_G12C
epi_BRAF_V600E[,
  variant_B_s := factor(variant_B_s, levels = gene_order_KRAS_G12C, ordered = T)
]
epi_BRAF_V600E[,
  data_type := factor(data_type, levels = c("nonSmo", "Smo"), ordered = T)
]

## edit CI bar which has NA values,Jeff think that CI should be zero if it is NA
epi_BRAF_V600E$ci_low_95_ces_B_on_A_forplot <- ifelse(
  is.na(epi_BRAF_V600E$ci_low_95_ces_B_on_A),
  0,
  epi_BRAF_V600E$ci_low_95_ces_B_on_A
)
epi_BRAF_V600E$errorBarColor <- ifelse(
  epi_BRAF_V600E$nAB %in% c("0", "1"),
  "Gray",
  "Black"
)
epi_BRAF_V600E$data_type2 <- as.vector(epi_BRAF_V600E$data_type)
epi_BRAF_V600E$data_type2[epi_BRAF_V600E$nAB == "0"] <- "Gray"


p_V600E <- ggplot(data = epi_BRAF_V600E, aes(x = data_type, y = ces_B_on_A)) +
  geom_errorbar(
    aes(
      ymin = ci_low_95_ces_B_on_A_forplot,
      ymax = ci_high_95_ces_B_on_A,
      color = errorBarColor
    ),
    width = 0,
    na.rm = T
  ) +
  geom_point(
    aes(fill = data_type2, color = errorBarColor),
    shape = 21,
    size = 2
  ) +
  scale_fill_manual(
    values = c("Smo" = "red", "nonSmo" = "blue", "Gray" = "gray")
  ) +
  scale_color_manual(values = c("Gray" = "gray", "Black" = "black")) +
  labs(x = "Tumor suppressor gene") +
  ylab(bquote(bold(atop("Cancer effect", "(BRAF V600E background)")))) +
  scale_y_log10(
    labels = function(x) format(x, big.mark = ",", scientific = F)
  ) +
  scale_x_discrete(labels = c("Never-smoker", "Smoker")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
    axis.text.y = element_text(size = 8, face = "bold"),
    axis.title.y = element_text(size = 8, face = "bold"),
    axis.title.x = element_text(size = 8, face = "bold"),
    legend.position = "none"
  ) +
  facet_wrap(~variant_B_s, nrow = 1) +
  theme(
    strip.text.x = element_text(
      angle = 45,
      hjust = 0.5,
      vjust = 0.5,
      size = 8,
      face = "bold"
    )
  ) #+

### for EGFR_L858R ####
epi_EGFR_L858R <- epi_data["EGFR_L858R", on = "variant_A_type"]
epi_EGFR_L858R_S <- epi_data[
  .("EGFR_L858R", "Smo"),
  on = .(variant_A_type, data_type)
]
epi_EGFR_L858R_nonS <- epi_data[
  .("EGFR_L858R", "nonSmo"),
  on = .(variant_A_type, data_type)
]
## use the gene_order_KRAS_G12C
epi_EGFR_L858R[,
  variant_B_s := factor(variant_B_s, levels = gene_order_KRAS_G12C, ordered = T)
]
epi_EGFR_L858R[,
  data_type := factor(data_type, levels = c("nonSmo", "Smo"), ordered = T)
]

## edit CI bar which has NA values, Jeff think that CI should be zero if it is NA
epi_EGFR_L858R$ci_low_95_ces_B_on_A_forplot <- ifelse(
  is.na(epi_EGFR_L858R$ci_low_95_ces_B_on_A),
  0,
  epi_EGFR_L858R$ci_low_95_ces_B_on_A
) ## epi_EGFR_L858R$ces_B_on_A
epi_EGFR_L858R$errorBarColor <- ifelse(
  epi_EGFR_L858R$nAB %in% c("0", "1"),
  "Gray",
  "Black"
)
epi_EGFR_L858R$data_type2 <- as.vector(epi_EGFR_L858R$data_type)
epi_EGFR_L858R$data_type2[epi_EGFR_L858R$nAB == "0"] <- "Gray"


p_L858R <- ggplot(data = epi_EGFR_L858R, aes(x = data_type, y = ces_B_on_A)) +
  geom_errorbar(
    aes(
      ymin = ci_low_95_ces_B_on_A_forplot,
      ymax = ci_high_95_ces_B_on_A,
      color = errorBarColor
    ),
    width = 0,
    na.rm = T
  ) +
  geom_point(
    aes(fill = data_type2, color = errorBarColor),
    shape = 21,
    size = 2
  ) +
  scale_fill_manual(
    values = c("Smo" = "red", "nonSmo" = "blue", "Gray" = "gray")
  ) +
  scale_color_manual(values = c("Gray" = "gray", "Black" = "black")) +
  labs(x = "Tumor suppressor gene") +
  ylab(bquote(bold(atop("Cancer effect", "(EGFR L858R background)")))) +
  scale_y_log10(
    labels = function(x) format(x, big.mark = ",", scientific = F)
  ) +
  scale_x_discrete(labels = c("Never-smoker", "Smoker")) + #
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
    axis.text.y = element_text(size = 8, face = "bold"),
    axis.title.y = element_text(size = 8, face = "bold"),
    axis.title.x = element_text(size = 8, face = "bold"),
    legend.position = "none"
  ) +
  facet_wrap(~variant_B_s, nrow = 1) +
  theme(
    strip.text.x = element_text(
      angle = 45,
      hjust = 0.5,
      vjust = 0.5,
      size = 8,
      face = "bold"
    )
  )


### plot CES_B_on_A of TSG under a specific variant background in one figure; only NS, S; sort by NS; NA is 0    ####
jpeg(
  file = "ces_B_on_A_variant_TSG_NAis0.jpeg",
  width = 13,
  height = 13,
  units = 'in',
  res = 300
)
p_G12C +
  p_G12D +
  p_V600E +
  p_L858R +
  plot_annotation(tag_levels = 'A') +
  plot_layout(nrow = 4, guides = "collect")
dev.off()




# scatter plot: compare relative tumor size/number with CES_B_on_A-----
format_p_math <- function(p) {
  if (p >= 0.001) {
    sprintf("%.3f", p)
  } else {
    s <- formatC(p, format = "e", digits = 3)
    parts <- strsplit(s, "e")[[1]]
    mantissa <- sprintf("%.3f", as.numeric(parts[1]))
    exponent <- as.integer(parts[2])
    bquote(.(mantissa) %*% 10^.(exponent))
  }
}

draw_scatter <- function(
  epi_data,
  MiceTumorData,
  MiceTumorData_type = "size", ## "size" or "number"
  variant_onco,
  MiceTumorlimit,
  data_type,
  max_overlaps,
  label_y,
  label_x
) {
  ## MiceTumorlimit = c(0.5,32)
  library(ggplot2)
  variants <- deparse(substitute(variant_onco)) # Get the name of the variable passed as variant_onco
  variants <- stringr::str_replace(variants, "_", " ")
  epi_TS_data <- dplyr::inner_join(
    epi_data,
    MiceTumorData,
    by = c("variant_B" = "TSG")
  )

  pearson_corr <- cor.test(
    epi_TS_data$ces_B_on_A,
    epi_TS_data$RelativeTumor_95,
    method = c("pearson")
  ) ## can only assign one method each time
  pearson_r <- pearson_corr$estimate
  pearson_p <- pearson_corr$p.value

  spearman_corr <- cor.test(
    epi_TS_data$ces_B_on_A,
    epi_TS_data$RelativeTumor_95,
    method = c("spearman")
  ) ## can only assign one method each time
  spearman_r <- spearman_corr$estimate
  spearman_p <- spearman_corr$p.value

  if (data_type == "S") {
    data_type_label <- "smokers"
  } else if (data_type == "nonS") {
    data_type_label <- "never-smokers"
  }

  pp <- ggplot(
    data = epi_TS_data[order(-ces_B_on_A)],
    aes(x = ces_B_on_A, y = RelativeTumor_95)
  ) + ##reorder(variant_B, log10(ces_B_on_A))
    geom_point(aes(color = variant_B)) + ## position = "jitter", size = 1
    labs(
      y = paste0(
        "Relative tumor",
        MiceTumorData_type,
        " in ",
        stringr::str_split(variants, " ", simplify = T)[2],
        ";Cas9"
      )
    ) +
    xlab(bquote(
      .(data_type_label) * ";" * ~CES ~ of ~ TSG ~ on ~ .(variants) ~
        scriptstyle(~ ~ (log[10]))
    )) +
    scale_x_log10(
      labels = function(x) format(x, big.mark = ",", scientific = F)
    ) +
    scale_y_continuous(trans = "log2", limits = MiceTumorlimit) + #tgutil::scale_y_log2()+
    geom_errorbar(
      data = epi_TS_data[
        !is.na(epi_TS_data$CI_low_tumor) &
          !is.na(epi_TS_data$CI_high_tumor),
      ],
      aes(ymin = CI_low_tumor, ymax = CI_high_tumor),
      width = 0,
      color = "black"
    ) +
    # geom_errorbar(data = epi_TS_data[!is.na(epi_TS_data$ci_low_95_ces_B_on_A) & !is.na(epi_TS_data$ci_high_95_ces_B_on_A), ],
    #               aes(xmin = ci_low_95_ces_B_on_A, xmax = ci_high_95_ces_B_on_A), width = 0, color = "black")+
    geom_smooth(method = "lm", se = T) +
    ggrepel::geom_label_repel(
      aes(label = variant_B),
      max.overlaps = max_overlaps,
      color = "grey"
    ) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(
      title = paste(
        expression(rho),
        "=",
        round(spearman_r, 2),
        ";",
        "p =",
        format(spearman_p, scientific = TRUE, digits = 2),
        "\n",
        "r",
        "=",
        round(pearson_r, 2),
        ";",
        "p =",
        format(pearson_p, scientific = TRUE, digits = 2)
      )
    )
  pp_wo_smooth <- ggplot(
    data = epi_TS_data[order(-ces_B_on_A)],
    aes(x = ces_B_on_A, y = RelativeTumor_95)
  ) + ##reorder(variant_B, log10(ces_B_on_A))

    labs(
      y = paste0(
        "Relative tumor ",
        MiceTumorData_type,
        # " in ",
        " (",
        stringr::str_split(variants, " ", simplify = T)[2],
        ";Cas9)"
      )
    ) +
    xlab(bquote(
      "TSG selection in" ~ .(variants) ~ .(data_type_label)
      # "(" * .(data_type_label) * ")" #~
      #scriptstyle(~ (log[10]))
    )) +
    scale_x_log10(
      labels = function(x) format(x, big.mark = ",", scientific = F)
    ) +
    scale_y_continuous(trans = "log2", limits = MiceTumorlimit) + #tgutil::scale_y_log2()+
    geom_errorbar(
      data = epi_TS_data[
        !is.na(epi_TS_data$CI_low_tumor) &
          !is.na(epi_TS_data$CI_high_tumor),
      ],
      aes(ymin = CI_low_tumor, ymax = CI_high_tumor),
      width = 0,
      color = "grey"
    ) +
    geom_errorbar(
      data = epi_TS_data[
        !is.na(epi_TS_data$ci_low_95_ces_B_on_A) &
          !is.na(epi_TS_data$ci_high_95_ces_B_on_A),
      ],
      aes(xmin = ci_low_95_ces_B_on_A, xmax = ci_high_95_ces_B_on_A),
      width = 0,
      color = "grey"
    ) +
    geom_point(aes(fill = variant_B), shape = 21, color = "black") + ## position = "jitter", size = 1
    #geom_smooth(method = "lm", se = T)+
    ggrepel::geom_label_repel(
      aes(label = variant_B),
      max.overlaps = max_overlaps,
      color = "black", # text
      #fill = "grey",        # fill in the rectangle
      # segment.colour = "grey",# border + connecting line

      size = 4.5
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 18, hjust = 0.8), # , face = "bold"# angle = 45,hjust = 1
      axis.text.y = element_text(size = 18), # ,face = "bold"
      axis.title.y = element_text(size = 18), #, face = "bold"
      axis.title.x = element_text(size = 18)
    ) + # , face = "bold"
    
    annotate(
      "text",
      x = label_x,
      y = label_y,
      label = bquote(
        atop(
          "Spearman's " *
            rho *
            ": " *
            .(round(spearman_r, 3)) *
            "; " *
            italic(P) *
            " = " *
            .(format_p_math(spearman_p)),

          "Pearson's " *
            italic(r) *
            ": " *
            .(round(pearson_r, 3)) *
            "; " *
            italic(P) *
            " = " *
            .(format_p_math(pearson_p))
        )
      ),
      hjust = -0.1,
      vjust = 1,
      size = 5
    )

  return(
    result = list(
      pp,
      pp_wo_smooth,
      epi_TS_data,
      pearson_r,
      pearson_p,
      spearman_r,
      spearman_p
    )
  )
}

### relative tumor size ####
## KRAS_G12D; rm near 0 CES values
G12D_tumorSize <- fread(paste0(location_data, "G12D_data_tumorSize_final.csv"))

epi_KRAS_G12D_nonS_TS_ana_rm0 <- draw_scatter(
  epi_KRAS_G12D_nonS[ces_B_on_A >= 1],
  G12D_tumorSize,
  MiceTumorData_type = "size",
  variant_onco = KRAS_G12D,
  MiceTumorlimit = c(0.5, 32),
  data_type = "nonS",
  max_overlaps = 5,
  label_y = 32,
  label_x = 10
)
epi_KRAS_G12D_nonS_TS_ana_rm0[[2]]
epi_KRAS_G12D_S_TS_ana_rm0 <- draw_scatter(
  epi_KRAS_G12D_S[ces_B_on_A >= 1],
  G12D_tumorSize,
  MiceTumorData_type = "size",
  variant_onco = KRAS_G12D,
  MiceTumorlimit = c(0.5, 32),
  data_type = "S",
  max_overlaps = 5,
  label_y = 32,
  label_x = 0
)
epi_KRAS_G12D_S_TS_ana_rm0[[2]]

# (epi_KRAS_G12D_nonS_TS_ana_rm0[[2]] / epi_KRAS_G12D_S_TS_ana_rm0[[2]]) +
#   plot_annotation(tag_levels = 'A')
## KRAS_G12C; rm near 0 CES values
G12C_tumorSize <- fread(paste0(location_data, "G12C_data_tumorSize_final.csv"))
epi_KRAS_G12C_S_TS_ana_rm0 <- draw_scatter(
  epi_KRAS_G12C_S[ces_B_on_A >= 1],
  G12C_tumorSize,
  MiceTumorData_type = "size",
  variant_onco = KRAS_G12C,
  MiceTumorlimit = c(0.5, 32),
  data_type = "S",
  max_overlaps = 5,
  label_y = 32,
  label_x = 0
)
epi_KRAS_G12C_S_TS_ana_rm0[[2]]
epi_KRAS_G12C_nonS_TS_ana_rm0 <- draw_scatter(
  epi_KRAS_G12C_nonS[ces_B_on_A >= 1],
  G12C_tumorSize,
  MiceTumorData_type = "size",
  variant_onco = KRAS_G12C,
  MiceTumorlimit = c(0.5, 32),
  data_type = "nonS",
  max_overlaps = 5,
  label_y = 32,
  label_x = 0
)
epi_KRAS_G12C_nonS_TS_ana_rm0[[2]]
## BRAF_V600E; rm near 0 CES values
BRAF_V600E_tumorSize <- fread(paste0(
  location_data,
  "BRAF_V600E_data_tumorSize_final.csv"
))
epi_BRAF_V600E_S_TS_ana_rm0 <- draw_scatter(
  epi_BRAF_V600E_S[ces_B_on_A >= 1],
  BRAF_V600E_tumorSize,
  MiceTumorData_type = "size",
  variant_onco = BRAF_V600E,
  MiceTumorlimit = c(0.25, 8),
  data_type = "S",
  max_overlaps = 5,
  label_y = 8,
  label_x = 0
)
epi_BRAF_V600E_S_TS_ana_rm0[[2]]
epi_BRAF_V600E_nonS_TS_ana_rm0 <- draw_scatter(
  epi_BRAF_V600E_nonS[ces_B_on_A >= 1],
  BRAF_V600E_tumorSize,
  MiceTumorData_type = "size",
  variant_onco = BRAF_V600E,
  MiceTumorlimit = c(0.25, 8),
  data_type = "nonS",
  max_overlaps = 5,
  label_y = 8,
  label_x = 0
)
epi_BRAF_V600E_nonS_TS_ana_rm0[[2]]
## EGFR L858R; rm near 0 CES values
EGFR_L858R_tumorSize <- fread(paste0(
  location_data,
  "EGFR_L858R_data_tumorSize_final.csv"
))
epi_EGFR_L858R_S_TS_ana_rm0 <- draw_scatter(
  epi_EGFR_L858R_S[ces_B_on_A >= 1],
  EGFR_L858R_tumorSize,
  MiceTumorData_type = "size",
  variant_onco = EGFR_L858R,
  MiceTumorlimit = c(0.25, 8),
  data_type = "S",
  max_overlaps = 5,
  label_y = 8,
  label_x = 0
)
epi_EGFR_L858R_S_TS_ana_rm0[[2]]
epi_EGFR_L858R_nonS_TS_ana_rm0 <- draw_scatter(
  epi_EGFR_L858R_nonS[ces_B_on_A >= 1],
  EGFR_L858R_tumorSize,
  MiceTumorData_type = "size",
  variant_onco = EGFR_L858R,
  MiceTumorlimit = c(0.25, 8),
  data_type = "nonS",
  max_overlaps = 5,
  label_y = 8,
  label_x = 1
)
epi_EGFR_L858R_nonS_TS_ana_rm0[[2]]
svg(file = "p_MiceTumorSize_vs_epi.svg", height = 24, width = 16)
(epi_KRAS_G12C_nonS_TS_ana_rm0[[2]] | epi_KRAS_G12C_S_TS_ana_rm0[[2]]) /
  (epi_KRAS_G12D_nonS_TS_ana_rm0[[2]] | epi_KRAS_G12D_S_TS_ana_rm0[[2]]) /
  (epi_BRAF_V600E_nonS_TS_ana_rm0[[2]] | epi_BRAF_V600E_S_TS_ana_rm0[[2]]) /
  (epi_EGFR_L858R_nonS_TS_ana_rm0[[2]] | epi_EGFR_L858R_S_TS_ana_rm0[[2]]) /
  plot_annotation(tag_levels = waiver())
dev.off()

# svg(file = "p_MiceTumorSize_vs_epi_v2.svg", height = 16, width = 11)
jpeg(
  file = "p_MicetumorSize_vs_epi.jpeg",
  height = 24,
  width = 16,
  unit = "in",
  res = 600
)
(epi_KRAS_G12C_nonS_TS_ana_rm0[[2]] | epi_KRAS_G12C_S_TS_ana_rm0[[2]]) /
  (epi_KRAS_G12D_nonS_TS_ana_rm0[[2]] | epi_KRAS_G12D_S_TS_ana_rm0[[2]]) /
  (epi_BRAF_V600E_nonS_TS_ana_rm0[[2]] | epi_BRAF_V600E_S_TS_ana_rm0[[2]]) /
  (epi_EGFR_L858R_nonS_TS_ana_rm0[[2]] | epi_EGFR_L858R_S_TS_ana_rm0[[2]]) /
  plot_annotation(tag_levels = waiver())
dev.off()

### relative tumor number ####
## KRAS_G12D; rm near 0 CES values
G12D_tumorNumber <- fread(paste0(
  location_data,
  "G12D_data_tumorNumber_final.csv"
))

epi_KRAS_G12D_nonS_TS_ana_rm0_TN <- draw_scatter(
  epi_KRAS_G12D_nonS[ces_B_on_A >= 1],
  G12D_tumorNumber,
  MiceTumorData_type = "number",
  variant_onco = KRAS_G12D,
  MiceTumorlimit = c(0.12, 8),
  data_type = "nonS",
  max_overlaps = 5,
  label_y = 8,
  label_x = 10
)
epi_KRAS_G12D_nonS_TS_ana_rm0_TN[[2]]
epi_KRAS_G12D_S_TS_ana_rm0_TN <- draw_scatter(
  epi_KRAS_G12D_S[ces_B_on_A >= 1],
  G12D_tumorNumber,
  MiceTumorData_type = "number",
  variant_onco = KRAS_G12D,
  MiceTumorlimit = c(0.12, 8),
  data_type = "S",
  max_overlaps = 5,
  label_y = 8,
  label_x = 0
)
epi_KRAS_G12D_S_TS_ana_rm0_TN[[2]]

## KRAS_G12C; rm near 0 CES values
G12C_tumorNumber <- fread(paste0(
  location_data,
  "G12C_data_tumorNumber_final.csv"
))
epi_KRAS_G12C_S_TS_ana_rm0_TN <- draw_scatter(
  epi_KRAS_G12C_S[ces_B_on_A >= 1],
  G12C_tumorNumber,
  MiceTumorData_type = "number",
  variant_onco = KRAS_G12C,
  MiceTumorlimit = c(0.12, 8),
  data_type = "S",
  max_overlaps = 5,
  label_y = 8,
  label_x = 0
)
epi_KRAS_G12C_S_TS_ana_rm0_TN[[2]]
epi_KRAS_G12C_nonS_TS_ana_rm0_TN <- draw_scatter(
  epi_KRAS_G12C_nonS[ces_B_on_A >= 1],
  G12C_tumorNumber,
  MiceTumorData_type = "number",
  variant_onco = KRAS_G12C,
  MiceTumorlimit = c(0.12, 8),
  data_type = "nonS",
  max_overlaps = 5,
  label_y = 8,
  label_x = 0
)
epi_KRAS_G12C_nonS_TS_ana_rm0_TN[[2]]
## BRAF_V600E; rm near 0 CES values
BRAF_V600E_tumorNumber <- fread(paste0(
  location_data,
  "BRAF_V600E_data_tumorNumber_final.csv"
))
epi_BRAF_V600E_S_TS_ana_rm0_TN <- draw_scatter(
  epi_BRAF_V600E_S[ces_B_on_A >= 1],
  BRAF_V600E_tumorNumber,
  MiceTumorData_type = "number",
  variant_onco = BRAF_V600E,
  MiceTumorlimit = c(0.12, 8),
  data_type = "S",
  max_overlaps = 5,
  label_y = 8,
  label_x = 0
)
epi_BRAF_V600E_S_TS_ana_rm0_TN[[2]]
epi_BRAF_V600E_nonS_TS_ana_rm0_TN <- draw_scatter(
  epi_BRAF_V600E_nonS[ces_B_on_A >= 1],
  BRAF_V600E_tumorNumber,
  MiceTumorData_type = "number",
  variant_onco = BRAF_V600E,
  MiceTumorlimit = c(0.12, 8),
  data_type = "nonS",
  max_overlaps = 5,
  label_y = 8,
  label_x = 0
)
epi_BRAF_V600E_nonS_TS_ana_rm0_TN[[2]]
## EGFR L858R; rm near 0 CES values
EGFR_L858R_tumorNumber <- fread(paste0(
  location_data,
  "EGFR_L858R_data_tumorNumber_final.csv"
))
epi_EGFR_L858R_S_TS_ana_rm0_TN <- draw_scatter(
  epi_EGFR_L858R_S[ces_B_on_A >= 1],
  EGFR_L858R_tumorNumber,
  MiceTumorData_type = "number",
  variant_onco = EGFR_L858R,
  MiceTumorlimit = c(0.12, 8),
  data_type = "S",
  max_overlaps = 5,
  label_y = 8,
  label_x = 0
)
epi_EGFR_L858R_S_TS_ana_rm0_TN[[2]]
epi_EGFR_L858R_nonS_TS_ana_rm0_TN <- draw_scatter(
  epi_EGFR_L858R_nonS[ces_B_on_A >= 1],
  EGFR_L858R_tumorNumber,
  MiceTumorData_type = "number",
  variant_onco = EGFR_L858R,
  MiceTumorlimit = c(0.12, 8),
  data_type = "nonS",
  max_overlaps = 5,
  label_y = 8,
  label_x = 0
)
epi_EGFR_L858R_nonS_TS_ana_rm0_TN[[2]]
svg(file = "p_MicetumorNumber_vs_epi.svg", height = 24, width = 16)
(epi_KRAS_G12C_nonS_TS_ana_rm0_TN[[2]] | epi_KRAS_G12C_S_TS_ana_rm0_TN[[2]]) /
  (epi_KRAS_G12D_nonS_TS_ana_rm0_TN[[2]] | epi_KRAS_G12D_S_TS_ana_rm0_TN[[2]]) /
  (epi_BRAF_V600E_nonS_TS_ana_rm0_TN[[2]] |
    epi_BRAF_V600E_S_TS_ana_rm0_TN[[2]]) /
  (epi_EGFR_L858R_nonS_TS_ana_rm0_TN[[2]] |
    epi_EGFR_L858R_S_TS_ana_rm0_TN[[2]]) /
  plot_annotation(tag_levels = waiver())
dev.off()

jpeg(
  file = "p_MicetumorNumber_vs_epi.jpeg",
  height = 24,
  width = 16,
  unit = "in",
  res = 600
)
(epi_KRAS_G12C_nonS_TS_ana_rm0_TN[[2]] | epi_KRAS_G12C_S_TS_ana_rm0_TN[[2]]) /
  (epi_KRAS_G12D_nonS_TS_ana_rm0_TN[[2]] | epi_KRAS_G12D_S_TS_ana_rm0_TN[[2]]) /
  (epi_BRAF_V600E_nonS_TS_ana_rm0_TN[[2]] |
    epi_BRAF_V600E_S_TS_ana_rm0_TN[[2]]) /
  (epi_EGFR_L858R_nonS_TS_ana_rm0_TN[[2]] |
    epi_EGFR_L858R_S_TS_ana_rm0_TN[[2]]) /
  plot_annotation(tag_levels = waiver())
dev.off()

### relative tumor size; keep all point ####
## KRAS_G12D; keep near 0 CES values
G12D_tumorSize <- fread(paste0(location_data, "G12D_data_tumorSize_final.csv"))

epi_KRAS_G12D_nonS_TS_ana <- draw_scatter(
  epi_KRAS_G12D_nonS,
  G12D_tumorSize,
  MiceTumorData_type = "size",
  variant_onco = KRAS_G12D,
  MiceTumorlimit = c(0.5, 32),
  data_type = "nonS",
  max_overlaps = 5,
  label_y = 32,
  label_x = 10
)
epi_KRAS_G12D_nonS_TS_ana[[2]]
epi_KRAS_G12D_S_TS_ana <- draw_scatter(
  epi_KRAS_G12D_S,
  G12D_tumorSize,
  MiceTumorData_type = "size",
  variant_onco = KRAS_G12D,
  MiceTumorlimit = c(0.5, 32),
  data_type = "S",
  max_overlaps = 5,
  label_y = 32,
  label_x = 0
)
epi_KRAS_G12D_S_TS_ana[[2]]

## KRAS_G12C; rm near 0 CES values
G12C_tumorSize <- fread(paste0(location_data, "G12C_data_tumorSize_final.csv"))
epi_KRAS_G12C_S_TS_ana <- draw_scatter(
  epi_KRAS_G12C_S,
  G12C_tumorSize,
  MiceTumorData_type = "size",
  variant_onco = KRAS_G12C,
  MiceTumorlimit = c(0.5, 32),
  data_type = "S",
  max_overlaps = 5,
  label_y = 32,
  label_x = 0
)
epi_KRAS_G12C_S_TS_ana[[2]]
epi_KRAS_G12C_nonS_TS_ana <- draw_scatter(
  epi_KRAS_G12C_nonS,
  G12C_tumorSize,
  MiceTumorData_type = "size",
  variant_onco = KRAS_G12C,
  MiceTumorlimit = c(0.5, 32),
  data_type = "nonS",
  max_overlaps = 5,
  label_y = 32,
  label_x = 0
)
epi_KRAS_G12C_nonS_TS_ana[[2]]
## BRAF_V600E; rm near 0 CES values
BRAF_V600E_tumorSize <- fread(paste0(
  location_data,
  "BRAF_V600E_data_tumorSize_final.csv"
))
epi_BRAF_V600E_S_TS_ana <- draw_scatter(
  epi_BRAF_V600E_S,
  BRAF_V600E_tumorSize,
  MiceTumorData_type = "size",
  variant_onco = BRAF_V600E,
  MiceTumorlimit = c(0.25, 8),
  data_type = "S",
  max_overlaps = 5,
  label_y = 8,
  label_x = 0
)
epi_BRAF_V600E_S_TS_ana[[2]]
epi_BRAF_V600E_nonS_TS_ana <- draw_scatter(
  epi_BRAF_V600E_nonS,
  BRAF_V600E_tumorSize,
  MiceTumorData_type = "size",
  variant_onco = BRAF_V600E,
  MiceTumorlimit = c(0.25, 8),
  data_type = "nonS",
  max_overlaps = 5,
  label_y = 8,
  label_x = 0
)
epi_BRAF_V600E_nonS_TS_ana[[2]]
## EGFR L858R; rm near 0 CES values
EGFR_L858R_tumorSize <- fread(paste0(
  location_data,
  "EGFR_L858R_data_tumorSize_final.csv"
))
epi_EGFR_L858R_S_TS_ana <- draw_scatter(
  epi_EGFR_L858R_S,
  EGFR_L858R_tumorSize,
  MiceTumorData_type = "size",
  variant_onco = EGFR_L858R,
  MiceTumorlimit = c(0.25, 8),
  data_type = "S",
  max_overlaps = 5,
  label_y = 8,
  label_x = 0
)
epi_EGFR_L858R_S_TS_ana[[2]]
epi_EGFR_L858R_nonS_TS_ana <- draw_scatter(
  epi_EGFR_L858R_nonS,
  EGFR_L858R_tumorSize,
  MiceTumorData_type = "size",
  variant_onco = EGFR_L858R,
  MiceTumorlimit = c(0.25, 8),
  data_type = "nonS",
  max_overlaps = 5,
  label_y = 8,
  label_x = 1
)
epi_EGFR_L858R_nonS_TS_ana[[2]]
svg(file = "p_MiceTumorSize_vs_epi_keep0.svg", height = 24, width = 16)
(epi_KRAS_G12C_nonS_TS_ana[[2]] | epi_KRAS_G12C_S_TS_ana[[2]]) /
  (epi_KRAS_G12D_nonS_TS_ana[[2]] | epi_KRAS_G12D_S_TS_ana[[2]]) /
  (epi_BRAF_V600E_nonS_TS_ana[[2]] | epi_BRAF_V600E_S_TS_ana[[2]]) /
  (epi_EGFR_L858R_nonS_TS_ana[[2]] | epi_EGFR_L858R_S_TS_ana[[2]]) /
  plot_annotation(tag_levels = waiver())
dev.off()

# svg(file = "p_MiceTumorSize_vs_epi_v2.svg", height = 16, width = 11)
jpeg(
  file = "p_MicetumorSize_vs_epi_keep0.jpeg",
  height = 24,
  width = 16,
  unit = "in",
  res = 600
)
(epi_KRAS_G12C_nonS_TS_ana[[2]] | epi_KRAS_G12C_S_TS_ana[[2]]) /
  (epi_KRAS_G12D_nonS_TS_ana[[2]] | epi_KRAS_G12D_S_TS_ana[[2]]) /
  (epi_BRAF_V600E_nonS_TS_ana[[2]] | epi_BRAF_V600E_S_TS_ana[[2]]) /
  (epi_EGFR_L858R_nonS_TS_ana[[2]] | epi_EGFR_L858R_S_TS_ana[[2]]) /
  plot_annotation(tag_levels = waiver())
dev.off()

### relative tumor number; keep all point ####
## KRAS_G12D; keep near 0 CES values
G12D_tumorNumber <- fread(paste0(
  location_data,
  "G12D_data_tumorNumber_final.csv"
))

epi_KRAS_G12D_nonS_TS_ana_TN <- draw_scatter(
  epi_KRAS_G12D_nonS,
  G12D_tumorNumber,
  MiceTumorData_type = "number",
  variant_onco = KRAS_G12D,
  MiceTumorlimit = c(0.12, 8),
  data_type = "nonS",
  max_overlaps = 5,
  label_y = 8,
  label_x = 10
)
epi_KRAS_G12D_nonS_TS_ana_TN[[2]]
epi_KRAS_G12D_S_TS_ana_TN <- draw_scatter(
  epi_KRAS_G12D_S,
  G12D_tumorNumber,
  MiceTumorData_type = "number",
  variant_onco = KRAS_G12D,
  MiceTumorlimit = c(0.12, 8),
  data_type = "S",
  max_overlaps = 5,
  label_y = 8,
  label_x = 0
)
epi_KRAS_G12D_S_TS_ana_TN[[2]]

## KRAS_G12C; rm near 0 CES values
G12C_tumorNumber <- fread(paste0(
  location_data,
  "G12C_data_tumorNumber_final.csv"
))
epi_KRAS_G12C_S_TS_ana_TN <- draw_scatter(
  epi_KRAS_G12C_S,
  G12C_tumorNumber,
  MiceTumorData_type = "number",
  variant_onco = KRAS_G12C,
  MiceTumorlimit = c(0.12, 8),
  data_type = "S",
  max_overlaps = 5,
  label_y = 8,
  label_x = 0
)
epi_KRAS_G12C_S_TS_ana_TN[[2]]
epi_KRAS_G12C_nonS_TS_ana_TN <- draw_scatter(
  epi_KRAS_G12C_nonS,
  G12C_tumorNumber,
  MiceTumorData_type = "number",
  variant_onco = KRAS_G12C,
  MiceTumorlimit = c(0.12, 8),
  data_type = "nonS",
  max_overlaps = 5,
  label_y = 8,
  label_x = 0
)
epi_KRAS_G12C_nonS_TS_ana_TN[[2]]
## BRAF_V600E; rm near 0 CES values
BRAF_V600E_tumorNumber <- fread(paste0(
  location_data,
  "BRAF_V600E_data_tumorNumber_final.csv"
))
epi_BRAF_V600E_S_TS_ana_TN <- draw_scatter(
  epi_BRAF_V600E_S,
  BRAF_V600E_tumorNumber,
  MiceTumorData_type = "number",
  variant_onco = BRAF_V600E,
  MiceTumorlimit = c(0.12, 8),
  data_type = "S",
  max_overlaps = 5,
  label_y = 8,
  label_x = 0
)
epi_BRAF_V600E_S_TS_ana_TN[[2]]
epi_BRAF_V600E_nonS_TS_ana_TN <- draw_scatter(
  epi_BRAF_V600E_nonS,
  BRAF_V600E_tumorNumber,
  MiceTumorData_type = "number",
  variant_onco = BRAF_V600E,
  MiceTumorlimit = c(0.12, 8),
  data_type = "nonS",
  max_overlaps = 5,
  label_y = 8,
  label_x = 0
)
epi_BRAF_V600E_nonS_TS_ana_TN[[2]]
## EGFR L858R; rm near 0 CES values
EGFR_L858R_tumorNumber <- fread(paste0(
  location_data,
  "EGFR_L858R_data_tumorNumber_final.csv"
))
epi_EGFR_L858R_S_TS_ana_TN <- draw_scatter(
  epi_EGFR_L858R_S,
  EGFR_L858R_tumorNumber,
  MiceTumorData_type = "number",
  variant_onco = EGFR_L858R,
  MiceTumorlimit = c(0.12, 8),
  data_type = "S",
  max_overlaps = 5,
  label_y = 8,
  label_x = 0
)
epi_EGFR_L858R_S_TS_ana_TN[[2]]
epi_EGFR_L858R_nonS_TS_ana_TN <- draw_scatter(
  epi_EGFR_L858R_nonS,
  EGFR_L858R_tumorNumber,
  MiceTumorData_type = "number",
  variant_onco = EGFR_L858R,
  MiceTumorlimit = c(0.12, 8),
  data_type = "nonS",
  max_overlaps = 5,
  label_y = 8,
  label_x = 0
)
epi_EGFR_L858R_nonS_TS_ana_TN[[2]]
svg(file = "p_MicetumorNumber_vs_epi_keep0.svg", height = 24, width = 16)
(epi_KRAS_G12C_nonS_TS_ana_TN[[2]] | epi_KRAS_G12C_S_TS_ana_TN[[2]]) /
  (epi_KRAS_G12D_nonS_TS_ana_TN[[2]] | epi_KRAS_G12D_S_TS_ana_TN[[2]]) /
  (epi_BRAF_V600E_nonS_TS_ana_TN[[2]] |
    epi_BRAF_V600E_S_TS_ana_TN[[2]]) /
  (epi_EGFR_L858R_nonS_TS_ana_TN[[2]] |
    epi_EGFR_L858R_S_TS_ana_TN[[2]]) /
  plot_annotation(tag_levels = waiver())
dev.off()

jpeg(
  file = "p_MicetumorNumber_vs_epi_keep0.jpeg",
  height = 24,
  width = 16,
  unit = "in",
  res = 600
)
(epi_KRAS_G12C_nonS_TS_ana_TN[[2]] | epi_KRAS_G12C_S_TS_ana_TN[[2]]) /
  (epi_KRAS_G12D_nonS_TS_ana_TN[[2]] | epi_KRAS_G12D_S_TS_ana_TN[[2]]) /
  (epi_BRAF_V600E_nonS_TS_ana_TN[[2]] |
    epi_BRAF_V600E_S_TS_ana_TN[[2]]) /
  (epi_EGFR_L858R_nonS_TS_ana_TN[[2]] |
    epi_EGFR_L858R_S_TS_ana_TN[[2]]) /
  plot_annotation(tag_levels = waiver())
dev.off()


