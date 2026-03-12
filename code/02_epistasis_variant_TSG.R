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
# variant_A     variant_B   ces_A0    ces_B0 ces_A_on_B ces_B_on_A p_A_change p_B_change  p_epistasis expected_nAB_epistasis expected_nAB_null AB_epistatic_ratio ces_A_null ces_B_null   nA0   nB0   nAB   n00
# <char>        <char>    <num>     <num>      <num>      <num>      <num>      <num>        <num>                  <num>             <num>              <num>      <num>      <num> <int> <int> <int> <int>
#   1:      BRAF           NF1 26603.35  32.33761   942.5758      0.001 0.02766277 0.03177044 0.0007687001           1.287935e-01         6.8484275       1.880629e-02   25620.12   32.09869    99   497     0  6586
# 2:      BRAF         CMTR2 19259.62 686.60999     0.0010      0.001 0.78046350 0.81724238 0.6378346428           3.164870e-07         0.4079435       7.758108e-07   18943.23  683.12476    15    40     0  1414
# 3:      BRAF          MSH2 25391.23  50.55111     0.0010      0.001 0.50264115 0.50864243 0.2556700461           1.352591e-05         1.3399209       1.009456e-05   25214.86   50.20612    99    99     0  7099
# 4:      BRAF CDKN2A.p14arf 25909.93  10.29680     0.0010      0.001 0.81158911 0.81141951 0.6585327624           2.065599e-05         0.4200404       4.917621e-05   25853.56   10.22410   120    30     0  8482
# 5:      EGFR CDKN2A.p14arf 94217.87  10.64417   992.5697      0.001 0.28557740 0.30321398 0.0854395625           1.325079e-02         2.3573431       5.621071e-03   94037.78   10.22410   676    30     0  7926
# n_total ci_low_95_ces_A0 ci_high_95_ces_A0 ci_low_95_ces_B0 ci_high_95_ces_B0 ci_low_95_ces_A_on_B ci_high_95_ces_A_on_B ci_low_95_ces_B_on_A ci_high_95_ces_B_on_A                      run_type data_type
# <int>            <num>             <num>            <num>             <num>                <num>                 <num>                <num>                 <num>                        <char>    <char>
#   1:    7182         21702.43          32195.92        29.574408          35.26279                   NA              14406.74                   NA              17.21931             BRAF_V600E_vs_NF1       Pan
# 2:    1469         11082.92          30715.43       495.212276         921.93211                   NA             149211.93                   NA            6658.46458           BRAF_V600E_vs_CMTR2       Pan
# 3:    7297         20712.42          30726.71        41.236169          61.17357                   NA              70804.68                   NA             143.43837            BRAF_V600E_vs_MSH2       Pan
# 4:    8632         21537.04          30812.98         7.036941          14.43134                   NA             240374.58                   NA              94.51096 BRAF_V600E_vs_CDKN2A.p16INK4a       Pan
# 5:    8632         87290.72         101502.33         7.275141          14.91995                   NA             144266.25                   NA              16.60463 EGFR_L858R_vs_CDKN2A.p16INK4a       Pan
epi_final_pan[nAB == "1"]
# variant_A       variant_B    ces_A0    ces_B0 ces_A_on_B  ces_B_on_A p_A_change p_B_change p_epistasis expected_nAB_epistasis expected_nAB_null AB_epistatic_ratio ces_A_null ces_B_null   nA0   nB0   nAB   n00
# <char>          <char>     <num>     <num>      <num>       <num>      <num>      <num>       <num>                  <num>             <num>              <num>      <num>      <num> <int> <int> <int> <int>
#   1:      KRAS           FBXW7  88101.18  36.53048  1425.7583   13.147406  0.5132914  0.5153277 0.073741042              0.9973779         5.1054818          0.1953543   87474.99   35.80307   367    83     1  5601
# 2:      KRAS           CMTR2  41281.53 676.33268  1662.7592 1134.966128  1.0000000  1.0000000 1.000000000              0.9785996         1.1444596          0.8550757   40864.47  683.12476    43    39     1  1386
# 3:      KRAS   CDKN2A.p14arf  50900.64  10.63118   967.0665    4.736560  0.7316014  0.6946555 0.175569863              1.0029539         4.0430088          0.2480712   50808.86   10.22410  1154    29     1  7448
# 4:      BRAF           FBXW7  26099.24  35.63838  1107.0788   59.688208  1.0000000  1.0000000 0.996242250              0.9999047         1.1653620          0.8580207   25921.65   35.80307    83    83     1  5885
# 5:      BRAF           RNF43  23530.72  37.58202  1091.9836   82.831535  1.0000000  1.0000000 0.997428356              0.9693933         0.8668666          1.1182727   23407.52   37.87366    81    69     1  6354
# 6:      BRAF           RBM10  23935.27  99.75628  1599.7886   30.649553  0.4922832  0.5101819 0.069841803              0.9977945         5.0655143          0.1969779   23025.34   99.30098    64   409     1  4763
# 7:      BRAF             NF2  25318.74  99.57811  1036.7676  226.877708  1.0000000  1.0000000 0.999267883              0.9689917         0.8417035          1.1512269   25214.86  100.45819    98    62     1  7136
# 8:      BRAF           PALB2  25439.69  52.44701  1094.1697   58.618293  0.9486069  0.9955446 0.792701939              0.9922611         1.6977387          0.5844604   25214.86   52.48703    98   125     1  7073
# 9:      BRAF   CDKN2A.p14arf  25548.47  24.71247  1190.5917   17.497663  0.9720338  0.9594761 0.523277189              0.9929644         2.5909126          0.3832489   25225.27   24.66171    98   190     1  7005
# 10:      BRAF CDKN2A.p16INK4a  25487.29  30.54260  1237.0489   17.098397  0.8531336  0.8658082 0.331874213              1.0018653         3.2251653          0.3106400   25083.75   30.44990    98   238     1  6998
# 11:      EGFR            TSC1 100732.01  46.43883   963.2401    9.991092  0.1771249  0.1501405 0.004462914              0.9847883         8.2797966          0.1189387  100031.09   44.86841   608    99     1  6613
# 12:      EGFR           CMTR2 111603.89 698.22283  2872.6316  362.949784  0.9197585  0.8022751 0.298766747              1.0445320         3.6405988          0.2869121  110259.70  683.12476   137    39     1  1292
# 13:      EGFR             NF2 100825.81 103.39192  1044.7693   35.749803  0.4955304  0.4728220 0.061964334              0.9769771         5.2286293          0.1868515  100372.22  100.45819   608    62     1  6626
# n_total ci_low_95_ces_A0 ci_high_95_ces_A0 ci_low_95_ces_B0 ci_high_95_ces_B0 ci_low_95_ces_A_on_B ci_high_95_ces_A_on_B ci_low_95_ces_B_on_A ci_high_95_ces_B_on_A                      run_type data_type
# <int>            <num>             <num>            <num>             <num>                <num>                 <num>                <num>                 <num>                        <char>    <char>
#   1:    6052         79402.01          97415.81        29.224206          44.96274                   NA             114678.18           0.19262011              59.98923            KRAS_G12D_vs_FBXW7       Pan
# 2:    1469         30247.33          54716.18       485.664846         911.39152                   NA             322990.90          38.85966943            5161.56752            KRAS_G12D_vs_CMTR2       Pan
# 3:    8632         48018.48          53896.11         7.216165          14.98297                   NA              85720.16           0.07991724              21.54929  KRAS_G12C_vs_CDKN2A.p16INK4a       Pan
# 4:    6052         20901.49          32078.88        28.508312          43.86088                   NA             155833.59           1.97444457             269.81534           BRAF_V600E_vs_FBXW7       Pan
# 5:    6505         18795.70          28996.99        29.398152          47.16110                   NA             182457.13           3.07584684             373.56039           BRAF_V600E_vs_RNF43       Pan
# 6:    5237         18565.50          30252.68        90.394467         109.74032                   NA              31101.05                   NA             159.14171           BRAF_V600E_vs_RBM10       Pan
# 7:    7297         20652.43          30639.84        76.803929         126.46269                   NA             205603.68           9.08036918            1018.84168             BRAF_V600E_vs_NF2       Pan
# 8:    7297         20751.25          30788.89        43.780195          62.18389                   NA              93261.62           1.21047797             267.14579           BRAF_V600E_vs_PALB2       Pan
# 9:    7294         20834.33          30915.15        21.361351          28.39277                   NA              67852.39                   NA              82.00298   BRAF_V600E_vs_CDKN2A.p14arf       Pan
# 10:    7335         20787.24          30847.42        26.824589          34.58906                   NA              53042.08                   NA              80.79118 BRAF_V600E_vs_CDKN2A.p16INK4a       Pan
# 11:    7321         92938.72         108949.62        37.881183          56.19728                   NA              82984.66           0.15138992              45.54734            EGFR_L858R_vs_TSC1       Pan
# 12:    1469         93986.25         131292.33       501.363990         940.88356                   NA             233596.43           3.75250371            1666.20235           EGFR_L858R_vs_CMTR2       Pan
# 13:    7297         93025.50         109050.70        79.744728         131.30655                   NA             127446.59           1.02945323             161.20324             EGFR_L858R_vs_NF2       Pan
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


# ## plot significant; use bar instead of point; combine two together ####
# epi_signDiff <- rbind(
#   #epi_EGFR_L858R[variant_B %in% c("TP53")],
#   epi_KRAS_G12C[variant_B %in% c("NF1")]
# )
# variant_label_signDiff <- epi_signDiff[,
#   .(max_y = max(ci_high_95_ces_B_on_A)),
#   by = variant_B
# ]
# # p_signDiff_signifi_bar <- ggplot(
# #   data = epi_signDiff,
# #   aes(x = variant_B, y = ces_B_on_A, fill = data_type)
# # ) +
# #   ylab("Cancer effect size") +
# #   geom_bar(
# #     stat = "identity",
# #     color = "black",
# #     position = position_dodge(width = 0.9)
# #   ) + ## color = errorBarColor,shape = 21, size =6
# #   geom_errorbar(
# #     aes(ymin = ci_low_95_ces_B_on_A_forplot, ymax = ci_high_95_ces_B_on_A),
# #     width = 0.3,
# #     na.rm = T,
# #     linewidth = 0.6,
# #     position = position_dodge(width = 0.9)
# #   ) + ## , color = errorBarColor
# #   # Add variant_name labels above the bars
# #   geom_text(
# #     data = variant_label_signDiff,
# #     aes(x = variant_B, y = max(variant_label$max_y) * 3, label = variant_B),
# #     inherit.aes = F,
# #     vjust = 0,
# #     size = 5,
# #     fontface = "bold"
# #   ) +
# #   scale_x_discrete(labels = rep(c("Never-smoker", "Smoker"), length.out = 2)) +
# #   scale_fill_manual(values = c("Smo" = "#D2B48C", "nonSmo" = "pink")) + ## ,"Gray" = "gray" ## "All" = "black", ## brown looks like red,so use light brown hex code #D2B48C instead
# #   scale_y_continuous(
# #     trans = scales::pseudo_log_trans(base = 10),
# #     breaks = c(0.01, 1, 100, 10000, 1000000),
# #     labels = scales::label_comma(accuracy = 1)
# #   ) +
# #   theme_classic() + ##theme_classic()
# #   theme(
# #     legend.direction = "horizontal",
# #     legend.position = "bottom",
# #     legend.box = "vertical",
# #     legend.text = element_text(size = 18),
# #     legend.title = element_text(size = 18),
# #     axis.text.x = element_blank(),
# #     axis.text.y = element_text(size = 18, face = "bold"),
# #     axis.title.y = element_text(size = 18, face = "bold"),
# #     axis.title.x = element_blank()
# #   ) +
# #   guides(fill = F)
# #
# # svg(file = "CES_select_VariantGeneLevel.svg", height = 10.5, width = 5)
# # p_signDiff_signifi_bar
# # dev.off()
#
# ## use bar instead of point; combine two together' do a square-root transformation ####
# epi_signDiff[,
#   ci_low_95_ces_B_on_A_forplot_v2 := ifelse(
#     is.na(ci_low_95_ces_B_on_A),
#     0.001,
#     ci_low_95_ces_B_on_A
#   )
# ]
# min(epi_signDiff[, ci_low_95_ces_B_on_A_forplot_v2]) ## 19.2795
# max(epi_signDiff[, ci_high_95_ces_B_on_A]) ## 3212.489
# p_signDiff_signifi_bar <- ggplot(
#   data = epi_signDiff,
#   aes(x = variant_B, y = ces_B_on_A, fill = data_type)
# ) +
#   ylab("Cancer effect size") +
#   geom_bar(
#     stat = "identity",
#     color = "black",
#     position = position_dodge(width = 0.9)
#   ) + ## color = errorBarColor,shape = 21, size =6
#   geom_errorbar(
#     aes(ymin = ci_low_95_ces_B_on_A_forplot_v2, ymax = ci_high_95_ces_B_on_A),
#     width = 0.3,
#     na.rm = T,
#     linewidth = 0.6,
#     position = position_dodge(width = 0.9)
#   ) + ## , color = errorBarColor
#   # Add variant_name labels above the bars
#   geom_text(
#     data = variant_label_signDiff,
#     aes(x = variant_B, y = max(variant_label$max_y) * 1.4, label = variant_B),
#     inherit.aes = F,
#     vjust = 0,
#     size = 5,
#     fontface = "bold"
#   ) +
#   scale_x_discrete(labels = rep(c("Never-smoker", "Smoker"), length.out = 2)) +
#   scale_fill_manual(values = c("Smo" = "#D2B48C", "nonSmo" = "pink")) + ## ,"Gray" = "gray" ## "All" = "black", ## brown looks like red,so use light brown hex code #D2B48C instead
#   #scale_color_manual(values = c("Gray" = "gray","Black" = "black"))+
#   # scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = F))+
#   # scale_x_discrete(labels = c("KRAS G12C Never-smoker","KRAS G12C Smoker"))+ ## change the x-axis label
#   scale_y_continuous(
#     trans = "sqrt",
#     breaks = c(
#       5000,
#       10000,
#       20000,
#       50000,
#       100000,
#       200000,
#       400000,
#       800000,
#       1200000
#     ),
#     labels = scales::label_comma(accuracy = 1),
#     expand = expansion(mult = c(0.05, 0.05))
#   ) +
#   theme_classic() + ##theme_classic()
#   theme(
#     legend.direction = "horizontal",
#     legend.position = "bottom",
#     legend.box = "vertical",
#     legend.text = element_text(size = 18),
#     legend.title = element_text(size = 18),
#     axis.text.x = element_blank(), #element_text( size=18,face = "bold",angle = 30, hjust = 1), ## , angle = 30
#     axis.text.y = element_text(size = 18, face = "bold"),
#     axis.title.y = element_text(size = 18, face = "bold"),
#     axis.title.x = element_blank()
#   ) +
#   guides(fill = F)
#
# svg(
#   file = "CES_select_VariantGeneLevel_siginifi_byNonOverlapping_signDiff_onlyNS_S_bar_v4.svg",
#   height = 10.5,
#   width = 5
# ) ## 6inch = 1800 pixels. 1278 pixels, 1000
# p_signDiff_signifi_bar
# dev.off()

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
    # labs(
    #   # title = paste(
    #   #   expression(rho),
    #   #   ":",
    #   #   round(spearman_r, 2),
    #   #   ";",
    #   #   "P =",
    #   #   format(spearman_p, scientific = TRUE, digits = 2),
    #   #   "\n",
    #   #   "Pearson's",
    #   #   ":",
    #   #   round(pearson_r, 2),
    #   #   ";",
    #   #   "P =",
    #   #   format(pearson_p, scientific = F)
    #   # )
    #   title = bquote(
    #     atop(
    #     "Spearman's " *rho * ": " * .(round(spearman_r, 2)) * "; " * italic(P) * " = " *.(format(spearman_p, digits = 3)),
    #       "Pearson's " * italic(r) * ":" * .(round(pearson_r, 2)) * "; " * italic(P) * " = " * .(format(pearson_p, digits = 3)))
    #     )
    # )+
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

# # check the underlying mutation rate of CMTR2, MSH2, NF1
# head(cesa$gene_rates)
# cesa$gene_rates[gene %in% c("CMTR2", "MSH2", "NF1")]
