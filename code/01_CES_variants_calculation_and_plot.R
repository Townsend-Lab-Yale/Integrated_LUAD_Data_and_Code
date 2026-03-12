# rm(list=ls())
if (!require("cancereffectsizeR")) {
  ## use cancereffectsizeR v2.10.2 here
  remotes::install_github(
    "Townsend-Lab-Yale/cancereffectsizeR@v2.10.2",
    dependencies = T,
    force = T
  )
  packageVersion("cancereffectsizeR") ## v2.10.2
}
if (!require("ces.refset.hg19")) {
  ## use ces.refset.hg19 as the reference dataset here
  remotes::install_github(
    "Townsend-Lab-Yale/ces.refset.hg19@*release",
    dependencies = T,
    force = T
  )
  require("ces.refset.hg19")
  packageVersion("ces.refset.hg19") ## v1.1.3
}

##  import other libraries
library(data.table)
library(dplyr)
library(rtracklayer)
library(stringr)
library(ggplot2)


location_output <- "../integrated_data/"
location_data <- "../data/"
rdata_output <- "./R_Data/"
## read maf file
maf_file <- read.csv(paste0(location_output, "merged_luad_maf.txt"))
colnames(maf_file)[2] <- 'Tumor_Sample_Barcode'
colnames(maf_file)[7] <- 'Tumor_Allele'
maf_list <- split(maf_file, maf_file$Source)

liftover_file = paste0(location_data, "hg38ToHg19.over.chain")

## preload all maf file ####
### WGS data
NCI_maf <- preload_maf(
  maf = maf_list$NCI,
  refset = ces.refset.hg19,
  keep_extra_columns = T
) ## chain_file = liftover_file,
NCI_maf <- NCI_maf[is.na(problem)]
NCI_maf <- NCI_maf[
  germline_variant_site == F &
    (repetitive_region == F | cosmic_site_tier %in% 1:3)
]

### WES+WGS
Broad_maf <- preload_maf(
  maf = maf_list$Broad,
  refset = ces.refset.hg19,
  keep_extra_columns = T
) ## chain_file = liftover_file,
Broad_maf <- Broad_maf[is.na(problem)]
Broad_maf <- Broad_maf[
  germline_variant_site == F &
    (repetitive_region == F | cosmic_site_tier %in% 1:3)
]
Broad_maf <- Broad_maf[
  !Unique_Patient_Identifier %in% c("LUAD-B01169", "LUAD-D01382")
]

### WES
MSK2015_maf <- preload_maf(
  maf = maf_list$MSK2015,
  refset = ces.refset.hg19,
  keep_extra_columns = T
) ## chain_file = liftover_file,
MSK2015_maf <- MSK2015_maf[is.na(problem)]
MSK2015_maf <- MSK2015_maf[
  germline_variant_site == F &
    (repetitive_region == F | cosmic_site_tier %in% 1:3)
]

OncoSG_maf <- preload_maf(
  maf = maf_list$OncoSG,
  refset = ces.refset.hg19,
  keep_extra_columns = T
) ## chain_file = liftover_file,
OncoSG_maf <- OncoSG_maf[is.na(problem)]
OncoSG_maf <- OncoSG_maf[
  germline_variant_site == F &
    (repetitive_region == F | cosmic_site_tier %in% 1:3)
]

TCGA_maf <- preload_maf(
  maf = maf_list$TCGA,
  refset = ces.refset.hg19,
  chain_file = liftover_file,
  keep_extra_columns = T
)
TCGA_maf <- TCGA_maf[is.na(problem)]
TCGA_maf <- TCGA_maf[
  germline_variant_site == F &
    (repetitive_region == F | cosmic_site_tier %in% 1:3)
]

TracerX_maf <- preload_maf(
  maf = maf_list$TracerX,
  refset = ces.refset.hg19,
  keep_extra_columns = T
) ## chain_file = liftover_file,
TracerX_maf <- TracerX_maf[is.na(problem)]
TracerX_maf <- TracerX_maf[
  germline_variant_site == F &
    (repetitive_region == F | cosmic_site_tier %in% 1:3)
]

CPTAC_maf <- cancereffectsizeR::preload_maf(
  maf = maf_list$CPTAC,
  refset = ces.refset.hg19,
  keep_extra_columns = T
)
CPTAC_maf <- CPTAC_maf[is.na(problem)]
CPTAC_maf <- CPTAC_maf[
  germline_variant_site == F &
    (repetitive_region == F | cosmic_site_tier %in% 1:3)
]

Yale_maf <- cancereffectsizeR::preload_maf(
  maf = maf_list$Yale,
  refset = ces.refset.hg19,
  keep_extra_columns = T
)
Yale_maf <- Yale_maf[is.na(problem)]
Yale_maf <- Yale_maf[
  germline_variant_site == F &
    (repetitive_region == F | cosmic_site_tier %in% 1:3)
]

### WES; TGS
MSK2018_maf <- preload_maf(
  maf = maf_list$MSK2018,
  refset = ces.refset.hg19,
  keep_extra_columns = T
)
MSK2018_maf <- MSK2018_maf[is.na(problem)]

### TGS data
FMAD_maf <- preload_maf(
  maf = maf_list$`FM-AD`,
  ces.refset.hg19,
  chain_file = liftover_file,
  keep_extra_columns = T
)
FMAD_maf <- FMAD_maf[is.na(problem)]

Genie_maf <- preload_maf(
  maf = maf_list$Genie,
  refset = ces.refset.hg19,
  keep_extra_columns = T
)
Genie_maf <- Genie_maf[is.na(problem)]

MSK2017_maf <- preload_maf(
  maf = maf_list$MSK2017,
  refset = ces.refset.hg19,
  keep_extra_columns = T
)
MSK2017_maf <- MSK2017_maf[is.na(problem)]

TSP_maf <- preload_maf(
  maf = maf_list$TSP,
  refset = ces.refset.hg19,
  keep_extra_columns = T
)
TSP_maf <- TSP_maf[is.na(problem)]
# why remove
TSP_maf <- TSP_maf[
  !Unique_Patient_Identifier %in%
    c("luad_tsp_16929", "luad_tsp_16901", "luad_tsp_16875", "luad_tsp_16915")
]


## extract necessary coverage info
### read in information necessary for loading maf files (exome/genome and coverage intervals for TGS)
broad_exome_or_genome <- fread(paste0(
  location_data,
  "luad_broad/data_clinical_sample.txt"
))[-(1:4), c('Sample Identifier', 'Platform')]
broad_exome_or_genome$Platform <- sapply(
  broad_exome_or_genome$Platform,
  function(x) {
    if (str_detect(x, 'WGS')) {
      return('WGS')
    } else return('WES')
  }
)
msk2017_panels_used <- fread(paste0(
  location_data,
  "lung_msk_2017/data_clinical_sample.txt"
))[-(1:4), c('Sample Identifier', 'Gene Panel')]
msk2018_panels_used <- fread(paste0(
  location_data,
  "nsclc_pd1_msk_2018/data_clinical_sample.txt"
))[-(1:4), c('Sample Identifier', 'Gene Panel')]
genie_panels_used <- fread(paste0(
  location_data,
  "genie_9/data_clinical_sample.txt"
))[-(1:4), c('Sample Identifier', 'Sequence Assay ID')]

### read in genes included in each panel and creating GRANGES object (bed files)  to pass into covered_regions parameter of load_maf
### once the granges are exported once, these functions don't need to be run anymore
# gene_granges <- rtracklayer::import(paste0(location_data, "gencode.v38lift37.basic.annotation.gtf")) ## too large
load(paste0(rdata_output, "gencode.v38lift37.basic.annotation.gtf.Rdata"))

location_bed <- paste0(location_data, 'bed_files/')
if (!dir.exists(location_bed)) {
  dir.create(location_bed)
}

location_gene_panels <- paste0(location_data, 'gene_panels/')

if (!file.exists(paste0(location_bed, "fmad_targets.bed"))) {
  fmad_genes <- unique(
    fread(paste0(location_gene_panels, "foundation_one.txt"))$Hugo_Symbol
  )
  fmad_granges <- gene_granges[gene_granges$gene_name %in% fmad_genes, ]
  fmad_granges <- fmad_granges[fmad_granges$type %in% c('CDS', 'stop_codon'), ]
  fmad_gr_clean <- cancereffectsizeR:::clean_granges_for_cesa(
    cesa = CESAnalysis(ces.refset.hg19),
    gr = fmad_granges
  )
  export(fmad_gr_clean, paste0(location_bed, "fmad_targets.bed"))
}

if (!file.exists(paste0(location_bed, "tsp_targets.bed"))) {
  tsp_genes <- unique(
    fread(paste0(location_gene_panels, "tsp.txt"))$Hugo_Symbol
  )
  tsp_granges <- gene_granges[gene_granges$gene_name %in% tsp_genes, ]
  tsp_granges <- tsp_granges[tsp_granges$type %in% c('CDS', 'stop_codon'), ]
  tsp_gr_clean <- cancereffectsizeR:::clean_granges_for_cesa(
    cesa = CESAnalysis(ces.refset.hg19),
    gr = tsp_granges
  )
  export(tsp_gr_clean, paste0(location_bed, "tsp_targets.bed"))
}

if (!file.exists(paste0(location_bed, "msk341_targets.bed"))) {
  msk_341_genes <- unique(
    fread(paste0(location_gene_panels, "msk341.txt"))$Hugo_Symbol
  )
  msk341_granges <- gene_granges[gene_granges$gene_name %in% msk_341_genes, ]
  msk341_granges <- msk341_granges[
    msk341_granges$type %in% c('CDS', 'stop_codon'),
  ]
  msk341_gr_clean <- cancereffectsizeR:::clean_granges_for_cesa(
    cesa = CESAnalysis(ces.refset.hg19),
    gr = msk341_granges
  )
  export(msk341_gr_clean, paste0(location_bed, "msk341_targets.bed"))
}

if (!file.exists(paste0(location_bed, "msk410_targets.bed"))) {
  msk_410_genes <- unique(
    fread(paste0(location_gene_panels, "msk410.txt"))$Hugo_Symbol
  )
  msk410_granges <- gene_granges[gene_granges$gene_name %in% msk_410_genes, ]
  msk410_granges <- msk410_granges[
    msk410_granges$type %in% c('CDS', 'stop_codon'),
  ]
  msk410_gr_clean <- cancereffectsizeR:::clean_granges_for_cesa(
    cesa = CESAnalysis(ces.refset.hg19),
    gr = msk410_granges
  )
  export(msk410_gr_clean, paste0(location_bed, "msk410_targets.bed"))
}

if (!file.exists(paste0(location_bed, "msk468_targets.bed"))) {
  msk_468_genes <- unique(
    fread(paste0(location_gene_panels, "msk468.txt"))$Hugo_Symbol
  )
  msk468_granges <- gene_granges[gene_granges$gene_name %in% msk_468_genes, ]
  msk468_granges <- msk468_granges[
    msk468_granges$type %in% c('CDS', 'stop_codon'),
  ]
  msk468_gr_clean <- cancereffectsizeR:::clean_granges_for_cesa(
    cesa = CESAnalysis(ces.refset.hg19),
    gr = msk468_granges
  )
  export(msk468_gr_clean, paste0(location_bed, "msk468_targets.bed"))
}

### this should be run everytime so that we don't have to save a bunch of grange files
genie_panel_coverage_df <- fread(paste0(
  location_data,
  "genie_9/genomic_information.txt"
))[, c(
  'Chromosome',
  'Start_Position',
  'End_Position',
  'Hugo_Symbol',
  'Feature_Type',
  'SEQ_ASSAY_ID'
)]
genie_granges_list <- makeGRangesListFromDataFrame(
  genie_panel_coverage_df,
  ignore.strand = T,
  seqnames.field = 'Chromosome',
  start.field = 'Start_Position',
  end.field = 'End_Position',
  split.field = 'SEQ_ASSAY_ID'
)
seqlevels(genie_granges_list, pruning.mode = "fine") <- c(1:22, 'X', 'Y')

genie_panel_genes <- genie_panel_coverage_df[, c('Hugo_Symbol', 'SEQ_ASSAY_ID')]
genie_panel_genes <- genie_panel_genes[!duplicated(genie_panel_genes)]
genie_panel_genes_list <- split(
  genie_panel_genes$Hugo_Symbol,
  genie_panel_genes$SEQ_ASSAY_ID
)

### split maf files when different sequencing methods are used
Broad_maf <- merge(
  Broad_maf,
  broad_exome_or_genome,
  by.x = 'Unique_Patient_Identifier',
  by.y = 'Sample Identifier'
)
Broad_maf <- split(Broad_maf, Broad_maf$Platform)

MSK2017_maf <- merge(
  MSK2017_maf,
  msk2017_panels_used,
  by.x = 'Unique_Patient_Identifier',
  by.y = 'Sample Identifier'
)
MSK2017_maf <- split(MSK2017_maf, MSK2017_maf$`Gene Panel`)

MSK2018_maf <- merge(
  MSK2018_maf,
  msk2018_panels_used,
  by.x = 'Unique_Patient_Identifier',
  by.y = 'Sample Identifier'
)
MSK2018_maf <- split(MSK2018_maf, MSK2018_maf$`Gene Panel`)

Genie_maf <- merge(
  Genie_maf,
  genie_panels_used,
  by.x = 'Unique_Patient_Identifier',
  by.y = 'Sample Identifier'
)
Genie_maf <- split(Genie_maf, Genie_maf$`Sequence Assay ID`) ## Sequence Assay ID.x


## load maf files into cesa ojbect ####
cesa <- CESAnalysis(ces.refset.hg19)
### consider using covered regions padding when variants are outside the intervals
#### WES data
cesa <- load_maf(cesa, maf = Broad_maf$WES, maf_name = 'Broad_WES')
cesa <- load_maf(cesa, maf = MSK2015_maf, maf_name = 'MSK2015')
cesa <- load_maf(cesa, maf = OncoSG_maf, maf_name = 'OncoSG')
cesa <- load_maf(cesa, maf = TCGA_maf, maf_name = 'TCGA')
cesa <- load_maf(cesa, maf = TracerX_maf, maf_name = 'TracerX')
cesa <- load_maf(cesa, maf = CPTAC_maf, maf_name = 'CPTAC') ## 108 samples
cesa <- load_maf(cesa, maf = Yale_maf, maf_name = 'Yale') ## 108 samples

#### TGS data
cesa <- load_maf(
  cesa,
  maf = FMAD_maf,
  coverage = 'targeted',
  covered_regions = paste0(location_bed, "fmad_targets.bed"),
  covered_regions_name = 'fmad_regions',
  covered_regions_padding = 100,
  maf_name = 'FMAD'
) #padding based on 23 variants having distance from interval between 10 and 100.

for (i in 1:length(Genie_maf)) {
  cesa <- load_maf(
    cesa,
    maf = Genie_maf[i][[1]],
    coverage = 'targeted',
    covered_regions = genie_granges_list[names(Genie_maf)[i]][[1]],
    covered_regions_name = paste0(names(Genie_maf)[i], '_regions'),
    covered_regions_padding = 100,
    maf_name = names(Genie_maf)[i]
  )
}

cesa <- load_maf(
  cesa,
  maf = MSK2017_maf$IMPACT341,
  coverage = 'targeted',
  covered_regions = paste0(location_bed, "msk341_targets.bed"),
  covered_regions_name = 'msk341_regions',
  covered_regions_padding = 100,
  maf_name = 'MSK2017_IMPACT341'
)
cesa <- load_maf(
  cesa,
  maf = MSK2017_maf$IMPACT410,
  coverage = 'targeted',
  covered_regions = paste0(location_bed, "msk410_targets.bed"),
  covered_regions_name = 'msk410_regions',
  covered_regions_padding = 100,
  maf_name = 'MSK2017_IMPACT410'
)
cesa <- load_maf(
  cesa,
  maf = MSK2018_maf$IMPACT341,
  coverage = 'targeted',
  covered_regions = paste0(location_bed, "msk341_targets.bed"),
  covered_regions_name = 'msk341_regions',
  covered_regions_padding = 100,
  maf_name = 'MSK2018_IMPACT341'
)
cesa <- load_maf(
  cesa,
  maf = MSK2018_maf$IMPACT410,
  coverage = 'targeted',
  covered_regions = paste0(location_bed, "msk410_targets.bed"),
  covered_regions_name = 'msk410_regions',
  covered_regions_padding = 100,
  maf_name = 'MSK2018_IMPACT410'
)
cesa <- load_maf(
  cesa,
  maf = MSK2018_maf$IMPACT468,
  coverage = 'targeted',
  covered_regions = paste0(location_bed, "msk468_targets.bed"),
  covered_regions_name = 'msk468_regions',
  covered_regions_padding = 100,
  maf_name = 'MSK2018_IMPACT468'
)
cesa <- load_maf(
  cesa,
  maf = TSP_maf,
  coverage = 'targeted',
  covered_regions = paste0(location_bed, "tsp_targets.bed"),
  covered_regions_name = 'tsp_regions',
  covered_regions_padding = 100,
  maf_name = 'TSP'
)

#### WGS data
cesa <- load_maf(
  cesa,
  maf = Broad_maf$WGS,
  coverage = 'genome',
  maf_name = 'Broad_WGS'
)
cesa <- load_maf(cesa, maf = NCI_maf, coverage = 'genome', maf_name = 'NCI')

save_cesa(cesa, paste0(rdata_output, "load_maf_cesa_WES_TGS_WGS.rds"))

cesa_smoking_w_panel <- cesa
cesa_nonsmoking_w_panel <- cesa


nrow(cesa$samples) ## 9459
write.csv(
  cesa$samples,
  file = paste0(location_output, "sample_sequencingType_source_info.csv")
)
maf_mutation_gene <- cesa$maf
sapply(cesa$maf, is.list)
# Unique_Patient_Identifier                Chromosome            Start_Position          Reference_Allele              Tumor_Allele              variant_type
# FALSE                     FALSE                     FALSE                     FALSE                     FALSE                     FALSE
# variant_id                     genes                  top_gene           top_consequence               prelift_chr             prelift_start
# FALSE                      TRUE                     FALSE                     FALSE                     FALSE                     FALSE
# liftover_strand_flip
# FALSE
maf_mutation_gene[] <- lapply(maf_mutation_gene, function(col) {
  if (is.list(col)) {
    sapply(col, toString)
  } else {
    col
  }
})
write.csv(
  maf_mutation_gene,
  file = paste0(location_output, "merged_luad_maf_TransvertToGene.csv")
)
# head(maf_mutation_gene)
# head(maf_file)
## calculate mutation rates for all samples, smokers, and never-smokers  ####
treated_samples_signature_exclusions <- suggest_cosmic_signature_exclusions(
  cancer_type = "LUAD",
  treatment_naive = F
)
untreated_samples_signature_exclusions <- suggest_cosmic_signature_exclusions(
  cancer_type = "LUAD",
  treatment_naive = T
)

clin_df = fread(paste0(location_output, 'merged_luad_clinical.txt'))


## use all sample; include TGS data ###
all_samples = cesa$samples[, Unique_Patient_Identifier] ## coverage %in% c('exome','genome')
all_treated_samples_for_sigs = all_samples[
  all_samples %in% clin_df[(Treatment)]$`Sample ID`
]
all_untreated_samples_for_sigs = all_samples[
  all_samples %in% clin_df[(!Treatment)]$`Sample ID`
]
all_treatment_not_indicated_samples_for_sigs = all_samples[
  all_samples %in% clin_df[is.na(Treatment)]$`Sample ID`
]
sample_with_clinInfo <- c(
  all_treated_samples_for_sigs,
  all_untreated_samples_for_sigs,
  all_treatment_not_indicated_samples_for_sigs
)
### calculate trinucleotide mutation rates for pan data, needed for determinting smoking history from SBS4 weights####
cesa <- trinuc_mutation_rates(
  cesa,
  signature_set = ces.refset.hg19$signatures$COSMIC_v3.2,
  signature_exclusions = treated_samples_signature_exclusions,
  samples = c(
    all_treated_samples_for_sigs,
    all_treatment_not_indicated_samples_for_sigs
  ),
  cores = 4
)
cesa <- trinuc_mutation_rates(
  cesa,
  signature_set = ces.refset.hg19$signatures$COSMIC_v3.2,
  signature_exclusions = untreated_samples_signature_exclusions,
  samples = all_untreated_samples_for_sigs,
  cores = 4
)

### use mutational signature convolution on WES/WGS and clinical data for TGS to split the samples into smokers and nonsmokers ####
# Determining smoking history from SBS4 weights...
bio_weights <- cesa$mutational_signatures$biological_weights
bio_weights_unblended <- bio_weights[bio_weights$group_avg_blended == F]

snv_counts <- cesa$maf[
  variant_type == 'snv',
  .N,
  by = "Unique_Patient_Identifier"
]
good_samples <- snv_counts[N > 50, Unique_Patient_Identifier] ## only 1165 samples among 9027 samples have >50 SNV


#### NSLC_NCI patients will be added to the nonsmoking_samples list
NSLC_NCI_patients = unique(maf_list$NCI$Tumor_Sample_Barcode)
good_samples <- good_samples[!good_samples %in% NSLC_NCI_patients]
good_sample_weights <- bio_weights_unblended[
  Unique_Patient_Identifier %in% good_samples,
]

#### smoking samples are any samples with >0 SBS4 signature weights
smoking_samples <- good_sample_weights[SBS4 > 0, Unique_Patient_Identifier]
nonsmoking_samples <- good_sample_weights[SBS4 == 0, Unique_Patient_Identifier]
#### We are confident that these patients are never-smokers, and the publication indicated that they had low smoking signature despite some having a history of secondary smoking.
nonsmoking_samples <- c(nonsmoking_samples, NSLC_NCI_patients)

## set "A169" as nonsmoking samples manually
length(smoking_samples) ## 760
length(nonsmoking_samples) ## 494
smoking_samples <- smoking_samples[smoking_samples != "A169"]
nonsmoking_samples <- c(nonsmoking_samples, "A169")
length(smoking_samples) ## 759
length(nonsmoking_samples) ## 495

#### use clinical data for TGS to split the samples into smokers and nonsmokers
maf_clinical = fread(paste0(location_output, 'merged_final.txt'))
panel_smoking_samples = unique(maf_clinical[
  Source %in% c('MSK2017', 'MSK2018')
][Smoker == T, `Sample ID`])
panel_nonsmoking_samples = unique(maf_clinical[
  Source %in% c('MSK2017', 'MSK2018')
][Smoker == F, `Sample ID`])
length(sample_with_clinInfo) ## 9230


sample_smo_wP_forCesa <- cesa_smoking_w_panel$samples[
  Unique_Patient_Identifier %in%
    c(
      intersect(smoking_samples, sample_with_clinInfo),
      intersect(panel_smoking_samples, sample_with_clinInfo)
    ),
  Unique_Patient_Identifier
] ## 1066
sample_nonsmo_wP_forCesa <- cesa_nonsmoking_w_panel$samples[
  Unique_Patient_Identifier %in%
    c(
      intersect(nonsmoking_samples, sample_with_clinInfo),
      intersect(panel_nonsmoking_samples, sample_with_clinInfo)
    ),
  Unique_Patient_Identifier
] ## 656

sample_pan_forCsea <- cesa$samples[
  Unique_Patient_Identifier %in% sample_with_clinInfo,
  Unique_Patient_Identifier
] ## 9230

### calculate trinucleotide mutation rates for smoker####
cesa_smoking_w_panel <- trinuc_mutation_rates(
  cesa_smoking_w_panel,
  signature_set = ces.refset.hg19$signatures$COSMIC_v3.2,
  signature_exclusions = treated_samples_signature_exclusions,
  samples = intersect(
    c(
      all_treated_samples_for_sigs,
      all_treatment_not_indicated_samples_for_sigs
    ),
    sample_smo_wP_forCesa
  ),
  cores = 4
)
cesa_smoking_w_panel <- trinuc_mutation_rates(
  cesa_smoking_w_panel,
  signature_set = ces.refset.hg19$signatures$COSMIC_v3.2,
  signature_exclusions = untreated_samples_signature_exclusions,
  samples = intersect(all_untreated_samples_for_sigs, sample_smo_wP_forCesa),
  cores = 4
)


### calculate trinucleotide mutation rates for never-smoker####
cesa_nonsmoking_w_panel <- trinuc_mutation_rates(
  cesa_nonsmoking_w_panel,
  signature_set = ces.refset.hg19$signatures$COSMIC_v3.2,
  signature_exclusions = treated_samples_signature_exclusions,
  samples = intersect(
    c(
      all_treated_samples_for_sigs,
      all_treatment_not_indicated_samples_for_sigs
    ),
    sample_nonsmo_wP_forCesa
  ),
  cores = 4
)
cesa_nonsmoking_w_panel <- trinuc_mutation_rates(
  cesa_nonsmoking_w_panel,
  signature_set = ces.refset.hg19$signatures$COSMIC_v3.2,
  signature_exclusions = untreated_samples_signature_exclusions,
  samples = intersect(all_untreated_samples_for_sigs, sample_nonsmo_wP_forCesa),
  cores = 4
)


### calculate gene mutation rates for smokers ####
cesa_smoking_w_panel = gene_mutation_rates(
  cesa_smoking_w_panel,
  covariates = 'lung',
  samples = sample_smo_wP_forCesa
)


### calculate gene mutation rates for never-smokers ####
cesa_nonsmoking_w_panel = gene_mutation_rates(
  cesa_nonsmoking_w_panel,
  covariates = 'lung',
  samples = sample_nonsmo_wP_forCesa
)

### calculate gene mutation rates for pan data ####
cesa <- gene_mutation_rates(
  cesa,
  covariates = "lung",
  samples = sample_pan_forCsea
)


## calculate cancer effect size ####
### calculate cancer effect size for smokers #####
cesa_smoking_w_panel <- ces_variant(
  cesa = cesa_smoking_w_panel,
  run_name = "recurrents",
  samples = sample_smo_wP_forCesa
)
length(unique(cesa_smoking_w_panel$samples[
  Unique_Patient_Identifier %in% sample_smo_wP_forCesa,
  Unique_Patient_Identifier
])) ## 1066 samples
table(
  cesa_smoking_w_panel$samples[
    Unique_Patient_Identifier %in% sample_smo_wP_forCesa
  ]$coverage,
  exclude = NULL
)
# exome   genome targeted
# 739       20      307
### calculate cancer effect size for never-smokers #####
cesa_nonsmoking_w_panel <- ces_variant(
  cesa = cesa_nonsmoking_w_panel,
  run_name = "recurrents",
  samples = sample_nonsmo_wP_forCesa
)
length(unique(cesa_nonsmoking_w_panel$samples[
  Unique_Patient_Identifier %in% sample_nonsmo_wP_forCesa,
  Unique_Patient_Identifier
])) ## 656 samples
table(
  cesa_nonsmoking_w_panel$samples[
    Unique_Patient_Identifier %in% sample_nonsmo_wP_forCesa
  ]$coverage,
  exclude = NULL
)
# exome   genome targeted
# 304      190      162

### calculate cancer effect size for all samples #####
cesa <- ces_variant(
  cesa = cesa,
  run_name = "recurrents",
  samples = sample_pan_forCsea
) ## on

save_cesa(cesa_smoking_w_panel, paste0(rdata_output, "cesa_smoking.rds"))
save_cesa(cesa_nonsmoking_w_panel, paste0(rdata_output, "cesa_nonsmoking.rds"))
# cesa_smoking_w_panel <- load_cesa(paste0(rdata_output, "cesa_smoking.rds"))
# cesa_nonsmoking_w_panel <- load_cesa(paste0(rdata_output, "cesa_nonsmoking.rds"))

## Fig 1A: only NS, smoker, fill by group; use bar instead of point; pink & brown is better####
oncogene_variants <- c("KRAS_G12D", "KRAS_G12C", "BRAF_V600E", "EGFR_L858R")
target_effect <- rbind(
  cesa_smoking_w_panel$selection$recurrent[
    oncogene_variants,
    on = "variant_name"
  ][, group := "Smoker"][],
  cesa_nonsmoking_w_panel$selection$recurrents[
    oncogene_variants,
    on = "variant_name"
  ][, group := "Never-smoker"][]
)
target_effect[,
  group := factor(group, levels = c("Never-smoker", "Smoker"), ordered = T)
]
variant_order <- target_effect["Never-smoker", on = "group"][order(
  selection_intensity,
  decreasing = T
)][, variant_name]
target_effect[,
  variant_name := factor(variant_name, levels = variant_order, ordered = T)
]
target_effect[,
  prevalence := round(included_with_variant / (included_total + held_out), 3)
]


# svg(file = "p_ces_variant.svg", height = 11, width = 8)
# variant_label <- target_effect[,
#   .(max_y = max(selection_intensity)),
#   by = variant_name
# ]
# ggplot(
#   target_effect,
#   aes(x = variant_name, y = selection_intensity, fill = group)
# ) +
#   geom_bar(
#     stat = "identity",
#     color = "black",
#     position = position_dodge(width = 0.9)
#   ) + ##"dodge"
#   geom_errorbar(
#     aes(ymin = ci_low_95, ymax = ci_high_95),
#     color = "black",
#     na.rm = T,
#     width = 0.3,
#     linewidth = 0.6,
#     position = position_dodge(width = 0.9)
#   ) + ##position_dodge(width = 0.5)
#   # Add variant_name labels above the bars
#   geom_text(
#     data = variant_label,
#     aes(x = variant_name, y = max(max_y) * 3, label = variant_name),
#     inherit.aes = F,
#     vjust = 0,
#     size = 5,
#     fontface = "bold"
#   ) +
#   scale_x_discrete(labels = rep(c("Never-smoker", "Smoker"), length.out = 8)) +
#   scale_y_continuous(
#     trans = scales::pseudo_log_trans(base = 10),
#     breaks = c(0.01, 1, 100, 10000, 1000000),
#     labels = scales::label_comma(accuracy = 1),
#     expand = expansion(mult = c(0.05, 0.05))
#   ) +
#   labs(x = "", y = "Cancer effect size") +
#   scale_fill_manual(values = c("Smoker" = "#D2B48C", "Never-smoker" = "pink")) +
#   theme_classic() +
#   theme(
#     legend.direction = "horizontal",
#     legend.position = "bottom",
#     legend.box = "vertical",
#     legend.text = element_text(size = 18),
#     legend.title = element_text(size = 18),
#     axis.text.x = element_blank(),
#     axis.text.y = element_text(size = 18, face = "bold"),
#     axis.title.y = element_text(size = 18, face = "bold"),
#     axis.title.x = element_text(size = 18, face = "bold")
#   ) +
#   guides(fill = F)
#
# dev.off()

## Fig 1A: only NS, smoker, fill by group; use bar instead of point; pink & brown is better；do a square-root transformation ####
min(target_effect[!group %in% c("Smo_NS"), ci_low_95]) # 4340.395
max(target_effect[!group %in% c("Smo_NS"), ci_high_95]) #  1212397

variant_label <- target_effect[!group %in% c("Smo_NS")][,
  .(max_y = max(selection_intensity)),
  by = variant_name
]
variant_label[, variant_name2 := gsub("_", " ", variant_name)]
# svg(file = "p_ces_variant.svg", height = 11, width = 8)
#
# # the data pans a large range (4342 to 1212428), undo the log or pseudo-log transformation and instead apply a square-root transformation to reduce visual compression while preserving some of the dynamic range.
# ggplot(
#   target_effect[!group %in% c("Smo_NS")],
#   aes(x = variant_name, y = selection_intensity, fill = group)
# ) +
#   geom_bar(
#     stat = "identity",
#     color = "black",
#     position = position_dodge(width = 0.9)
#   ) + ##"dodge"  这个设置不同 group 之间 SI 点分离的距离
#   geom_errorbar(
#     aes(ymin = ci_low_95, ymax = ci_high_95),
#     color = "black",
#     na.rm = T,
#     width = 0.3,
#     linewidth = 0.6,
#     position = position_dodge(width = 0.9)
#   ) + ##position_dodge(width = 0.5)
#   # Add variant_name labels above the bars
#   geom_text(
#     data = variant_label,
#     aes(x = variant_name, y = max(max_y) * 1.4, label = variant_name2),
#     inherit.aes = F,
#     vjust = 0,
#     size = 5,
#     fontface = "bold"
#   ) +
#
#   # Modify x-axis to show "Never-smoker" and "Smoker" below the bars
#   #scale_x_discrete(labels = c("Never-smoker", "Smoker")) +
#   scale_x_discrete(labels = rep(c("Never-smoker", "Smoker"), length.out = 8)) +
#   # scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = F))+ ## expand = expansion(mult = c(0.01,0.05))
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
#   labs(x = "", y = "Cancer effect size") + ## "Somatic variant"
#   scale_fill_manual(values = c("Smoker" = "#D2B48C", "Never-smoker" = "pink")) +
#   theme_classic() + ## theme_minimal()
#   theme(
#     legend.direction = "horizontal",
#     legend.position = "bottom",
#     legend.box = "vertical",
#     legend.text = element_text(size = 18),
#     legend.title = element_text(size = 18),
#     axis.text.x = element_blank(), #element_text( size=18,face = "bold",angle = 30, hjust = 1), ## , angle = 30
#     axis.text.y = element_text(size = 18, face = "bold"),
#     axis.title.y = element_text(size = 18, face = "bold"),
#     axis.title.x = element_text(size = 18, face = "bold")
#   ) +
#   guides(fill = "none")
#
#
# dev.off()

# add tumor growth in mice information -----
## tumor growth vs scaled human CES data; bar plot; dual axis ########
tm_m <- fread(paste0(
  location_data,
  "fig3f_sgInertTumorBurden_webPlotDigitizer_final.csv"
))
# head(tm_m)
tm_m$Genotype <- recode(
  ## to keep consistent with target_effect
  tm_m$Genotype,
  "G12D" = "KRAS_G12D",
  "G12C" = "KRAS_G12C",
  "EGFR" = "EGFR_L858R",
  "BRAF" = "BRAF_V600E"
)
# comput the scaling factor
mean_mouse_growth <- mean(tm_m$TM_M)
mean_ces <- mean(target_effect$selection_intensity)

median_mouse_growth <- median(tm_m$TM_M)
median_ces <- median(target_effect$selection_intensity)

## scale mouse to human
scaling_factor_mean <- mean_ces / mean_mouse_growth ## 0.01642605
scaling_factor_median <- median_ces / median_mouse_growth ## 0.01053374. the scale factor is too small, the mouse bar will appear tiny. so we do the reverse way, scale human to mouse

## scale human to mouse
scaling_factor_mean_rs <- mean_mouse_growth / mean_ces ## 60.87892
scaling_factor_median_rs <- median_mouse_growth / median_ces ## 94.93306

# create a human_Scaled data frame, here is "a global median-matching across variants, not per-variant scaling."
target_effect_sim <- target_effect[, .(
  variant_name,
  group,
  selection_intensity,
  ci_low_95,
  ci_high_95
)]
setnames(
  target_effect_sim,
  old = c("selection_intensity", "ci_low_95", "ci_high_95"),
  new = c("selection_intensity_ori", "ci_low_95_ori", "ci_high_95_ori")
)
human_scaled_median <- target_effect_sim %>%
  mutate(
    group = paste0(group, "_scaledToMouse"),
    selection_intensity = selection_intensity_ori * scaling_factor_median_rs,
    ci_low_95 = ci_low_95_ori * scaling_factor_median_rs,
    ci_high_95 = ci_high_95_ori * scaling_factor_median_rs
  )

human_scaled_mean <- target_effect_sim %>%
  mutate(
    group = paste0(group, "_scaledToMouse"),
    selection_intensity = selection_intensity_ori * scaling_factor_mean_rs,
    ci_low_95 = ci_low_95_ori * scaling_factor_mean_rs,
    ci_high_95 = ci_high_95_ori * scaling_factor_mean_rs
  )


# mouse raw data
mouse_raw_median <- data.frame(
  group = "Mouse_raw",
  selection_intensity = tapply(tm_m$TM_M, tm_m$Genotype, median),
  ci_low_95 = tapply(tm_m$TM_M, tm_m$Genotype, function(x) quantile(x, 0.025)),
  ci_high_95 = tapply(tm_m$TM_M, tm_m$Genotype, function(x) quantile(x, 0.975))
)
mouse_raw_median$variant_name <- rownames(mouse_raw_median)

mouse_raw_mean <- data.frame(
  group = "Mouse_raw",
  selection_intensity = tapply(tm_m$TM_M, tm_m$Genotype, mean),
  ci_low_95 = tapply(tm_m$TM_M, tm_m$Genotype, function(x) quantile(x, 0.025)),
  ci_high_95 = tapply(tm_m$TM_M, tm_m$Genotype, function(x) quantile(x, 0.975))
)
mouse_raw_mean$variant_name <- rownames(mouse_raw_mean)

# # mouse scaled data
# # prefer to use mean value to do the scaling as human CES valus are from likelihood models (model-derived means)
# # but if use tm_mean, then the order is not consistent with the original paper. the original paper use median value as tumor burden are skewed distributions
# mouse_scaled_mean <- data.frame(
#   #variant_name = unique(tm_m$Genotype),  # e.g., "G12D"
#   group = "Mouse_mean", # new group label
#   selection_intensity = tapply(tm_m$TM_M, tm_m$Genotype, mean) *
#     scaling_factor_mean,
#   ci_low_95 = tapply(tm_m$TM_M, tm_m$Genotype, function(x) quantile(x, 0.025)) *
#     scaling_factor_mean,
#   ci_high_95 = tapply(
#     tm_m$TM_M,
#     tm_m$Genotype,
#     function(x) quantile(x, 0.975)
#   ) *
#     scaling_factor_mean
# )
# mouse_scaled_mean$variant_name <- rownames(mouse_scaled_mean)
#
# # median
# mouse_scaled_median <- data.frame(
#   #variant_name = unique(tm_m$Genotype),  # e.g., "G12D"
#   group = "Mouse_median", # new group label
#   selection_intensity = tapply(tm_m$TM_M, tm_m$Genotype, median) *
#     scaling_factor_median,
#   ci_low_95 = tapply(tm_m$TM_M, tm_m$Genotype, function(x) quantile(x, 0.025)) *
#     scaling_factor_median,
#   ci_high_95 = tapply(
#     tm_m$TM_M,
#     tm_m$Genotype,
#     function(x) quantile(x, 0.975)
#   ) *
#     scaling_factor_median
# )
# mouse_scaled_median$variant_name <- rownames(mouse_scaled_median)

# combine scaled data
combined_df_median <- rbind(
  human_scaled_median[, c(
    "group",
    "selection_intensity",
    "ci_low_95",
    "ci_high_95",
    "variant_name"
  )],
  mouse_raw_median
)
combined_df_mean <- rbind(
  human_scaled_mean[, c(
    "group",
    "selection_intensity",
    "ci_low_95",
    "ci_high_95",
    "variant_name"
  )],
  mouse_raw_mean
)

combined_df_median$variant_name <- factor(
  combined_df_median$variant_name,
  levels = c("KRAS_G12D", "BRAF_V600E", "KRAS_G12C", "EGFR_L858R")
)
combined_df_mean$variant_name <- factor(
  combined_df_mean$variant_name,
  levels = c("KRAS_G12D", "BRAF_V600E", "KRAS_G12C", "EGFR_L858R")
)
combined_df_median$variant_name_label <- stringr::str_replace(
  combined_df_median$variant_name,
  "_",
  " "
)
combined_df_median$variant_name_label <- factor(
  combined_df_median$variant_name_label,
  levels = c("KRAS G12D", "BRAF V600E", "KRAS G12C", "EGFR L858R")
)
combined_df_mean$variant_name_label <- stringr::str_replace(
  combined_df_mean$variant_name,
  "_",
  " "
)
combined_df_mean$variant_name_label <- factor(
  combined_df_mean$variant_name_label,
  levels = c("KRAS G12D", "BRAF V600E", "KRAS G12C", "EGFR L858R")
)

# str(combined_df_median)
# plot
library(ggplot2)

plot_dual_axis <- function(
  df,
  title_text = NULL,
  reverse_scaling_factor = reverse_scaling_factor
) {
  ggplot(
    df,
    aes(x = variant_name_label, y = selection_intensity, fill = group)
  ) +
    geom_bar(
      stat = "identity",
      # need to explicitly set bar width <= dodge width
      width = 0.8, ## if not explicitly assign, then the default value is 0.9, which is larger than the dodge width, leading to overlapping bars
      position = position_dodge(width = 0.8), # , preserve = "single"
      color = "black"
    ) +
    geom_errorbar(
      aes(ymin = ci_low_95, ymax = ci_high_95),
      width = 0.3,
      position = position_dodge(width = 0.8), # , preserve = "single"
      color = "black"
    ) +
    scale_y_continuous(
      name = "Tumor growth in mice (neoplastic cells / 1e5 ifu)",
      sec.axis = sec_axis(
        trans = ~ . / reverse_scaling_factor, ## scaled back to human CES
        labels = function(x) format(x, big.mark = ",", scientific = F),
        name = "Scaled selection coefficient" #,
        # breaks = c(1e4, 1e5, 5e5, 1e6, 2e6),
        # labels = label_comma(accuracy = 1),
        # trans = pseudo_log_trans(base = 10)
      ),
      # trans = scales::pseudo_log_trans(base = 10),
      # breaks = c(1e4, 1e5, 5e5, 1e6, 2e6),
      # labels = label_comma(accuracy = 1),
      labels = function(x) format(x, big.mark = ",", scientific = F) #,
      # expand = expansion(mult = c(0.05, 0.1))
    ) +
    scale_fill_manual(
      name = "Group",
      labels = c("Mice", "Never-smoker", "Smoker"),
      values = c(
        "Smoker_scaledToMouse" = "#D2B48C",
        "Never-smoker_scaledToMouse" = "pink",
        "Mouse_raw" = "#999999"
      )
    ) +
    labs(x = "", fill = "Group", title = title_text) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 18, face = "bold"), # angle = 45,hjust = 1
      axis.text.y = element_text(size = 18),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.title.y.right = element_text(
        size = 18,
        face = "bold" #,
        #color = "black"
      ),
      axis.title.x = element_text(size = 18, face = "bold"),
      legend.position = "none"
    )
}

svg(
  file = "p_MicetTumorGrowth_Median_vs_CES_scaled_v9.svg",
  height = 8,
  width = 12
)
plot_dual_axis(
  combined_df_median,
  # "Mouse (median) vs Human (CES scaled to mouse)",
  reverse_scaling_factor = scaling_factor_median_rs
)
dev.off()
# svg can not directly be inserted into google doc, so use jpeg instead
jpeg(
  file = "p_MicetTumorGrowth_Median_vs_CES_scaled_v9.jpeg",
  height = 8,
  width = 12,
  units = 'in',
  res = 600
)
plot_dual_axis(
  combined_df_median,
  # "Mouse (median) vs Human (CES scaled to mouse)",
  reverse_scaling_factor = scaling_factor_median_rs
)
dev.off()
svg(
  file = "p_MicetTumorGrowth_Mean_vs_CES_scaled_v9.svg",
  height = 8,
  width = 12
)
plot_dual_axis(
  combined_df_mean,
  # "Mouse (mean) vs Human (CES scaled to mouse)",
  reverse_scaling_factor = scaling_factor_mean_rs
)
dev.off()
# svg can not directly be inserted into google doc, so use jpeg instead
jpeg(
  file = "p_MicetTumorGrowth_Mean_vs_CES_scaled_v9.jpeg",
  height = 8,
  width = 12,
  units = 'in',
  res = 600
)
plot_dual_axis(
  combined_df_mean,
  # "Mouse (median) vs Human (CES scaled to mouse)",
  reverse_scaling_factor = scaling_factor_mean_rs
)
dev.off()

## scaled tumor growth vs human CES data; scatter plot #####
# scatter plot: compare relative tumor size/number with CES_B_on_A-----
draw_scatter_forCES <- function(
  human_data,
  MiceTumorData,
  MiceTumorData_type = "size", ## "size" or "number"

  #MiceTumorlimit,
  # data_type,
  max_overlaps = 10 ## default
) {
  ## MiceTumorlimit = c(0.5,32)
  library(ggplot2)
  #variants <- deparse(substitute(variant_onco)) # Get the name of the variable passed as variant_onco
  human_TS_data <- dplyr::inner_join(
    human_data,
    MiceTumorData,
    by = c("variant_name")
  )

  pearson_corr <- cor.test(
    human_TS_data$tumorGrowth,
    human_TS_data$selection_intensity,
    method = c("pearson")
  ) ## can only assign one method each time
  pearson_r <- pearson_corr$estimate
  pearson_p <- pearson_corr$p.value

  spearman_corr <- cor.test(
    human_TS_data$tumorGrowth,
    human_TS_data$selection_intensity,
    method = c("spearman")
  ) ## can only assign one method each time
  spearman_r <- spearman_corr$estimate
  spearman_p <- spearman_corr$p.value

  pp_wo_smooth <- ggplot(
    data = human_TS_data[order(-selection_intensity)],
    aes(x = selection_intensity, y = tumorGrowth)
  ) + ##reorder(variant_B, log10(ces_B_on_A))

    labs(
      y = "Tumor growth in mice (neoplastic cells / 1e5 ifu)"
    ) +
    xlab(bquote(
      ~Cancer ~ Effect ~ Size ~ scriptstyle(~ ~ (log[10]))
    )) +
    scale_x_log10(
      labels = function(x) format(x, big.mark = ",", scientific = F)
    ) +
    scale_y_log10(
      labels = function(x) format(x, big.mark = ",", scientific = F)
    ) +
    #scale_y_continuous(trans = "log2", limits = MiceTumorlimit) + #tgutil::scale_y_log2()+
    geom_errorbar(
      data = human_TS_data[
        !is.na(human_TS_data$ci_low_95_TMm) &
          !is.na(human_TS_data$ci_high_95_TMm),
      ],
      aes(ymin = ci_low_95_TMm, ymax = ci_high_95_TMm),
      width = 0,
      color = "grey"
    ) +
    geom_errorbar(
      data = human_TS_data[
        !is.na(human_TS_data$ci_low_95) &
          !is.na(human_TS_data$ci_high_95),
      ],
      aes(xmin = ci_low_95, xmax = ci_high_95),
      width = 0,
      color = "grey"
    ) +
    geom_point(aes(color = group, shape = group), size = 4) + ## position = "jitter", size = 1
    #geom_smooth(method = "lm", se = T)+
    ggrepel::geom_label_repel(
      aes(label = variant_name, color = variant_name), # fill = group,
      max.overlaps = max_overlaps,
      # text
      #fill = "grey",        # fill in the rectangle
      # segment.colour = "grey",# border + connecting line

      size = 3
    ) +
    theme_classic() +
    #theme(legend.position = "none") +
    guides(color = "none", label = "none")
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

  return(
    result = list(
      pp_wo_smooth,
      human_TS_data,
      pearson_r,
      pearson_p,
      spearman_r,
      spearman_p
    )
  )
}
target_effect_forScatter <- target_effect[, .(
  variant_name,
  group,
  selection_intensity,
  ci_low_95,
  ci_high_95
)]

mouse_raw_median_forScatter <- data.frame(
  tumorGrowth = tapply(tm_m$TM_M, tm_m$Genotype, median),
  ci_low_95_TMm = tapply(
    tm_m$TM_M,
    tm_m$Genotype,
    function(x) quantile(x, 0.025)
  ),
  ci_high_95_TMm = tapply(
    tm_m$TM_M,
    tm_m$Genotype,
    function(x) quantile(x, 0.975)
  )
)
mouse_raw_median_forScatter$variant_name <- rownames(
  mouse_raw_median_forScatter
)

svg(
  file = "p_MicetTumorGrowth_Median_vs_CES_scatter.svg",
  height = 11,
  width = 8
)
p_scatter <- draw_scatter_forCES(
  target_effect_forScatter,
  mouse_raw_median_forScatter
)
p_scatter[[1]]
dev.off()

## calculate the Pearson and Spearman correlations separately for smoker and never-smoker ####

compute_corr <- function(df, label) {
  pearson_corr <- cor.test(
    df$tumorGrowth,
    df$selection_intensity,
    method = "pearson"
  )

  spearman_corr <- cor.test(
    df$tumorGrowth,
    df$selection_intensity,
    method = "spearman"
  )

  data.frame(
    group = label,
    n = nrow(df),
    pearson_r = unname(pearson_corr$estimate),
    pearson_p = pearson_corr$p.value,
    spearman_r = unname(spearman_corr$estimate),
    spearman_p = spearman_corr$p.value
  )
}

human_TS_data <- p_scatter[[2]]
corr_results <- rbind(
  compute_corr(
    subset(human_TS_data, group == "Smoker"),
    "Smoker"
  ),
  compute_corr(
    subset(human_TS_data, group == "Never-smoker"),
    "Never-smoker"
  ),
  compute_corr(
    human_TS_data,
    "NS_S"
  )
)

corr_results

## for testing
cesa$coverage_ranges$exome$exome ## get the default exome covered region
