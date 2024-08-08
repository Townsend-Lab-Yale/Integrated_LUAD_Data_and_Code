if(!require("cancereffectsizeR")){
  ## use cancereffectsizeR v2.8.1 here
  remotes::install_github("Townsend-Lab-Yale/cancereffectsizeR@v2.8.1", dependencies = T, force = T) 
  require("cancereffectsizeR")
}
if(!require("ces.refset.hg19")){
  ## use ces.refset.hg19 as the reference dataset here
  remotes::install_github("Townsend-Lab-Yale/ces.refset.hg19@*release", dependencies = T, force = T) 
  require("ces.refset.hg19")
}
##  import other libraries
library(data.table)
library(dplyr)
library(rtracklayer)
library(stringr)
library(ggplot2)


location_output <- "../integrated_data/"
location_data <- "../data/"
## read maf file
maf_file <- read.csv(paste0(location_output,"merged_luad_maf.txt"))
colnames(maf_file)[2] <- 'Tumor_Sample_Barcode'
colnames(maf_file)[7] <- 'Tumor_Allele'
maf_list <- split(maf_file, maf_file$Source) 

liftover_file = paste0(location_data, "hg38ToHg19.over.chain")

## preload all maf file ####
### WES+WGS
Broad_maf <- preload_maf(maf = maf_list$Broad, refset = ces.refset.hg19, chain_file = liftover_file, keep_extra_columns = T)
Broad_maf <- Broad_maf[is.na(problem)]
Broad_maf <- Broad_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]
Broad_maf <- Broad_maf[!Unique_Patient_Identifier %in% c("LUAD-B01169","LUAD-D01382")]
### TGS data
FMAD_maf <- preload_maf(maf = maf_list$`FM-AD`, ces.refset.hg19, chain_file = liftover_file, keep_extra_columns = T)
FMAD_maf <- FMAD_maf[is.na(problem)]

Genie_maf <- preload_maf(maf = maf_list$Genie, refset = ces.refset.hg19, chain_file = liftover_file, keep_extra_columns = T)
Genie_maf <- Genie_maf[is.na(problem)]

### WES
MSK2015_maf <- preload_maf(maf = maf_list$MSK2015, refset = ces.refset.hg19, chain_file = liftover_file, keep_extra_columns = T)
MSK2015_maf <- MSK2015_maf[is.na(problem)]
MSK2015_maf <- MSK2015_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

### TGS data
MSK2017_maf <- preload_maf(maf = maf_list$MSK2017, refset = ces.refset.hg19, chain_file = liftover_file, keep_extra_columns = T)
MSK2017_maf <- MSK2017_maf[is.na(problem)]

### WES; TGS
MSK2018_maf <- preload_maf(maf = maf_list$MSK2018, refset = ces.refset.hg19, chain_file = liftover_file, keep_extra_columns = T)
MSK2018_maf <- MSK2018_maf[is.na(problem)]

### WES data
OncoSG_maf <- preload_maf(maf = maf_list$OncoSG, refset = ces.refset.hg19, chain_file = liftover_file, keep_extra_columns = T)
OncoSG_maf <- OncoSG_maf[is.na(problem)]
OncoSG_maf <- OncoSG_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

### WES data
TCGA_maf <- preload_maf(maf = maf_list$TCGA, refset = ces.refset.hg19, chain_file = liftover_file, keep_extra_columns = T)
TCGA_maf <- TCGA_maf[is.na(problem)]
TCGA_maf <- TCGA_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

### WES data
TracerX_maf <- preload_maf(maf = maf_list$TracerX, refset = ces.refset.hg19, chain_file = liftover_file, keep_extra_columns = T)
TracerX_maf <- TracerX_maf[is.na(problem)]
TracerX_maf <- TracerX_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]

### TGS data
TSP_maf <- preload_maf(maf = maf_list$TSP, refset = ces.refset.hg19, chain_file = liftover_file, keep_extra_columns = T)
TSP_maf <- TSP_maf[is.na(problem)]
TSP_maf <- TSP_maf[!Unique_Patient_Identifier %in% c("luad_tsp_16929", "luad_tsp_16901", "luad_tsp_16875","luad_tsp_16915")]

### WGS data
NCI_maf <- preload_maf(maf = maf_list$NCI, refset = ces.refset.hg19, chain_file = liftover_file, keep_extra_columns = T)
NCI_maf <- NCI_maf[is.na(problem)]
NCI_maf <- NCI_maf[germline_variant_site == F & (repetitive_region == F | cosmic_site_tier %in% 1:3)]


## extract necessary coverage info
### read in information nessary for loading maf files (exome/genome and coverage intervals for TGS)
broad_exome_or_genome <- fread(paste0(location_data,"luad_broad/data_clinical_sample.txt"))[-(1:4),c('Sample Identifier', 'Platform')]
broad_exome_or_genome$Platform <- sapply(broad_exome_or_genome$Platform, function(x){if(str_detect(x, 'WGS')){return('WGS')} else return('WES')})
msk2017_panels_used <- fread(paste0(location_data,"lung_msk_2017/data_clinical_sample.txt"))[-(1:4),c('Sample Identifier', 'Gene Panel')]
msk2018_panels_used <- fread(paste0(location_data,"nsclc_pd1_msk_2018/data_clinical_sample.txt"))[-(1:4),c('Sample Identifier', 'Gene Panel')]
genie_panels_used <- fread(paste0(location_data,"genie_9/data_clinical_sample.txt"))[-(1:4),c('Sample Identifier', 'Sequence Assay ID')]

### read in genes included in each panel and creating GRANGES object (bed files)  to pass into covered_regions parameter of load_maf 
### once the granges are exported once, these functions don't need to be run anymore
# gene_granges <- rtracklayer::import(paste0(location_data, "gencode.v38lift37.basic.annotation.gtf"))
load("gencode.v38lift37.basic.annotation.gtf.Rdata")

location_bed <- paste0(location_data,'bed_files/')
if(!dir.exists(location_bed)){
  dir.create(location_bed)
}

location_gene_panels <- paste0(location_data,'gene_panels/')

if(!file.exists(paste0(location_bed,"fmad_targets.bed"))){
  fmad_genes <- unique(fread(paste0(location_gene_panels,"foundation_one.txt"))$Hugo_Symbol)
  fmad_granges <- gene_granges[gene_granges$gene_name %in% fmad_genes, ]
  fmad_granges <- fmad_granges[fmad_granges$type %in% c('CDS','stop_codon'),]
  fmad_gr_clean <- cancereffectsizeR:::clean_granges_for_cesa(cesa = CESAnalysis(ces.refset.hg19), gr = fmad_granges)
  export(fmad_gr_clean, paste0(location_bed,"fmad_targets.bed"))
}

if(!file.exists(paste0(location_bed,"tsp_targets.bed"))){
  tsp_genes <- unique(fread(paste0(location_gene_panels,"tsp.txt"))$Hugo_Symbol)
  tsp_granges <- gene_granges[gene_granges$gene_name %in% tsp_genes, ]
  tsp_granges <- tsp_granges[tsp_granges$type %in% c('CDS','stop_codon'),]
  tsp_gr_clean <- cancereffectsizeR:::clean_granges_for_cesa(cesa = CESAnalysis(ces.refset.hg19), gr = tsp_granges)
  export(tsp_gr_clean, paste0(location_bed,"tsp_targets.bed"))
}

if(!file.exists(paste0(location_bed,"msk341_targets.bed"))){
  msk_341_genes <- unique(fread(paste0(location_gene_panels,"msk341.txt"))$Hugo_Symbol)
  msk341_granges <- gene_granges[gene_granges$gene_name %in% msk_341_genes, ]
  msk341_granges <- msk341_granges[msk341_granges$type %in% c('CDS','stop_codon'),]
  msk341_gr_clean <- cancereffectsizeR:::clean_granges_for_cesa(cesa = CESAnalysis(ces.refset.hg19), gr = msk341_granges)
  export(msk341_gr_clean, paste0(location_bed,"msk341_targets.bed"))
}

if(!file.exists(paste0(location_bed,"msk410_targets.bed"))){
  msk_410_genes <- unique(fread(paste0(location_gene_panels,"msk410.txt"))$Hugo_Symbol)
  msk410_granges <- gene_granges[gene_granges$gene_name %in% msk_410_genes, ]
  msk410_granges <- msk410_granges[msk410_granges$type %in% c('CDS','stop_codon'),]
  msk410_gr_clean <- cancereffectsizeR:::clean_granges_for_cesa(cesa = CESAnalysis(ces.refset.hg19), gr = msk410_granges)
  export(msk410_gr_clean, paste0(location_bed,"msk410_targets.bed"))
}

if(!file.exists(paste0(location_bed,"msk468_targets.bed"))){
  msk_468_genes <- unique(fread(paste0(location_gene_panels,"msk468.txt"))$Hugo_Symbol)
  msk468_granges <- gene_granges[gene_granges$gene_name %in% msk_468_genes, ]
  msk468_granges <- msk468_granges[msk468_granges$type %in% c('CDS','stop_codon'),]
  msk468_gr_clean <- cancereffectsizeR:::clean_granges_for_cesa(cesa = CESAnalysis(ces.refset.hg19), gr = msk468_granges)
  export(msk468_gr_clean, paste0(location_bed,"msk468_targets.bed"))
}

### this should be run everytime so that we don't have to save a bunch of grange files
genie_panel_genes <- fread(paste0(location_data,"genie_9/genomic_information.txt"))[,c('Chromosome', 'Start_Position', 'End_Position', 'Hugo_Symbol', 'Feature_Type', 'SEQ_ASSAY_ID')]
genie_granges_list <- makeGRangesListFromDataFrame(genie_panel_genes, ignore.strand = T, seqnames.field = 'Chromosome', start.field = 'Start_Position', end.field = 'End_Position', split.field = 'SEQ_ASSAY_ID')
seqlevels(genie_granges_list, pruning.mode = "fine") <- c(1:22,'X','Y')

### some panels in GENIE do not cover TP53 or KRAS so they must be removed
genie_panel_genes_2 <- genie_panel_genes[,c('Hugo_Symbol','SEQ_ASSAY_ID')]
genie_panel_genes_2 <- genie_panel_genes_2[!duplicated(genie_panel_genes_2)]
genie_panel_genes_list <- split(genie_panel_genes_2$Hugo_Symbol, genie_panel_genes_2$SEQ_ASSAY_ID)

Genie_maf <- merge(Genie_maf, genie_panels_used, by.x = 'Unique_Patient_Identifier', by.y = 'Sample Identifier')

panels_to_remove <- c()
for(panel in names(genie_panel_genes_list)){
  if(!('TP53' %in% genie_panel_genes_list[[panel]]) | !('KRAS' %in% genie_panel_genes_list[[panel]])){
    panels_to_remove <- c(panels_to_remove, panel)
  }
}

Genie_maf <- Genie_maf[!'Sequence Assay ID' %in% panels_to_remove]

### split maf files when different sequencing methods are used
Broad_maf <- merge(Broad_maf, broad_exome_or_genome, by.x = 'Unique_Patient_Identifier', by.y = 'Sample Identifier')
Broad_maf <- split(Broad_maf, Broad_maf$Platform)

MSK2017_maf <- merge(MSK2017_maf, msk2017_panels_used, by.x = 'Unique_Patient_Identifier', by.y = 'Sample Identifier')
MSK2017_maf <- split(MSK2017_maf, MSK2017_maf$`Gene Panel`)

MSK2018_maf <- merge(MSK2018_maf, msk2018_panels_used, by.x = 'Unique_Patient_Identifier', by.y = 'Sample Identifier')
MSK2018_maf <- split(MSK2018_maf, MSK2018_maf$`Gene Panel`)

Genie_maf <- merge(Genie_maf, genie_panels_used, by.x = 'Unique_Patient_Identifier', by.y = 'Sample Identifier')
Genie_maf <- split(Genie_maf, Genie_maf$`Sequence Assay ID.x`)


## load maf files into cesa ojbect ####
cesa <- CESAnalysis(ces.refset.hg19)
### consider using covered regions padding when variants are outside the intervals
#### WES data
cesa <- load_maf(cesa, maf = Broad_maf$WES)
cesa <- load_maf(cesa, maf = MSK2015_maf)
cesa <- load_maf(cesa, maf = OncoSG_maf)
cesa <- load_maf(cesa, maf = TCGA_maf)
cesa <- load_maf(cesa, maf = TracerX_maf)

#### TGS data
cesa <- load_maf(cesa, maf = FMAD_maf, coverage = 'targeted', 
                 covered_regions = paste0(location_bed,"fmad_targets.bed"), 
                 covered_regions_name = 'fmad_regions', covered_regions_padding = 100) #padding based on 23 variants having distance from interval between 10 and 100.

for(i in 1:length(Genie_maf)){
  cesa <- load_maf(cesa, maf = Genie_maf[i][[1]], coverage = 'targeted', 
                   covered_regions = genie_granges_list[names(Genie_maf)[i]][[1]], 
                   covered_regions_name = paste0(names(Genie_maf)[i], '_regions'), covered_regions_padding = 100)
}

cesa <- load_maf(cesa, maf = MSK2017_maf$IMPACT341, coverage = 'targeted', covered_regions = paste0(location_bed,"msk341_targets.bed"), covered_regions_name = 'msk341_regions', covered_regions_padding = 100)
cesa <- load_maf(cesa, maf = MSK2017_maf$IMPACT410, coverage = 'targeted', covered_regions = paste0(location_bed,"msk410_targets.bed"), covered_regions_name = 'msk410_regions', covered_regions_padding = 100)
cesa <- load_maf(cesa, maf = MSK2018_maf$IMPACT341, coverage = 'targeted', covered_regions = paste0(location_bed,"msk341_targets.bed"), covered_regions_name = 'msk341_regions', covered_regions_padding = 100)
cesa <- load_maf(cesa, maf = MSK2018_maf$IMPACT410, coverage = 'targeted', covered_regions = paste0(location_bed,"msk410_targets.bed"), covered_regions_name = 'msk410_regions', covered_regions_padding = 100)
cesa <- load_maf(cesa, maf = MSK2018_maf$IMPACT468, coverage = 'targeted', covered_regions = paste0(location_bed,"msk468_targets.bed"), covered_regions_name = 'msk468_regions', covered_regions_padding = 100)
cesa <- load_maf(cesa, maf = TSP_maf, coverage = 'targeted', covered_regions = paste0(location_bed,"tsp_targets.bed"), covered_regions_name = 'tsp_regions', covered_regions_padding = 100)

#### WGS data
cesa <- load_maf(cesa, maf = Broad_maf$WGS, coverage = 'genome')
cesa <- load_maf(cesa, maf = NCI_maf, coverage = 'genome')

save_cesa(cesa,"load_maf_cesa_WES_TGS_WGS.rds")
save.image()


## calculate mutation rates for all samples, smokers, and never-smokers  ####
### use mutational signature convolution on WES/WGS and clinical data for TGS to split the samples into smokers and nonsmokers ####
signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "LUAD", treatment_naive = F)
cesa <- trinuc_mutation_rates(cesa,
                              signature_set = ces.refset.hg19$signatures$COSMIC_v3.2,
                              signature_exclusions = signature_exclusions,  
                             # samples = cesa$samples[coverage %in% c('exome','genome'), Unique_Patient_Identifier], ## keep TGS samples for trinuc_mutation_rates() as Jeff Mandell recommended
                              cores=4)
bio_weights <- cesa$mutational_signatures$biological_weights
bio_weights_unblended <- bio_weights[bio_weights$group_avg_blended == F]

snv_counts <- cesa$maf[variant_type == 'snv', .N, by = "Unique_Patient_Identifier"]
good_samples <- snv_counts[N > 50, Unique_Patient_Identifier]


#### NSLC_NCI patients will be added to the nonsmoking_samples list
NSLC_NCI_patients = unique(maf_list$NCI$Tumor_Sample_Barcode)
good_samples <- good_samples[! good_samples %in% NSLC_NCI_patients]
good_sample_weights <- bio_weights_unblended[Unique_Patient_Identifier %in% good_samples,]

#### smoking samples are any samples with >0 signature weights 
smoking_samples <- good_sample_weights[SBS4 > 0, Unique_Patient_Identifier]
nonsmoking_samples <- good_sample_weights[SBS4 == 0, Unique_Patient_Identifier]
#### We are confident that these patients are never-smokers, and the publication indicated that they had low smoking signature despite some having a history of secondary smoking.
nonsmoking_samples <- c(nonsmoking_samples,
                        NSLC_NCI_patients)

#### include panel data
maf_clinical = fread(paste0(location_output, 'merged_final.txt'))
panel_smoking_samples = unique(maf_clinical[Source %in% c('MSK2017','MSK2018')][Smoker == T, `Sample ID`])
panel_nonsmoking_samples = unique(maf_clinical[Source %in% c('MSK2017','MSK2018')][Smoker == F, `Sample ID`])

### calculate mutation rates for all samples####
cesa <- gene_mutation_rates(cesa, covariates = "lung")

### calculate mutation rates for smokers ####
#### mutation rates for smokers; WES+WGS+Panel
cesa_smoking_w_panel = clear_gene_rates(cesa)
cesa_smoking_w_panel = gene_mutation_rates(cesa_smoking_w_panel, covariates = 'lung',
                                           samples = cesa_smoking_w_panel$samples[Unique_Patient_Identifier %in% c(smoking_samples, panel_smoking_samples), Unique_Patient_Identifier])
### calculate mutation rates for never-smokers ####
#### mutation rates for smokers; WES+WGS+Panel
cesa_nonsmoking_w_panel = clear_gene_rates(cesa)
cesa_nonsmoking_w_panel = gene_mutation_rates(cesa_nonsmoking_w_panel, covariates = 'lung',
                                              samples = cesa_nonsmoking_w_panel$samples[Unique_Patient_Identifier %in% c(nonsmoking_samples, panel_nonsmoking_samples), Unique_Patient_Identifier])



## calculate cancer effect size ####
### calculate cancer effect size for all samples #####
cesa <- ces_variant(cesa = cesa, run_name = "recurrents")
length(unique(cesa$samples$Unique_Patient_Identifier)) ## 9242 samples
table(cesa$samples$coverage, exclude = NULL)
# exome   genome targeted 
# 1113      211     7918 

### calculate cancer effect size for smokers #####
cesa_smoking_w_panel <- ces_variant(cesa = cesa_smoking_w_panel, run_name = "recurrents",
                                           samples = cesa_smoking_w_panel$samples[Unique_Patient_Identifier %in% c(smoking_samples, panel_smoking_samples), Unique_Patient_Identifier])
length(unique(cesa_smoking_w_panel$samples[Unique_Patient_Identifier %in% c(smoking_samples, panel_smoking_samples), Unique_Patient_Identifier])) ## 1082 samples
table(cesa_smoking_w_panel$samples[Unique_Patient_Identifier %in% c(smoking_samples, panel_smoking_samples)]$coverage, exclude = NULL)
# exome   genome targeted 
# 669       19      394
### calculate cancer effect size for never-smokers #####
cesa_nonsmoking_w_panel <- ces_variant(cesa = cesa_nonsmoking_w_panel, run_name = "recurrents",
                                              samples = cesa_nonsmoking_w_panel$samples[Unique_Patient_Identifier %in% c(nonsmoking_samples, panel_nonsmoking_samples), Unique_Patient_Identifier])
length(unique(cesa_nonsmoking_w_panel$samples[Unique_Patient_Identifier %in% c(nonsmoking_samples, panel_nonsmoking_samples), Unique_Patient_Identifier])) ## 642 samples
table(cesa_nonsmoking_w_panel$samples[Unique_Patient_Identifier %in% c(nonsmoking_samples, panel_nonsmoking_samples)]$coverage, exclude = NULL)
# exome   genome targeted 
# 286      192      164

# save_cesa(cesa, "cesa_allSamples.rds")
# save_cesa(cesa,"cesa_smoking.rds")
# save_cesa(cesa, "cesa_nonsmoking.rds")
load_cesa("cesa_allSamples.rds")
load_cesa("cesa_smoking.rds")
load_cesa("cesa_nonsmoking.rds")

## ploting cancer effect size ####
oncogene_variants <- c("KRAS_G12D","KRAS_G12C","BRAF_V600E","EGFR_L858R")
target_effect <- rbind(cesa$selection$recurrents[oncogene_variants, on = "variant_name"][, group := "All"][],
      cesa_smoking_w_panel$selection$recurrents[oncogene_variants, on = "variant_name"][, group:= "Smoker"][],
      cesa_nonsmoking_w_panel$selection$recurrents[oncogene_variants, on = "variant_name"][,group := "Never-smoker"][])
target_effect[,group := factor(group, levels = c("All","Smoker","Never-smoker"), ordered =T)]
variant_order <- target_effect["All", on = "group"][order(selection_intensity, decreasing = T)][,variant_name]
target_effect[,variant_name:= factor(variant_name, levels = variant_order, ordered =T)]
target_effect[,prevalence := round(included_with_variant/(included_total+held_out),3)]

svg(file = "p_ces_variant_allInOne_facetByGroup.svg", height = 7, width = 9)
ggplot(target_effect, aes(x = variant_name, y = selection_intensity, fill = group)) +
  geom_errorbar(aes(ymin = ci_low_95, ymax = ci_high_95), color = "azure4", na.rm = T, width = 0.4, linewidth = 0.4)+
  geom_point(shape =21, color = "gray20", aes(size =prevalence))+ 
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = F))+
  labs(x = "",  y= "Cancer effect size",size = "Prevalence")+ ## Somatic variant
  scale_fill_manual(values =c("All" = "black", "Smoker" = "red","Never-smoker" = "blue"))+
  facet_wrap(~group, nrow = 1, strip.position = "top")+
  theme_classic()+
  theme(legend.direction = "horizontal",
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.title= element_text(size = 16),
        axis.text.x = element_text( size=10,face = "bold", angle = 30, hjust = 1),
        axis.text.y = element_text(size =16,face = "bold"),
        axis.title.y = element_text( size = 18,face= "bold"),
        axis.title.x = element_text( size = 18,face = "bold"),
        strip.text.x = element_text(angle = 0, hjust =0.5, vjust = 0.5,size =16, face = "bold"),
        panel.spacing = unit(1,"lines"))+ 
  guides(fill = "none")+
  scale_size_area(labels = (scales::label_percent(accuracy = 0.1)), max_size = 8)
dev.off()










