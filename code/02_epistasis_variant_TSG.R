library(ggplot2)
library(ggpubr)
library(patchwork)

## epistasis analysis between a variant and a gene; using select_variants() to include more samples #####
get_epistasis_cesa <- function(cesa_object, variant, gene_name, specific_samples){ 
  variants_geneName <- cancereffectsizeR::select_variants(cesa = cesa_object, genes = gene_name) 
  variants_variant <- cancereffectsizeR::select_variants(cesa_object, variant_ids = variant) 
  variants_target = rbind(variants_variant,variants_geneName)
  comp_target <- define_compound_variants(cesa = cesa_object, variant_table = variants_target, by = "gene", merge_distance = Inf)

  run_name <- paste0(variant,"_vs_",gene_name)
  cesa_object <- ces_epistasis(cesa = cesa_object, variants = comp_target, run_name = run_name, samples = specific_samples)
  epi <- cesa_object$epistasis[[run_name]]
  set(epi, j= "run_type", value = run_name) 
  result <- list(epi, cesa_object, comp_target, variants_target)
  return(result)
}

tsg <- c("PTEN","STK11","Nf1","Setd2","Stag2","Tsc1","Rb1","Fbxw7","Rnf43","Rbm10","Mga","Keap1","Kmt2d","Smad4","Cmtr2","Tp53","Nf2","Apc","Ptprd","Atm","Msh2","Arid2","Kdm6a","Atrx","BRCA2","Palb2")
tsg <- toupper(tsg)
tsg <- c(tsg, "CDKN2A.p14arf","CDKN2A.p16INK4a")

# Pre-allocate epi_final with column names from the first iteration
KRAS_G12D_vs_PTEN <- get_epistasis_cesa(cesa_object = cesa_smoking_w_panel, variant = oncogene_variants[1], gene_name = tsg[1], specific_samples = sample_smo_wP_forCesa)
epi_initial <- KRAS_G12D_vs_PTEN[[1]]
epi_final <-  data.table(matrix(NA, nrow = 0, ncol = length(epi_initial))) ## create an empty data table with 0 rows and length(epi_initial) columns
setnames(epi_final, names(epi_initial)) ## then set the column names
epi_final_smoking_wP <- epi_final
epi_final_nonsmoking_wP <- epi_final


for(i in 1:length(oncogene_variants)){
  for(j in 1:length(tsg)){
    
    epi_smoking_wP <- get_epistasis_cesa(cesa_object = cesa_smoking_w_panel, variant = oncogene_variants[i], gene_name = tsg[j], 
                                         specific_samples = sample_smo_wP_forCesa)[[1]]
    epi_final_smoking_wP <- rbindlist(list(epi_final_smoking_wP, epi_smoking_wP), fill = F)
    
    epi_nonsmoking_wP <- get_epistasis_cesa(cesa_object = cesa_nonsmoking_w_panel, variant = oncogene_variants[i], gene_name = tsg[j], 
                                            specific_samples = sample_nonsmo_wP_forCesa)[[1]]
    epi_final_nonsmoking_wP <- rbindlist(list(epi_final_nonsmoking_wP, epi_nonsmoking_wP), fill = F) 
    }
}



epi_data <- rbind(epi_final_smoking_wP[, data_type := "Smo"],
                  epi_final_nonsmoking_wP[, data_type := "nonSmo"])
epi_data <- epi_data[,c("variant_A_type","remaining_column") := tstrsplit(run_type, split = "_vs", keep = 1)]
epi_data[,gene_Name :=stringr::str_split(run_type, "_", simplify = T)[,4]]

## remove all rows with "CDKN2A.p16INK4a" in "gene_Name" column and "CDKN2A.p14arf" in "variant_B" column   
## because select_variants() function also select all variant information of CDKN2A.p14arf" when the "gene_name" is "CDKN2A.p16INK4a", which cause a duplication.
epi_data <- epi_data[!.("CDKN2A.p16INK4a","CDKN2A.p14arf"), on=.(gene_Name,variant_B)]

## add variant_B_s column, to shorten CDKN2A.p16INK4a and CDKN2A.p14arf to p16INK4a and p14arf
epi_data[,variant_B_s := variant_B]
epi_data[stringr::str_detect(variant_B_s,"\\."),variant_B_s := stringr::str_split(variant_B_s,"\\.", simplify = T)[,2]]

save(epi_data, file = paste0(rdata_output, "epistasis_variant_TSG_output.Rdata"))


### for KRAS_G12C ####
epi_KRAS_G12C <- epi_data["KRAS_G12C", on = "variant_A_type"]
epi_KRAS_G12C_S <- epi_data[.("KRAS_G12C","Smo"), on=.(variant_A_type, data_type)]
epi_KRAS_G12C_nonS <- epi_data[.("KRAS_G12C","nonSmo"), on=.(variant_A_type, data_type)]

## sorting gene by ces_B_on_A of TSG under KRSA_G12C background in NS, then smokers
variant_order_KRAS_G12C <- epi_KRAS_G12C_nonS[ces_B_on_A > 0.001000][order(-ces_B_on_A)]$variant_B_s
variant_lf_KRAS_G12C <- epi_KRAS_G12C_nonS[ces_B_on_A == 0.001000]$variant_B_s
variant_lf_order_KRAS_G12C <- epi_KRAS_G12C_S[variant_lf_KRAS_G12C, on = "variant_B_s"][ces_B_on_A > 0.001000][order(-ces_B_on_A)]$variant_B_s
variant_lf_2_KRAS_G12C <- epi_KRAS_G12C_S[variant_lf_KRAS_G12C, on = "variant_B_s"][ces_B_on_A == 0.001000]$variant_B_s

gene_order_KRAS_G12C <- c(variant_order_KRAS_G12C,variant_lf_order_KRAS_G12C,variant_lf_2_KRAS_G12C)
gene_order_KRAS_G12C

epi_KRAS_G12C[,variant_B_s := factor(variant_B_s, levels = gene_order_KRAS_G12C, ordered = T )] 
epi_KRAS_G12C[,data_type := factor(data_type, levels = c("nonSmo","Smo"), ordered = T )] 

## edit CI bar which has NA values,CI should be zero if it is NA
epi_KRAS_G12C$ci_low_95_ces_B_on_A_forplot <- ifelse(is.na(epi_KRAS_G12C$ci_low_95_ces_B_on_A),0, epi_KRAS_G12C$ci_low_95_ces_B_on_A)
epi_KRAS_G12C$errorBarColor <- ifelse(epi_KRAS_G12C$nAB %in% c("0","1"), "Gray","Black")
epi_KRAS_G12C$data_type2 <-  as.vector(epi_KRAS_G12C$data_type)
epi_KRAS_G12C$data_type2[epi_KRAS_G12C$nAB == "0"] <- "Gray"


p_G12C <- ggplot(data= epi_KRAS_G12C, aes(x= data_type , y = ces_B_on_A))+ 
   geom_errorbar(aes(ymin = ci_low_95_ces_B_on_A_forplot, ymax = ci_high_95_ces_B_on_A, color = errorBarColor), width = 0,  na.rm = T) + 
  geom_point(aes(fill = data_type2, color = errorBarColor),shape = 21, size =2) + 
  scale_fill_manual(values =c("Smo" = "red","nonSmo" = "blue","Gray" = "gray")) + ## "All" = "black", 
  scale_color_manual(values = c("Gray" = "gray","Black" = "black"))+
  labs(x = "Tumor suppressor gene")+
  ylab(bquote(bold(atop("Cancer effect",
                        "(KRAS G12C background)"))))+ 
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = F))+
  scale_x_discrete(labels = c("Never-smoker","Smoker"))+ ## change the x-axis label ## "All",
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size= 8, face = "bold"),
        axis.text.y = element_text( size= 8, face = "bold"),
        axis.title.y = element_text(size = 8, face= "bold"),
        axis.title.x = element_text(size = 8, face = "bold"),
        legend.position = "none") +
  facet_wrap(~variant_B_s, nrow=1)+
  theme(strip.text.x = element_text(angle = 45, hjust =0.5, vjust = 0.5,size =8, face = "bold")) 
  
### for KRAS_G12D ####
epi_KRAS_G12D <- epi_data["KRAS_G12D", on = "variant_A_type"]

## use the gene_order_KRAS_G12C
epi_KRAS_G12D[,variant_B_s := factor(variant_B_s, levels = gene_order_KRAS_G12C, ordered = T )] 
epi_KRAS_G12D[,data_type := factor(data_type, levels = c("nonSmo","Smo"), ordered = T )] 

## edit CI bar which has NA values,Jeff think that CI should be zero if it is NA
epi_KRAS_G12D$ci_low_95_ces_B_on_A_forplot <- ifelse(is.na(epi_KRAS_G12D$ci_low_95_ces_B_on_A),0, epi_KRAS_G12D$ci_low_95_ces_B_on_A)
epi_KRAS_G12D$errorBarColor <- ifelse(epi_KRAS_G12D$nAB %in% c("0","1"), "Gray","Black")
epi_KRAS_G12D$data_type2 <-  as.vector(epi_KRAS_G12D$data_type)
epi_KRAS_G12D$data_type2[epi_KRAS_G12D$nAB == "0"] <- "Gray"


p_G12D <- ggplot(data= epi_KRAS_G12D, aes(x= data_type , y = ces_B_on_A))+ 
  geom_errorbar(aes(ymin = ci_low_95_ces_B_on_A_forplot, ymax = ci_high_95_ces_B_on_A, color = errorBarColor), width = 0,  na.rm = T) + 
  geom_point(aes(fill = data_type2, color = errorBarColor),shape = 21, size =2) + 
  scale_fill_manual(values =c("Smo" = "red","nonSmo" = "blue","Gray" = "gray")) + 
  scale_color_manual(values = c("Gray" = "gray","Black" = "black"))+
  labs(x = "Tumor suppressor gene")+
  ylab(bquote(bold(atop("Cancer effect",
                        "(KRAS G12D background)"))))+ 
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = F))+
  scale_x_discrete(labels = c("Never-smoker","Smoker"))+ 
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size= 8, face = "bold"),
        axis.text.y = element_text( size= 8, face = "bold"),
        axis.title.y = element_text(size = 8, face= "bold"),
        axis.title.x = element_text(size = 8, face = "bold"),
        legend.position = "none") +
  facet_wrap(~variant_B_s, nrow=1)+
  theme(strip.text.x = element_text(angle = 45, hjust =0.5, vjust = 0.5,size =8, face = "bold")) #+

### for BRAF_V600E ####
epi_BRAF_V600E <- epi_data["BRAF_V600E", on = "variant_A_type"]

## use the gene_order_KRAS_G12C
epi_BRAF_V600E[,variant_B_s := factor(variant_B_s, levels = gene_order_KRAS_G12C, ordered = T )] 
epi_BRAF_V600E[,data_type := factor(data_type, levels = c("nonSmo","Smo"), ordered = T )] 

## edit CI bar which has NA values,Jeff think that CI should be zero if it is NA
epi_BRAF_V600E$ci_low_95_ces_B_on_A_forplot <- ifelse(is.na(epi_BRAF_V600E$ci_low_95_ces_B_on_A),0, epi_BRAF_V600E$ci_low_95_ces_B_on_A)
epi_BRAF_V600E$errorBarColor <- ifelse(epi_BRAF_V600E$nAB %in% c("0","1"), "Gray","Black")
epi_BRAF_V600E$data_type2 <-  as.vector(epi_BRAF_V600E$data_type)
epi_BRAF_V600E$data_type2[epi_BRAF_V600E$nAB == "0"] <- "Gray"


p_V600E <- ggplot(data= epi_BRAF_V600E, aes(x= data_type , y = ces_B_on_A))+ 
  geom_errorbar(aes(ymin = ci_low_95_ces_B_on_A_forplot, ymax = ci_high_95_ces_B_on_A, color = errorBarColor), width = 0,  na.rm = T) + 
  geom_point(aes(fill = data_type2, color = errorBarColor),shape = 21, size =2) + 
  scale_fill_manual(values =c("Smo" = "red","nonSmo" = "blue","Gray" = "gray")) + 
  scale_color_manual(values = c("Gray" = "gray","Black" = "black"))+
  labs(x = "Tumor suppressor gene")+
  ylab(bquote(bold(atop("Cancer effect",
                        "(BRAF V600E background)"))))+ 
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = F))+
  scale_x_discrete(labels = c("Never-smoker","Smoker"))+ 
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size= 8, face = "bold"),
        axis.text.y = element_text( size= 8, face = "bold"),
        axis.title.y = element_text(size = 8, face= "bold"),
        axis.title.x = element_text(size = 8, face = "bold"),
        legend.position = "none") +
  facet_wrap(~variant_B_s, nrow=1)+
  theme(strip.text.x = element_text(angle = 45, hjust =0.5, vjust = 0.5,size =8, face = "bold")) #+

### for EGFR_L858R ####
epi_EGFR_L858R <- epi_data["EGFR_L858R", on = "variant_A_type"]

## use the gene_order_KRAS_G12C
epi_EGFR_L858R[,variant_B_s := factor(variant_B_s, levels = gene_order_KRAS_G12C, ordered = T )] 
epi_EGFR_L858R[,data_type := factor(data_type, levels = c("nonSmo","Smo"), ordered = T )] 

## edit CI bar which has NA values, Jeff think that CI should be zero if it is NA
epi_EGFR_L858R$ci_low_95_ces_B_on_A_forplot <- ifelse(is.na(epi_EGFR_L858R$ci_low_95_ces_B_on_A),0, epi_EGFR_L858R$ci_low_95_ces_B_on_A) ## epi_EGFR_L858R$ces_B_on_A
epi_EGFR_L858R$errorBarColor <- ifelse(epi_EGFR_L858R$nAB %in% c("0","1"), "Gray","Black")
epi_EGFR_L858R$data_type2 <-  as.vector(epi_EGFR_L858R$data_type)
epi_EGFR_L858R$data_type2[epi_EGFR_L858R$nAB == "0"] <- "Gray"


p_L858R <- ggplot(data= epi_EGFR_L858R, aes(x= data_type , y = ces_B_on_A))+ 
  geom_errorbar(aes(ymin = ci_low_95_ces_B_on_A_forplot, ymax = ci_high_95_ces_B_on_A, color = errorBarColor), width = 0,  na.rm = T) + 
  geom_point(aes(fill = data_type2, color = errorBarColor),shape = 21, size =2) + 
  scale_fill_manual(values =c("Smo" = "red","nonSmo" = "blue","Gray" = "gray")) + 
  scale_color_manual(values = c("Gray" = "gray","Black" = "black"))+
  labs(x = "Tumor suppressor gene")+
  ylab(bquote(bold(atop("Cancer effect",
                        "(EGFR L858R background)"))))+ 
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = F))+
  scale_x_discrete(labels = c("Never-smoker","Smoker"))+ #
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size= 8, face = "bold"),
        axis.text.y = element_text( size= 8, face = "bold"),
        axis.title.y = element_text(size = 8, face= "bold"),
        axis.title.x = element_text(size = 8, face = "bold"),
        legend.position = "none") +
  facet_wrap(~variant_B_s, nrow=1)+
  theme(strip.text.x = element_text(angle = 45, hjust =0.5, vjust = 0.5,size =8, face = "bold"))




### plot CES_B_on_A of TSG under a specific variant background in one figure; only NS, S; sort by NS; NA is 0    ####
jpeg(file = "ces_B_on_A_variant_TSG_NAis0.jpeg", width = 13, height = 13, units = 'in', res = 300) 
p_G12C + p_G12D + p_V600E + p_L858R + plot_annotation(tag_levels = 'A') + plot_layout(nrow = 4, guides = "collect") 
dev.off()


  
## plot significant; use bar instead of point; combine two together ####
epi_G12C_L858R <- rbind(
  epi_EGFR_L858R[variant_B %in% c("TP53") ],
  epi_KRAS_G12C[variant_B %in% c("NF1") ])
variant_label_G12C_L858R <- epi_G12C_L858R[, .(max_y = max(ci_high_95_ces_B_on_A)), by = variant_B]
p_G12C_L858R_signifi_bar <- ggplot(data= epi_G12C_L858R, aes(x= variant_B , y = ces_B_on_A, fill = data_type))+ 
  ylab("Cancer effect size")+ 
  geom_bar(stat = "identity", color = "black",position = position_dodge(width = 0.9)) + ## color = errorBarColor,shape = 21, size =6
  geom_errorbar(aes(ymin = ci_low_95_ces_B_on_A_forplot, ymax = ci_high_95_ces_B_on_A), width = 0.3,  na.rm = T, linewidth = 0.6, position = position_dodge(width = 0.9)) + ## , color = errorBarColor
  # Add variant_name labels above the bars
  geom_text(data = variant_label_G12C_L858R,
            aes(x = variant_B, y = max(variant_label$max_y) * 3, label = variant_B),
            inherit.aes = F,
            vjust = 0,
            size = 5,
            fontface = "bold") +
  scale_x_discrete(labels = rep(c("Never-smoker", "Smoker"), length.out = 2))+
  scale_fill_manual(values =c("Smo" = "#D2B48C","nonSmo" = "pink")) + ## ,"Gray" = "gray" ## "All" = "black", ## brown looks like red,so use light brown hex code #D2B48C instead 
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), 
                     breaks = c(0.01, 1, 100, 10000,1000000),
                     labels = scales::label_comma(accuracy = 1))+
  theme_classic()+ ##theme_classic()
  theme(legend.direction = "horizontal",
        legend.position = "bottom",
        legend.box = "vertical",
        legend.text = element_text(size = 18),
        legend.title= element_text(size = 18),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size =18,face = "bold"),
        axis.title.y = element_text( size = 18,face= "bold"),
        axis.title.x = element_blank()) +
  guides(fill = F)

svg(file = "CES_select_VariantGeneLevel.svg", height =10.5, width = 5) 
p_G12C_L858R_signifi_bar
dev.off()









