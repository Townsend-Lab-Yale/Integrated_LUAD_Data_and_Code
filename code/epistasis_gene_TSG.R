library(ggplot2)
library(ggpubr)
require(patchwork)
## epistasis analysis between a oncogene and a tsg gene; variants selection: select; noncoding; recurrent #####
get_gene_epistasis_cesa<- function(cesa_object,variants_type, gene_onco, gene_tsg, specific_samples){ 
  if (variants_type %in% c("select_")){
    ## Way1: using select_variants() to include amino-acid-change mutation and noncoding SNVs
    # gene_onco <- c("EGFR")
    # gene_tsg <- "TP53"
    # cesa_object <- cesa
    variants_geneOnco <- cancereffectsizeR::select_variants(cesa = cesa_object, genes = gene_onco) 
    variants_geneTsg <- cancereffectsizeR::select_variants(cesa = cesa_object, genes = gene_tsg) 
    variants_target <- rbind(variants_geneOnco, variants_geneTsg)
    desired_order <- c(gene_onco, gene_tsg)
    variants_target <- variants_target[order(match(gene, desired_order))]
    
  }else if(variants_type %in% c("recurrent_")){
    ## Way2: setting variants parameter as "recurrent" to include all recurrent variants (maf_prevalence >1)
    variants_target <- "recurrent"
    
  }else if(variants_type %in% c("nonsilent_")){
    ## Way3: setting variants parameter as "nonsilent" to include all nonsynonymous coding variants and any variants in essential splice sites 
    variants_target <- "nonsilent"
  }
  genes_target = c(gene_onco, gene_tsg)
  run_name <- paste0(variants_type,gene_onco,"_vs_", gene_tsg)
  cesa_object <- ces_gene_epistasis(cesa = cesa_object, genes = genes_target, variants = variants_target, run_name = run_name, samples = specific_samples)
  epi <- cesa_object$epistasis[[run_name]]
  set(epi, j= "run_type", value = run_name) 
  result <- list(epi, cesa_object, comp_target, variants_target)
  return(result)
}


#Pre-allocate epi_final with column names from the first iteration
KRAS_vs_PTEN <- get_gene_epistasis_cesa(cesa_object = cesa, variants_type= "recurrent_", gene_onco = "KRAS", gene_tsg = "PTEN", specific_samples = cesa$samples$Unique_Patient_Identifier)
epi_v3_initial <- KRAS_vs_PTEN[[1]]
epi_v3_initial[1,] <- "NA"
epi_v3_final <-  data.table(matrix(NA, nrow = 0, ncol = length(epi_v3_initial))) ## create an empty data table with 0 rows and length(epi_v3_initial) columns
setnames(epi_v3_final, names(epi_v3_initial)) ## then set the column names
epi_v3_final_smoking <- epi_v3_final
epi_v3_final_nonsmoking <- epi_v3_final

oncogene <- c("KRAS","BRAF","EGFR")
for(i in 1:length(oncogene)){
  for(j in 1:length(tsg)){ 
    ## For all samples
    print("For cesa, oncogene[i], tsg[j], variants_type: select ")
    epi_v3_select <- get_gene_epistasis_cesa(cesa_object = cesa, variants_type = "select_", gene_onco = oncogene[i], gene_tsg = tsg[j], specific_samples = cesa$samples$Unique_Patient_Identifier)[[1]] ## select
    epi_v3_final <- rbindlist(list(epi_v3_final, epi_v3_select), fill = F)
    if(nrow(cancereffectsizeR::select_variants(cesa, gene = tsg[j])[maf_prevalence >1]) > 0){
      epi_v3_recurrent <- get_gene_epistasis_cesa(cesa_object = cesa, variants_type = "recurrent_", gene_onco = oncogene[i], gene_tsg = tsg[j], specific_samples = cesa$samples$Unique_Patient_Identifier)[[1]] ## recurrent_
    }else{
      epi_v3_recurrent <- epi_v3_initial
      epi_v3_recurrent[,c("variant_A","variant_B", "run_type") := list(oncogene[i], tsg[j], paste0("recurrent_",oncogene[i],"_vs_", tsg[j]))]
    }
    epi_v3_final <- rbindlist(list(epi_v3_final, epi_v3_recurrent), fill = F)
    epi_v3_nonsilent <- get_gene_epistasis_cesa(cesa_object = cesa, variants_type = "nonsilent_", gene_onco = oncogene[i], gene_tsg = tsg[j], specific_samples = cesa$samples$Unique_Patient_Identifier)[[1]] ## nonsilent_
    epi_v3_final <- rbindlist(list(epi_v3_final, epi_v3_nonsilent), fill = F)
    
    ## For smokers
    epi_v3_select_smoking <- get_gene_epistasis_cesa(cesa_object = cesa_smoking_w_panel, variants_type = "select_", gene_onco = oncogene[i], gene_tsg = tsg[j], 
                                                     specific_samples = cesa_smoking_w_panel$samples[Unique_Patient_Identifier %in% c(smoking_samples, panel_smoking_samples), Unique_Patient_Identifier])[[1]] ## select
    epi_v3_final_smoking <- rbindlist(list(epi_v3_final_smoking, epi_v3_select_smoking), fill = F)
    if(nrow(cancereffectsizeR::select_variants(cesa_smoking_w_panel, gene = tsg[j])[maf_prevalence >1]) > 0){
      epi_v3_recurrent_smoking <- get_gene_epistasis_cesa(cesa_object = cesa_smoking_w_panel, variants_type = "recurrent_", gene_onco = oncogene[i], gene_tsg = tsg[j],
                                                          specific_samples = cesa_smoking_w_panel$samples[Unique_Patient_Identifier %in% c(smoking_samples, panel_smoking_samples), Unique_Patient_Identifier])[[1]] ## recurrent_
    }else{
      epi_v3_recurrent_smoking <- epi_v3_initial
      epi_v3_recurrent_smoking[,c("variant_A","variant_B", "run_type") := list(oncogene[i], tsg[j], paste0("recurrent_",oncogene[i],"_vs_", tsg[j]))]
    }
    epi_v3_final_smoking <- rbindlist(list(epi_v3_final_smoking, epi_v3_recurrent_smoking), fill = F)
    epi_v3_nonsilent_smoking <- get_gene_epistasis_cesa(cesa_object = cesa_smoking_w_panel, variants_type = "nonsilent_", gene_onco = oncogene[i], gene_tsg = tsg[j], 
                                                        specific_samples = cesa_smoking_w_panel$samples[Unique_Patient_Identifier %in% c(smoking_samples, panel_smoking_samples), Unique_Patient_Identifier])[[1]] ## nonsilent_
    epi_v3_final_smoking <- rbindlist(list(epi_v3_final_smoking, epi_v3_nonsilent_smoking), fill = F)
    
    ## For never-smoker
    epi_v3_select_nonsmoking <- get_gene_epistasis_cesa(cesa_object = cesa_nonsmoking_w_panel, variants_type = "select_", gene_onco = oncogene[i], gene_tsg = tsg[j], 
                                                     specific_samples = cesa_nonsmoking_w_panel$samples[Unique_Patient_Identifier %in% c(nonsmoking_samples, panel_nonsmoking_samples), Unique_Patient_Identifier])[[1]] ## select
    epi_v3_final_nonsmoking <- rbindlist(list(epi_v3_final_nonsmoking, epi_v3_select_nonsmoking), fill = F)
    if(nrow(cancereffectsizeR::select_variants(cesa_nonsmoking_w_panel, gene = tsg[j])[maf_prevalence >1]) > 0){
      epi_v3_recurrent_nonsmoking <- get_gene_epistasis_cesa(cesa_object = cesa_nonsmoking_w_panel, variants_type = "recurrent_", gene_onco = oncogene[i], gene_tsg = tsg[j],
                                                             specific_samples = cesa_nonsmoking_w_panel$samples[Unique_Patient_Identifier %in% c(nonsmoking_samples, panel_nonsmoking_samples), Unique_Patient_Identifier])[[1]] ## recurrent_
    }else{
      epi_v3_recurrent_nonsmoking <- epi_v3_initial
      epi_v3_recurrent_nonsmoking[,c("variant_A","variant_B", "run_type") := list(oncogene[i], tsg[j], paste0("recurrent_",oncogene[i],"_vs_", tsg[j]))]
    }
    epi_v3_final_nonsmoking <- rbindlist(list(epi_v3_final_nonsmoking, epi_v3_recurrent_nonsmoking), fill = F)
    epi_v3_nonsilent_nonsmoking <- get_gene_epistasis_cesa(cesa_object = cesa_nonsmoking_w_panel, variants_type = "nonsilent_", gene_onco = oncogene[i], gene_tsg = tsg[j], 
                                                        specific_samples = cesa_nonsmoking_w_panel$samples[Unique_Patient_Identifier %in% c(nonsmoking_samples, panel_nonsmoking_samples), Unique_Patient_Identifier])[[1]] ## nonsilent_
    epi_v3_final_nonsmoking <- rbindlist(list(epi_v3_final_nonsmoking, epi_v3_nonsilent_nonsmoking), fill = F)
  }
}


epi_v3_data <- rbind(epi_v3_final[, data_type := "All"],
                     epi_v3_final_smoking[, data_type := "Smo"],epi_v3_final_nonsmoking[, data_type := "nonSmo"])
# View(epi_v3_data)
epi_v3_data_p1 <- epi_v3_data[c("KRAS", "BRAF", "EGFR"), on = "variant_A"]
epi_v3_data_p2 <- epi_v3_data[c("KRAS", "BRAF", "EGFR"), on = "variant_B"] ## KRAS/BRAF/EGFR at variantB position
# View(epi_v3_data_p2)
## change the order to let KRAS/BRAF/EGFR at variantA position
epi_v3_data_p2_trans <- data.table(matrix(NA, nrow = nrow(epi_v3_data_p2) , ncol = length(epi_v3_data_p2)))
setnames(epi_v3_data_p2_trans, names(epi_v3_data_p2))
View(epi_v3_data_p2_trans)
colnames(epi_v3_data_p2_trans)
epi_v3_data_p2_trans$variant_A <- epi_v3_data_p2$variant_B
epi_v3_data_p2_trans$variant_B <- epi_v3_data_p2$variant_A
epi_v3_data_p2_trans$ces_A0 <- epi_v3_data_p2$ces_B0
epi_v3_data_p2_trans$ces_B0 <- epi_v3_data_p2$ces_A0
epi_v3_data_p2_trans$ces_A_on_B <- epi_v3_data_p2$ces_B_on_A
epi_v3_data_p2_trans$ces_B_on_A <- epi_v3_data_p2$ces_A_on_B
epi_v3_data_p2_trans$p_A_change <- epi_v3_data_p2$p_B_change
epi_v3_data_p2_trans$p_B_change <- epi_v3_data_p2$p_A_change
epi_v3_data_p2_trans$p_epistasis <- epi_v3_data_p2$p_epistasis
epi_v3_data_p2_trans$expected_nAB_epistasis <- epi_v3_data_p2$expected_nAB_epistasis
epi_v3_data_p2_trans$expected_nAB_null <- epi_v3_data_p2$expected_nAB_null
epi_v3_data_p2_trans$AB_epistatic_ratio <- epi_v3_data_p2$AB_epistatic_ratio
epi_v3_data_p2_trans$nA0 <- epi_v3_data_p2$nB0
epi_v3_data_p2_trans$nB0 <- epi_v3_data_p2$nA0
epi_v3_data_p2_trans$nAB <- epi_v3_data_p2$nAB
epi_v3_data_p2_trans$n00 <- epi_v3_data_p2$n00
epi_v3_data_p2_trans$n_total <- epi_v3_data_p2$n_total
epi_v3_data_p2_trans$ci_low_95_ces_A0 <- epi_v3_data_p2$ci_low_95_ces_B0
epi_v3_data_p2_trans$ci_high_95_ces_A0 <- epi_v3_data_p2$ci_high_95_ces_B0
epi_v3_data_p2_trans$ci_low_95_ces_B0 <- epi_v3_data_p2$ci_low_95_ces_A0
epi_v3_data_p2_trans$ci_high_95_ces_B0 <- epi_v3_data_p2$ci_high_95_ces_A0
epi_v3_data_p2_trans$ci_low_95_ces_A_on_B <- epi_v3_data_p2$ci_low_95_ces_B_on_A
epi_v3_data_p2_trans$ci_high_95_ces_A_on_B <- epi_v3_data_p2$ci_high_95_ces_B_on_A
epi_v3_data_p2_trans$ci_low_95_ces_B_on_A <- epi_v3_data_p2$ci_low_95_ces_A_on_B
epi_v3_data_p2_trans$ci_high_95_ces_B_on_A <- epi_v3_data_p2$ci_high_95_ces_A_on_B
epi_v3_data_p2_trans$run_type <- epi_v3_data_p2$run_type
epi_v3_data_p2_trans$data_type <- epi_v3_data_p2$data_type

## add variants_type column 
epi_v3_data2 <- rbind(epi_v3_data_p1, epi_v3_data_p2_trans)
epi_v3_data2[,"variants_type":= ifelse(stringr::str_detect(run_type,"select"),"select",ifelse(stringr::str_detect(run_type,"recurrent"),"recurrent","nonsilent"))]
table(epi_v3_data2$variants_type, exclude = NULL)
## add variant_B_s column, to shorten CDKN2A.p16INK4a and CDKN2A.p14arf to p16INK4a and p14arf
epi_v3_data2[,variant_B_s := variant_B]
epi_v3_data2[stringr::str_detect(variant_B_s,"\\."),variant_B_s := stringr::str_split(variant_B_s,"\\.", simplify = T)[,2]]
View(epi_v3_data2)
str(epi_v3_data2)
head(epi_v3_data2[ces_B_on_A == 0.001000])
## edit ci_low_95_ces_B_on_A NA to ces_B_on_A to show CI bar
# nrow(epi_v3_data2[is.na(ci_low_95_ces_B_on_A) & !is.na(ces_B_on_A)]) ## 431 rows of ci_low_95_ces_B_on_A are NAs
# nrow(epi_v3_data2[is.na(ces_B_on_A) & is.na(n00)]) ## 9 rows are those no samples covering one variants. 
# nrow(epi_v3_data2[is.na(ces_B_on_A) & !is.na(n00)]) ## 36 rows are those gene pairs had no eligible variants in one or both genes of jointly-covering samples. Epistatic selection intensities are all NAs for these pairs. 

epi_v3_data2[, "ci_low_95_ces_B_on_A_new" := ifelse(is.na(ci_low_95_ces_B_on_A) & !is.na(ces_B_on_A) ,ces_B_on_A, ci_low_95_ces_B_on_A)]

## Visualizing the ces_B_on_A #######
## here we use "select" variants (all amino-acid-change mutation and noncoding SNVs are included, so the sample sizes are more powerful than "nonsilent" and "recurrent")
variants_type_target <- "select" ## "nonsilent"; "recurrent"
### for KRAS ####
epi_v3_KRAS_all <- epi_v3_data2[.("KRAS","All",variants_type_target), on=.(variant_A, data_type, variants_type)]
epi_v3_KRAS_S <- epi_v3_data2[.("KRAS","Smo",variants_type_target), on=.(variant_A, data_type, variants_type)]
epi_v3_KRAS_nonS <- epi_v3_data2[.("KRAS","nonSmo",variants_type_target), on=.(variant_A, data_type, variants_type)]
epi_v3_KRAS <- epi_v3_data2[.("KRAS",variants_type_target), on=.(variant_A, variants_type)]
epi_v3_KRAS_NA <- epi_v3_KRAS[is.na(ces_B_on_A)]
nrow(epi_v3_KRAS_NA) ## there are 15个 gene pair have NA ces_B_on_A
epi_v3_KRAS <- epi_v3_KRAS[!is.na(ces_B_on_A)]
## sorting gene by ces_B_on_A of TSG under KRSA_G12C background in smokers, then neversmokers, then all samples 
variant_order_KRAS <- epi_v3_KRAS_S[ces_B_on_A > 0.001][order(-ces_B_on_A)]$variant_B_s
variant_lf_KRAS <- epi_v3_KRAS_S[(ces_B_on_A == 0.001000) | is.na(ces_B_on_A)]$variant_B_s
if(length(variant_lf_KRAS) >=1){
  variant_lf_order_KRAS <- epi_v3_KRAS_nonS[variant_lf_KRAS, on = "variant_B_s"][ces_B_on_A > 0.001000][order(-ces_B_on_A)]$variant_B_s
  variant_lf_2_KRAS <- variant_lf_KRAS[!(variant_lf_KRAS %in% variant_lf_order_KRAS)]
  variant_lf_2_order_KRAS <- epi_v3_KRAS_all[variant_lf_2_KRAS, on = "variant_B_s"][ces_B_on_A > 0.001000][order(-ces_B_on_A)]$variant_B_s
  variant_lf_3_KRAS <- variant_lf_KRAS[!(variant_lf_KRAS %in% c(variant_lf_order_KRAS, variant_lf_2_order_KRAS))]
  variant_lf_3_KRAS <- epi_v3_KRAS_all[variant_lf_2_KRAS, on = "variant_B_s"][!(ces_B_on_A > 0.001000)]$variant_B_s
  if (length(variant_lf_3_KRAS) >= 1){
    variant_lf_3_order_KRAS <- c("NF2","MSH2")## manually order
    gene_order_KRAS <- c(variant_order_KRAS,variant_lf_order_KRAS,variant_lf_2_order_KRAS,variant_lf_3_order_KRAS)
  }else{
    gene_order_KRAS <- c(variant_order_KRAS,variant_lf_order_KRAS,variant_lf_2_order_KRAS)
  }
}else{
  gene_order_KRAS <- variant_order_KRAS
}
length(gene_order_KRAS) ## 28
epi_v3_KRAS[,variant_B_s := factor(variant_B_s, levels = gene_order_KRAS, ordered = T )] 
epi_v3_KRAS[,data_type := factor(data_type, levels = c("All","Smo","nonSmo"), ordered = T )] 


epi_v3_KRAS$errorBarColor <- ifelse(epi_v3_KRAS$nAB %in% c("0","1"), "Gray","Black")
epi_v3_KRAS$data_type2 <-  as.vector(epi_v3_KRAS$data_type)
epi_v3_KRAS$data_type2[epi_v3_KRAS$nAB == "0"] <- "Gray"

# ## add epi_data_t_test ##
# epi_KRAS_t_test <- get_epi_data_t_test(epi_v3_KRAS,tsg)
# epi_KRAS_t_test$p_value <- ifelse(stringr::str_detect(epi_KRAS_t_test$t_test_p_mia_v2,"nAB"), NA, as.numeric(epi_KRAS_t_test$t_test_p_mia_v2)) 
# epi_KRAS_t_test[,gene_s := gene]
# epi_KRAS_t_test[stringr::str_detect(gene_s,"\\."),gene_s := stringr::str_split(gene_s,"\\.", simplify = T)[,2]]
# epi_KRAS_t_test <- epi_KRAS_t_test[order(match(epi_KRAS_t_test$gene_s, gene_order_KRAS))]
# epi_KRAS_t_test$p_value_anno <- ifelse(is.na(epi_KRAS_t_test$p_value), "", ifelse(epi_KRAS_t_test$p_value < 0.001, "***", ifelse(epi_KRAS_t_test$p_value < 0.01, "**", ifelse(epi_KRAS_t_test$p_value < 0.05, "*", "ns"))))
# setnames(epi_KRAS_t_test, old = "gene_s", new = "variant_B_s")
# epi_KRAS_t_test[,variant_B_s := factor(variant_B_s, levels = gene_order_KRAS, ordered = T)]

p_KRAS <- ggplot(data= epi_v3_KRAS, aes(x= data_type , y = ces_B_on_A))+ 
   geom_errorbar(aes(ymin = ci_low_95_ces_B_on_A_new, ymax = ci_high_95_ces_B_on_A, color = errorBarColor), width = 0,  na.rm = T) + 
  geom_point(aes(fill = data_type2, color = errorBarColor),shape = 21, size =2) + 
  scale_fill_manual(values =c("All" = "black", "Smo" = "red","nonSmo" = "blue","Gray" = "gray")) +
  scale_color_manual(values = c("Gray" = "gray","Black" = "black"))+
  labs(x = "Tumor suppressor gene")+
  ylab(bquote(bold(atop("Cancer effect",
                        "(KRAS background)"))))+ 
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = F))+
  scale_x_discrete(labels = c("All","Smoker","Never-smoker"))+ ## change the x-axis label
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size= 6),
        axis.title.y = element_text(size = 8, face= "bold"),
        axis.title.x = element_text(size = 8, face = "bold"),
        legend.position = "none") +
  facet_wrap(~variant_B_s, nrow=1)+
  theme(strip.text.x = element_text(angle = 45, hjust =0.5, vjust = 0.5,size =7, face = "bold")) #+
  # geom_text(
  #   data = epi_KRAS_t_test,
  #   aes(x=2.5,y=max(epi_v3_KRAS$ci_high_95_ces_B_on_A)+1, label = p_value_anno),
  #   hjust = 0.5, vjust = 0.5 ) 

### for BRAF ####
# epi_v3_BRAF_all <- epi_v3_data2[.("BRAF","All",variants_type_target), on=.(variant_A, data_type, variants_type)]
# epi_v3_BRAF_S <- epi_v3_data2[.("BRAF","Smo",variants_type_target), on=.(variant_A, data_type, variants_type)]
# epi_v3_BRAF_nonS <- epi_v3_data2[.("BRAF","nonSmo",variants_type_target), on=.(variant_A, data_type, variants_type)]
epi_v3_BRAF <- epi_v3_data2[.("BRAF",variants_type_target), on=.(variant_A, variants_type)]
epi_v3_BRAF_NA <- epi_v3_BRAF[is.na(ces_B_on_A)]
nrow(epi_v3_BRAF_NA) ## there are 15个 gene pair have NA ces_B_on_A
epi_v3_BRAF <- epi_v3_BRAF[!is.na(ces_B_on_A)]
## use the gene_order_KRAS
epi_v3_BRAF[,variant_B_s := factor(variant_B_s, levels = gene_order_KRAS, ordered = T )] 
epi_v3_BRAF[,data_type := factor(data_type, levels = c("All","Smo","nonSmo"), ordered = T )] 

epi_v3_BRAF$errorBarColor <- ifelse(epi_v3_BRAF$nAB %in% c("0","1"), "Gray","Black")
epi_v3_BRAF$data_type2 <-  as.vector(epi_v3_BRAF$data_type)
epi_v3_BRAF$data_type2[epi_v3_BRAF$nAB == "0"] <- "Gray"

# ## add epi_data_t_test ##
# epi_BRAF_t_test <- get_epi_data_t_test(epi_v3_BRAF,tsg)
# epi_BRAF_t_test$p_value <- ifelse(stringr::str_detect(epi_BRAF_t_test$t_test_p_mia_v2,"nAB"), NA, as.numeric(epi_BRAF_t_test$t_test_p_mia_v2)) 
# epi_BRAF_t_test[,gene_s := gene]
# epi_BRAF_t_test[stringr::str_detect(gene_s,"\\."),gene_s := stringr::str_split(gene_s,"\\.", simplify = T)[,2]]
# epi_BRAF_t_test <- epi_BRAF_t_test[order(match(epi_BRAF_t_test$gene_s, gene_order_BRAF))]
# epi_BRAF_t_test$p_value_anno <- ifelse(is.na(epi_BRAF_t_test$p_value), "", ifelse(epi_BRAF_t_test$p_value < 0.001, "***", ifelse(epi_BRAF_t_test$p_value < 0.01, "**", ifelse(epi_BRAF_t_test$p_value < 0.05, "*", "ns"))))
# setnames(epi_BRAF_t_test, old = "gene_s", new = "variant_B_s")
# epi_BRAF_t_test[,variant_B_s := factor(variant_B_s, levels = gene_order_BRAF, ordered = T)]

p_BRAF <- ggplot(data= epi_v3_BRAF, aes(x= data_type , y = ces_B_on_A))+ 
  geom_errorbar(aes(ymin = ci_low_95_ces_B_on_A_new, ymax = ci_high_95_ces_B_on_A, color = errorBarColor), width = 0,  na.rm = T) + 
  geom_point(aes(fill = data_type2, color = errorBarColor),shape = 21, size =2) + 
  scale_fill_manual(values =c("All" = "black", "Smo" = "red","nonSmo" = "blue","Gray" = "gray")) +
  scale_color_manual(values = c("Gray" = "gray","Black" = "black"))+
  labs(x = "Tumor suppressor gene")+
  ylab(bquote(bold(atop("Cancer effect",
                        "(BRAF background)"))))+ 
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = F))+
  scale_x_discrete(labels = c("All","Smoker","Never-smoker"))+ ## change the x-axis label
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size= 6),
        axis.title.y = element_text(size = 8, face= "bold"),
        axis.title.x = element_text(size = 8, face = "bold"),
        legend.position = "none") +
  facet_wrap(~variant_B_s, nrow=1)+
  theme(strip.text.x = element_text(angle = 45, hjust =0.5, vjust = 0.5,size =7, face = "bold")) #+
# geom_text(
#   data = epi_BRAF_t_test,
#   aes(x=2.5,y=max(epi_v3_BRAF$ci_high_95_ces_B_on_A)+1, label = p_value_anno),
#   hjust = 0.5, vjust = 0.5 ) 

### for EGFR ####
# epi_v3_EGFR_all <- epi_v3_data2[.("EGFR","All",variants_type_target), on=.(variant_A, data_type, variants_type)]
# epi_v3_EGFR_S <- epi_v3_data2[.("EGFR","Smo",variants_type_target), on=.(variant_A, data_type, variants_type)]
# epi_v3_EGFR_nonS <- epi_v3_data2[.("EGFR","nonSmo",variants_type_target), on=.(variant_A, data_type, variants_type)]
epi_v3_EGFR <- epi_v3_data2[.("EGFR",variants_type_target), on=.(variant_A, variants_type)]
epi_v3_EGFR_NA <- epi_v3_EGFR[is.na(ces_B_on_A)]
nrow(epi_v3_EGFR_NA) ## there are 15个 gene pair have NA ces_B_on_A
epi_v3_EGFR <- epi_v3_EGFR[!is.na(ces_B_on_A)]
## use the gene_order_KRAS
epi_v3_EGFR[,variant_B_s := factor(variant_B_s, levels = gene_order_KRAS, ordered = T )] 
epi_v3_EGFR[,data_type := factor(data_type, levels = c("All","Smo","nonSmo"), ordered = T )] 

epi_v3_EGFR$errorBarColor <- ifelse(epi_v3_EGFR$nAB %in% c("0","1"), "Gray","Black")
epi_v3_EGFR$data_type2 <-  as.vector(epi_v3_EGFR$data_type)
epi_v3_EGFR$data_type2[epi_v3_EGFR$nAB == "0"] <- "Gray"

# ## add epi_data_t_test ##
# epi_EGFR_t_test <- get_epi_data_t_test(epi_v3_EGFR,tsg)
# epi_EGFR_t_test$p_value <- ifelse(stringr::str_detect(epi_EGFR_t_test$t_test_p_mia_v2,"nAB"), NA, as.numeric(epi_EGFR_t_test$t_test_p_mia_v2)) 
# epi_EGFR_t_test[,gene_s := gene]
# epi_EGFR_t_test[stringr::str_detect(gene_s,"\\."),gene_s := stringr::str_split(gene_s,"\\.", simplify = T)[,2]]
# epi_EGFR_t_test <- epi_EGFR_t_test[order(match(epi_EGFR_t_test$gene_s, gene_order_EGFR))]
# epi_EGFR_t_test$p_value_anno <- ifelse(is.na(epi_EGFR_t_test$p_value), "", ifelse(epi_EGFR_t_test$p_value < 0.001, "***", ifelse(epi_EGFR_t_test$p_value < 0.01, "**", ifelse(epi_EGFR_t_test$p_value < 0.05, "*", "ns"))))
# setnames(epi_EGFR_t_test, old = "gene_s", new = "variant_B_s")
# epi_EGFR_t_test[,variant_B_s := factor(variant_B_s, levels = gene_order_EGFR, ordered = T)]

p_EGFR <- ggplot(data= epi_v3_EGFR, aes(x= data_type , y = ces_B_on_A))+ 
  geom_errorbar(aes(ymin = ci_low_95_ces_B_on_A_new, ymax = ci_high_95_ces_B_on_A, color = errorBarColor), width = 0,  na.rm = T) + 
  geom_point(aes(fill = data_type2, color = errorBarColor),shape = 21, size =2) + 
  scale_fill_manual(values =c("All" = "black", "Smo" = "red","nonSmo" = "blue","Gray" = "gray")) +
  scale_color_manual(values = c("Gray" = "gray","Black" = "black"))+
  labs(x = "Tumor suppressor gene")+
  ylab(bquote(bold(atop("Cancer effect",
                        "(EGFR background)"))))+ 
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = F))+
  scale_x_discrete(labels = c("All","Smoker","Never-smoker"))+ ## change the x-axis label
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size= 6),
        axis.title.y = element_text(size = 8, face= "bold"),
        axis.title.x = element_text(size = 8, face = "bold"),
        legend.position = "none") +
  facet_wrap(~variant_B_s, nrow=1)+
  theme(strip.text.x = element_text(angle = 45, hjust =0.5, vjust = 0.5,size =7, face = "bold")) #+
# geom_text(
#   data = epi_EGFR_t_test,
#   aes(x=2.5,y=max(epi_v3_EGFR$ci_high_95_ces_B_on_A)+1, label = p_value_anno),
#   hjust = 0.5, vjust = 0.5 ) 

### plot CES_B_on_A of TSG under a specific oncogene background in one figure  ####
jpeg(file = "CES_B_on_A_gene_TSG.jpeg", width = 13, height = 13, units = 'in', res = 300) ## 6inch = 1800 pixels. 1278 pixels, 1000
p_KRAS + p_BRAF + p_EGFR + plot_annotation(tag_levels = 'A') + plot_layout(nrow = 3, guides = "collect") 
dev.off()

### plot CES_B_on_A of TSG under a specific oncogene background; only show those points with significant difference between smokers and never-smokers ####
# ## use t_test p value as the standards 
# epi_EGFR_t_test_signifi <- epi_EGFR_t_test[p_value_anno %in% c("*","**","***")]
# gene_order_signifi <- epi_EGFR_t_test_signifi[order(smo_mean, decreasing = T), as.character(variant_B_s)]
# epi_EGFR_t_test_signifi[,variant_B_s := as.character(variant_B_s) ] ## unfactor the variant_B_s
# epi_EGFR_t_test_signifi[,variant_B_s := factor(variant_B_s, levels = gene_order_signifi, ordered = T)]
# epi_EGFR_t_test_signifi[,p_value_anno_nonOverlap:= "*"]



## using non-overlapping CI as the standards; only show those points with non-overlapping CIs ####
gene_tsg_signifi_gene <- c("KEAP1","KMT2D","p16INK4a","ARID2","ATRX","ATM","APC")
epi_v3_EGFR_signifi <-  epi_v3_EGFR[variant_B_s %in% variant_tsg_signifi_gene]
gene_tsg_signifi <- data.frame(gene_tsg_signifi_gene = gene_tsg_signifi_gene,
                                  p_value_anno_nonOverlap = c("*"))

p_EGFR_select_signifi <- ggplot(data= epi_v3_EGFR_signifi, aes(x= data_type , y = ces_B_on_A))+ 
  geom_errorbar(aes(ymin = ci_low_95_ces_B_on_A_new, ymax = ci_high_95_ces_B_on_A, color = errorBarColor), width = 0,  na.rm = T) + 
  geom_point(aes(fill = data_type2, color = errorBarColor),shape = 21, size =6) + 
  scale_fill_manual(values =c("All" = "black", "Smo" = "red","nonSmo" = "blue","Gray" = "gray")) +
  scale_color_manual(values = c("Gray" = "gray","Black" = "black"))+
  labs(x = "Tumor suppressor gene")+
  ylab(bquote(bold(atop("Cancer effect",
                        "(EGFR background)"))))+ 
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = F))+
  scale_x_discrete(labels = c("All","Smoker","Never-smoker"))+ ## change the x-axis label
  theme_minimal()+
  mia_theme_forPoster+
  facet_wrap(~variant_B_s, nrow=1)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size=16,face = "bold"),
        axis.text.y = element_text(size =16,face = "bold"),
        axis.title.y = element_text( size = 18,face= "bold"),
        axis.title.x = element_text( size = 18,face = "bold"),
        strip.text.x = element_text(angle = 45, hjust =0.5, vjust = 0.5,size =16, face = "bold")) + 
  geom_text(
    data = gene_tsg_signifi,
    aes(x=2.5,y=max(epi_v3_EGFR_signifi$ci_high_95_ces_B_on_A, na.rm = T)+1, label = p_value_anno_nonOverlap),
    hjust = 0.5, vjust = 0, size = 6 ) ## size = 4,

svg(file = "CES_GeneTSG_signifi_byNonOverlapping.svg", height =9.60, width = 9.10 ) ## width = 13, height =10.626
p_EGFR_select_signifi
dev.off()



## sample size forLRT ####
KRAS_samples_forLRT <- epi_v3_KRAS[!"All", on="data_type",.(variant_A, variant_B, nA0, nB0, nAB, n00,data_type)]
write.csv(KRAS_samples_forLRT, file = "KRAS_samples_forLRT.csv", row.names = F)

BRAF_samples_forLRT <- epi_v3_BRAF[!"All", on="data_type",.(variant_A, variant_B, nA0, nB0, nAB, n00,data_type)]
write.csv(BRAF_samples_forLRT, file = "BRAF_samples_forLRT.csv", row.names = F)

EGFR_samples_forLRT <- epi_v3_EGFR[!"All", on="data_type",.(variant_A, variant_B, nA0, nB0, nAB, n00,data_type)]
write.csv(EGFR_samples_forLRT, file = "EGFR_samples_forLRT.csv", row.names = F)

## CES value only in All group; sample size forLRT ####
geneLevel_forLRT<- rbind(epi_v3_KRAS["All", on="data_type",.(variant_A, variant_B, ces_B_on_A, nA0, nB0, nAB, n00,data_type)],
                         epi_v3_BRAF["All", on="data_type",.(variant_A, variant_B, ces_B_on_A, nA0, nB0, nAB, n00,data_type)],
                         epi_v3_EGFR["All", on="data_type",.(variant_A, variant_B, ces_B_on_A, nA0, nB0, nAB, n00,data_type)])
View(geneLevel_forLRT)
write.csv(geneLevel_forLRT, file = "sampleSize_forLRT_TSGinOncogene_allSamples.csv", row.names = F)


##  save cesa oject ####
save_cesa(cesa = cesa, "cesa_allSamples.rds")
save_cesa(cesa = cesa_smoking_w_panel, "cesa_smoking.rds")
save_cesa(cesa = cesa_nonsmoking_w_panel, "cesa_nonsmoking.rds")

## save smokers samples and never-smokers samples ####
samples_smokers <-  cesa_smoking_w_panel$samples[Unique_Patient_Identifier %in% c(smoking_samples, panel_smoking_samples), Unique_Patient_Identifier]
samples_neversmokers <-  cesa_nonsmoking_w_panel$samples[Unique_Patient_Identifier %in% c(nonsmoking_samples, panel_nonsmoking_samples), Unique_Patient_Identifier]
write.csv(samples_smokers, file = "Unique_Patient_Identifier_smokers_1082.csv")
write.csv(samples_neversmokers, file = "Unique_Patient_Identifier_nonsmokers_642.csv")




