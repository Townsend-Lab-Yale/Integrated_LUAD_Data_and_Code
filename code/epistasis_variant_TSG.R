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
KRAS_G12D_vs_PTEN <- get_epistasis_cesa(cesa_object = cesa, variant = oncogene_variants[1], gene_name = tsg[1], specific_samples = cesa$samples$Unique_Patient_Identifier)
epi_initial <- KRAS_G12D_vs_PTEN[[1]]
epi_final <-  data.table(matrix(NA, nrow = 0, ncol = length(epi_initial))) ## create an empty data table with 0 rows and length(epi_initial) columns
setnames(epi_final, names(epi_initial)) ## then set the column names
epi_final_smoking <- epi_final
epi_final_nonsmoking <- epi_final

for(i in 1:length(oncogene_variants)){
  for(j in 1:length(tsg)){
    epi <- get_epistasis_cesa(cesa_object = cesa, variant = oncogene_variants[i], gene_name = tsg[j], specific_samples = cesa$samples$Unique_Patient_Identifier)[[1]]
    epi_final <- rbindlist(list(epi_final, epi), fill = F)
    
    epi_smoking_wP <- get_epistasis_cesa(cesa_object = cesa_smoking_w_panel, variant = oncogene_variants[i], gene_name = tsg[j], 
                                         specific_samples = cesa_smoking_w_panel$samples[Unique_Patient_Identifier %in% c(smoking_samples, panel_smoking_samples), Unique_Patient_Identifier])[[1]]
    epi_final_smoking <- rbindlist(list(epi_final_smoking, epi_smoking_wP), fill = F)
    
    epi_nonsmoking_wP <- get_epistasis_cesa(cesa_object = cesa_nonsmoking_w_panel, variant = oncogene_variants[i], gene_name = tsg[j], 
                                            specific_samples = cesa_nonsmoking_w_panel$samples[Unique_Patient_Identifier %in% c(nonsmoking_samples, panel_nonsmoking_samples), Unique_Patient_Identifier])[[1]]
    epi_final_nonsmoking <- rbindlist(list(epi_final_nonsmoking, epi_nonsmoking_wP), fill = F) 
    }
}



epi_data <- rbind(epi_final[, data_type := "All"],
                  epi_final_smoking[, data_type := "Smo"],
                  epi_final_nonsmoking[, data_type := "nonSmo"])
epi_data <- epi_data[,c("variant_A_type","remaining_column") := tstrsplit(run_type, split = "_vs", keep = 1)]
epi_data[,gene_Name :=stringr::str_split(run_type, "_", simplify = T)[,4]]

## remove all rows with "CDKN2A.p16INK4a" in "gene_Name" column and "CDKN2A.p14arf" in "variant_B" column   
## because select_variants() function also select all variant information of CDKN2A.p14arf" when the "gene_name" is "CDKN2A.p16INK4a", which cause a duplication.
epi_data <- epi_data[!.("CDKN2A.p16INK4a","CDKN2A.p14arf"), on=.(gene_Name,variant_B)]

## add variant_B_s column, to shorten CDKN2A.p16INK4a and CDKN2A.p14arf to p16INK4a and p14arf
epi_data[,variant_B_s := variant_B]
epi_data[stringr::str_detect(variant_B_s,"\\."),variant_B_s := stringr::str_split(variant_B_s,"\\.", simplify = T)[,2]]


save(epi_data, file = paste0(rdata_output, "epistasis_variant_TSG_output.Rdata"))

## Visualizing the ces_B_on_A #######


# ##  calculate the p value using Welch's t test ##
# t_test = function(mu1, mu2, sd1, sd2, n1, n2, mu0 = 0, detailed=F){
#   t_stat = (mu1 - mu2 - mu0) / sqrt((sd1^2/n1) + (sd2^2/n2))
#   df = ((sd1^2/n1 + sd2^2/n2)^2) / ((sd1^4/(n1^2*(n1 - 1))) + (sd2^4/(n2^2*(n2 - 1))))
#   p_val = 2 * pt(-abs(t_stat), df)
#   
#   if(!detailed){
#     return(p_val)
#   } else{
#     return(list(t_stat, p_val))
#   }
# }
# 
# epistasis_t_test_mia_v2 <-  function(wt_gamma_mle, mut_gamma_mle, wt_gamma_ci_low, wt_gamma_ci_high, mut_gamma_ci_low, mut_gamma_ci_high, wt_n, mut_n){ 
#   ## calculate the standard deviation based on mean,sample size and confidence interval
#   ## use z-distribution for sample size is greater than 30, use t-stat when sample size is smaller than or equal to 30
#   ci_level = 0.95
#   alpha = 1-ci_level
#   Zstat <- qnorm(1- alpha/2) 
#   
#   if (wt_n >30){
#   wt_sd= (wt_gamma_ci_high - wt_gamma_mle)*sqrt(wt_n)/Zstat
#   }else{ #<= 30 
#     df = wt_n - 1
#     t_stat = qt(1- alpha/2, df)
#     wt_sd= (wt_gamma_ci_high - wt_gamma_mle)*sqrt(wt_n)/t_stat
#   }
#   
#   if (mut_n >30){
#     mut_sd= (mut_gamma_ci_high - mut_gamma_mle)*sqrt(mut_n)/Zstat
#   }else{ #<= 30 
#     df = mut_n - 1
#     t_stat = qt(1- alpha/2, df)
#     mut_sd= (mut_gamma_ci_high - mut_gamma_mle)*sqrt(mut_n)/t_stat
#   }
#   return(t_test(wt_gamma_mle, mut_gamma_mle, wt_sd, mut_sd, wt_n, mut_n, detailed = T))
# }
# 
# get_epi_data_t_test <- function(epi_data, tsg){
#   epi_data_t_test <- data.table()
#   for (i in 1:length(tsg)){
#     tar_gene <- tsg[i] 
#     smo_u <- epi_data[.(tar_gene,"Smo"), on = .(variant_B, data_type),ces_B_on_A]
#     nonSmo_u <- epi_data[.(tar_gene,"nonSmo"), on = .(variant_B, data_type),ces_B_on_A]
#     smo_ci_low <- epi_data[.(tar_gene,"Smo"), on = .(variant_B, data_type), ci_low_95_ces_B_on_A]
#     smo_ci_high <- epi_data[.(tar_gene,"Smo"), on = .(variant_B, data_type), ci_high_95_ces_B_on_A]
#     nonSmo_ci_low <- epi_data[.(tar_gene,"nonSmo"), on = .(variant_B, data_type),ci_low_95_ces_B_on_A]
#     nonSmo_ci_high <- epi_data[.(tar_gene,"nonSmo"), on = .(variant_B, data_type),ci_high_95_ces_B_on_A]
#     smo_nAB <- epi_data[.(tar_gene,"Smo"), on = .(variant_B, data_type),nAB]
#     nonSmo_nAB <- epi_data[.(tar_gene,"nonSmo"), on = .(variant_B, data_type),nAB]
#     
#     new_row <- data.table(gene = tar_gene,smo_mean = smo_u, nonSmo_mean = nonSmo_u,
#                           smo_ci_low_mean = smo_ci_low,
#                           smo_ci_high_mean = smo_ci_high,
#                           nonSmo_ci_low_mean = nonSmo_ci_low,
#                           nonSmo_ci_high_mean = nonSmo_ci_high,
#                           smo_nAB_ = smo_nAB,nonSmo_nAB_ = nonSmo_nAB,
#                           t_test_p_mia_v2 = NULL, t_test_t_mia_v2 = NULL)
#     epi_data_t_test <- rbind(epi_data_t_test, new_row, fill = T)
#     if(is.na(smo_nAB) |is.na(nonSmo_nAB)){
#       epi_data_t_test$t_test_p_mia_v2[i] <- "nAB_NA"
#       epi_data_t_test$t_test_t_mia_v2[i] <- "nAB_NA"
#     }else{
#       if(smo_nAB == 0 | nonSmo_nAB == 0){
#         epi_data_t_test$t_test_p_mia_v2[i] <- "nAB_0"
#         epi_data_t_test$t_test_t_mia_v2[i] <- "nAB_0"
#       }else if (smo_nAB == 1 | nonSmo_nAB == 1){
#         epi_data_t_test$t_test_p_mia_v2[i] <- "nAB_1"
#         epi_data_t_test$t_test_t_mia_v2[i] <- "nAB_1"
#       }else{
#         t_test_mia_result_v2 <-  epistasis_t_test_mia_v2(smo_u,nonSmo_u,smo_ci_low,smo_ci_high,nonSmo_ci_low,nonSmo_ci_high,smo_nAB,nonSmo_nAB)
#         epi_data_t_test$t_test_p_mia_v2[i] <- t_test_mia_result_v2[[2]]
#         epi_data_t_test$t_test_t_mia_v2[i] <- t_test_mia_result_v2[[1]]
#       }
#       
#     }
#   }
#   return(epi_data_t_test)
# }

### for KRAS_G12C ####
epi_KRAS_G12C <- epi_data["KRAS_G12C", on = "variant_A_type"]
epi_KRAS_G12C_all <- epi_data[.("KRAS_G12C","All"), on=.(variant_A_type, data_type)]
epi_KRAS_G12C_S <- epi_data[.("KRAS_G12C","Smo"), on=.(variant_A_type, data_type)]
epi_KRAS_G12C_nonS <- epi_data[.("KRAS_G12C","nonSmo"), on=.(variant_A_type, data_type)]

## sorting gene by ces_B_on_A of TSG under KRSA_G12C background in smokers, then neversmokers, then all samples 
variant_order_KRAS_G12C <- epi_KRAS_G12C_S[ces_B_on_A > 0.001000][order(-ces_B_on_A)]$variant_B_s
variant_lf_KRAS_G12C <- epi_KRAS_G12C_S[ces_B_on_A == 0.001000]$variant_B_s
variant_lf_order_KRAS_G12C <- epi_KRAS_G12C_nonS[variant_lf_KRAS_G12C, on = "variant_B"][ces_B_on_A > 0.001000][order(-ces_B_on_A)]$variant_B_s
variant_lf_2_KRAS_G12C <- epi_KRAS_G12C_nonS[variant_lf_KRAS_G12C, on = "variant_B"][ces_B_on_A == 0.001000]$variant_B_s
variant_lf_2_order_KRAS_G12C <- epi_KRAS_G12C_all[variant_lf_2_KRAS_G12C, on = "variant_B"][ces_B_on_A > 0.001000][order(-ces_B_on_A)]$variant_B_s
variant_lf_3_KRAS_G12C <- epi_KRAS_G12C_all[variant_lf_2_KRAS_G12C, on = "variant_B"][ces_B_on_A == 0.001000]$variant_B_s

gene_order_KRAS_G12C <- c(variant_order_KRAS_G12C,variant_lf_order_KRAS_G12C,variant_lf_2_order_KRAS_G12C,variant_lf_3_KRAS_G12C)
gene_order_KRAS_G12C

epi_KRAS_G12C[,variant_B_s := factor(variant_B_s, levels = gene_order_KRAS_G12C, ordered = T )] 
epi_KRAS_G12C[,data_type := factor(data_type, levels = c("All","Smo","nonSmo"), ordered = T )] 

## edit CI bar which has NA values
epi_KRAS_G12C$ci_low_95_ces_B_on_A_forplot <- ifelse(is.na(epi_KRAS_G12C$ci_low_95_ces_B_on_A),epi_KRAS_G12C$ces_B_on_A, epi_KRAS_G12C$ci_low_95_ces_B_on_A)
epi_KRAS_G12C$errorBarColor <- ifelse(epi_KRAS_G12C$nAB %in% c("0","1"), "Gray","Black")
epi_KRAS_G12C$data_type2 <-  as.vector(epi_KRAS_G12C$data_type)
epi_KRAS_G12C$data_type2[epi_KRAS_G12C$nAB == "0"] <- "Gray"

# ## add epi_data_t_test ##
# epi_KRAS_G12C_t_test <- get_epi_data_t_test(epi_KRAS_G12C,tsg)
# epi_KRAS_G12C_t_test$p_value <- ifelse(stringr::str_detect(epi_KRAS_G12C_t_test$t_test_p_mia_v2,"nAB"), NA, as.numeric(epi_KRAS_G12C_t_test$t_test_p_mia_v2)) 
# epi_KRAS_G12C_t_test[,gene_s := gene]
# epi_KRAS_G12C_t_test[stringr::str_detect(gene_s,"\\."),gene_s := stringr::str_split(gene_s,"\\.", simplify = T)[,2]]
# epi_KRAS_G12C_t_test <- epi_KRAS_G12C_t_test[order(match(epi_KRAS_G12C_t_test$gene_s, gene_order_KRAS_G12C))]
# epi_KRAS_G12C_t_test$p_value_anno <- ifelse(is.na(epi_KRAS_G12C_t_test$p_value), "", ifelse(epi_KRAS_G12C_t_test$p_value < 0.01, "**", ifelse(epi_KRAS_G12C_t_test$p_value < 0.05, "*", "ns")))
# setnames(epi_KRAS_G12C_t_test, old = "gene_s", new = "variant_B_s")
# epi_KRAS_G12C_t_test[,variant_B_s := factor(variant_B_s, levels = gene_order_KRAS_G12C, ordered = T)]

p_G12C <- ggplot(data= epi_KRAS_G12C, aes(x= data_type , y = ces_B_on_A))+ 
   geom_errorbar(aes(ymin = ci_low_95_ces_B_on_A_forplot, ymax = ci_high_95_ces_B_on_A, color = errorBarColor), width = 0,  na.rm = T) + 
  geom_point(aes(fill = data_type2, color = errorBarColor),shape = 21, size =2) + 
  scale_fill_manual(values =c("All" = "black", "Smo" = "red","nonSmo" = "blue","Gray" = "gray")) +
  scale_color_manual(values = c("Gray" = "gray","Black" = "black"))+
  labs(x = "Tumor suppressor gene")+
  ylab(bquote(bold(atop("Cancer effect",
                        "(KRAS G12C background)"))))+ 
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
  #   data = epi_KRAS_G12C_t_test, 
  #   aes(x=2.5,y=max(epi_KRAS_G12C$ci_high_95_ces_B_on_A)+1, label = p_value_anno),
  #   hjust = 0.5, vjust = 0.5 ) 

### for KRAS_G12D ####
epi_KRAS_G12D <- epi_data["KRAS_G12D", on = "variant_A_type"]

## use the gene_order_KRAS_G12C
epi_KRAS_G12D[,variant_B_s := factor(variant_B_s, levels = gene_order_KRAS_G12C, ordered = T )] 
epi_KRAS_G12D[,data_type := factor(data_type, levels = c("All","Smo","nonSmo"), ordered = T )] 

## edit CI bar which has NA values
epi_KRAS_G12D$ci_low_95_ces_B_on_A_forplot <- ifelse(is.na(epi_KRAS_G12D$ci_low_95_ces_B_on_A),epi_KRAS_G12D$ces_B_on_A, epi_KRAS_G12D$ci_low_95_ces_B_on_A)
epi_KRAS_G12D$errorBarColor <- ifelse(epi_KRAS_G12D$nAB %in% c("0","1"), "Gray","Black")
epi_KRAS_G12D$data_type2 <-  as.vector(epi_KRAS_G12D$data_type)
epi_KRAS_G12D$data_type2[epi_KRAS_G12D$nAB == "0"] <- "Gray"

# ## add epi_data_t_test ##
# epi_KRAS_G12D_t_test <- get_epi_data_t_test(epi_KRAS_G12D,tsg)
# epi_KRAS_G12D_t_test$p_value <- ifelse(stringr::str_detect(epi_KRAS_G12D_t_test$t_test_p_mia_v2,"nAB"), NA, as.numeric(epi_KRAS_G12D_t_test$t_test_p_mia_v2)) 
# epi_KRAS_G12D_t_test[,gene_s := gene]
# epi_KRAS_G12D_t_test[stringr::str_detect(gene_s,"\\."),gene_s := stringr::str_split(gene_s,"\\.", simplify = T)[,2]]
# epi_KRAS_G12D_t_test <- epi_KRAS_G12D_t_test[order(match(epi_KRAS_G12D_t_test$gene_s, gene_order_KRAS_G12D))]
# epi_KRAS_G12D_t_test$p_value_anno <- ifelse(is.na(epi_KRAS_G12D_t_test$p_value), "", ifelse(epi_KRAS_G12D_t_test$p_value < 0.01, "**", ifelse(epi_KRAS_G12D_t_test$p_value < 0.05, "*", "ns")))
# setnames(epi_KRAS_G12D_t_test, old = "gene_s", new = "variant_B_s")
# epi_KRAS_G12D_t_test[,variant_B_s := factor(variant_B_s, levels = gene_order_KRAS_G12C, ordered = T)]

p_G12D <- ggplot(data= epi_KRAS_G12D, aes(x= data_type , y = ces_B_on_A))+ 
  geom_errorbar(aes(ymin = ci_low_95_ces_B_on_A_forplot, ymax = ci_high_95_ces_B_on_A, color = errorBarColor), width = 0,  na.rm = T) + 
  geom_point(aes(fill = data_type2, color = errorBarColor),shape = 21, size =2) + 
  scale_fill_manual(values =c("All" = "black", "Smo" = "red","nonSmo" = "blue","Gray" = "gray")) +
  scale_color_manual(values = c("Gray" = "gray","Black" = "black"))+
  labs(x = "Tumor suppressor gene")+
  ylab(bquote(bold(atop("Cancer effect",
                        "(KRAS G12D background)"))))+ 
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
#   data = epi_KRAS_G12D_t_test, 
#   aes(x=2.5,y=max(epi_KRAS_G12D$ci_high_95_ces_B_on_A)+1, label = p_value_anno),
#   hjust = 0.5, vjust = 0.5 ) 

### for BRAF_V600E ####
epi_BRAF_V600E <- epi_data["BRAF_V600E", on = "variant_A_type"]

## use the gene_order_KRAS_G12C
epi_BRAF_V600E[,variant_B_s := factor(variant_B_s, levels = gene_order_KRAS_G12C, ordered = T )] 
epi_BRAF_V600E[,data_type := factor(data_type, levels = c("All","Smo","nonSmo"), ordered = T )] 

## edit CI bar which has NA values
epi_BRAF_V600E$ci_low_95_ces_B_on_A_forplot <- ifelse(is.na(epi_BRAF_V600E$ci_low_95_ces_B_on_A),epi_BRAF_V600E$ces_B_on_A, epi_BRAF_V600E$ci_low_95_ces_B_on_A)
epi_BRAF_V600E$errorBarColor <- ifelse(epi_BRAF_V600E$nAB %in% c("0","1"), "Gray","Black")
epi_BRAF_V600E$data_type2 <-  as.vector(epi_BRAF_V600E$data_type)
epi_BRAF_V600E$data_type2[epi_BRAF_V600E$nAB == "0"] <- "Gray"

# ## add epi_data_t_test ##
# epi_BRAF_V600E_t_test <- get_epi_data_t_test(epi_BRAF_V600E,tsg)
# epi_BRAF_V600E_t_test$p_value <- ifelse(stringr::str_detect(epi_BRAF_V600E_t_test$t_test_p_mia_v2,"nAB"), NA, as.numeric(epi_BRAF_V600E_t_test$t_test_p_mia_v2)) 
# epi_BRAF_V600E_t_test[,gene_s := gene]
# epi_BRAF_V600E_t_test[stringr::str_detect(gene_s,"\\."),gene_s := stringr::str_split(gene_s,"\\.", simplify = T)[,2]]
# epi_BRAF_V600E_t_test <- epi_BRAF_V600E_t_test[order(match(epi_BRAF_V600E_t_test$gene_s, gene_order_BRAF_V600E))]
# epi_BRAF_V600E_t_test$p_value_anno <- ifelse(is.na(epi_BRAF_V600E_t_test$p_value), "", ifelse(epi_BRAF_V600E_t_test$p_value < 0.01, "**", ifelse(epi_BRAF_V600E_t_test$p_value < 0.05, "*", "ns")))
# setnames(epi_BRAF_V600E_t_test, old = "gene_s", new = "variant_B_s")
# epi_BRAF_V600E_t_test[,variant_B_s := factor(variant_B_s, levels = gene_order_KRAS_G12C, ordered = T)]

p_V600E <- ggplot(data= epi_BRAF_V600E, aes(x= data_type , y = ces_B_on_A))+ 
  geom_errorbar(aes(ymin = ci_low_95_ces_B_on_A_forplot, ymax = ci_high_95_ces_B_on_A, color = errorBarColor), width = 0,  na.rm = T) + 
  geom_point(aes(fill = data_type2, color = errorBarColor),shape = 21, size =2) + 
  scale_fill_manual(values =c("All" = "black", "Smo" = "red","nonSmo" = "blue","Gray" = "gray")) +
  scale_color_manual(values = c("Gray" = "gray","Black" = "black"))+
  labs(x = "Tumor suppressor gene")+
  ylab(bquote(bold(atop("Cancer effect",
                        "(KRAS G12D background)"))))+ 
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
#   data = epi_BRAF_V600E_t_test, 
#   aes(x=2.5,y=max(epi_BRAF_V600E$ci_high_95_ces_B_on_A)+1, label = p_value_anno),
#   hjust = 0.5, vjust = 0.5 ) 

### for EGFR_L858R ####
epi_EGFR_L858R <- epi_data["EGFR_L858R", on = "variant_A_type"]

## use the gene_order_KRAS_G12C
epi_EGFR_L858R[,variant_B_s := factor(variant_B_s, levels = gene_order_KRAS_G12C, ordered = T )] 
epi_EGFR_L858R[,data_type := factor(data_type, levels = c("All","Smo","nonSmo"), ordered = T )] 

## edit CI bar which has NA values
epi_EGFR_L858R$ci_low_95_ces_B_on_A_forplot <- ifelse(is.na(epi_EGFR_L858R$ci_low_95_ces_B_on_A),epi_EGFR_L858R$ces_B_on_A, epi_EGFR_L858R$ci_low_95_ces_B_on_A)
epi_EGFR_L858R$errorBarColor <- ifelse(epi_EGFR_L858R$nAB %in% c("0","1"), "Gray","Black")
epi_EGFR_L858R$data_type2 <-  as.vector(epi_EGFR_L858R$data_type)
epi_EGFR_L858R$data_type2[epi_EGFR_L858R$nAB == "0"] <- "Gray"

# ## add epi_data_t_test ##
# epi_EGFR_L858R_t_test <- get_epi_data_t_test(epi_EGFR_L858R,tsg)
# epi_EGFR_L858R_t_test$p_value <- ifelse(stringr::str_detect(epi_EGFR_L858R_t_test$t_test_p_mia_v2,"nAB"), NA, as.numeric(epi_EGFR_L858R_t_test$t_test_p_mia_v2)) 
# epi_EGFR_L858R_t_test[,gene_s := gene]
# epi_EGFR_L858R_t_test[stringr::str_detect(gene_s,"\\."),gene_s := stringr::str_split(gene_s,"\\.", simplify = T)[,2]]
# epi_EGFR_L858R_t_test <- epi_EGFR_L858R_t_test[order(match(epi_EGFR_L858R_t_test$gene_s, gene_order_EGFR_L858R))]
# epi_EGFR_L858R_t_test$p_value_anno <- ifelse(is.na(epi_EGFR_L858R_t_test$p_value), "", ifelse(epi_EGFR_L858R_t_test$p_value < 0.01, "**", ifelse(epi_EGFR_L858R_t_test$p_value < 0.05, "*", "ns")))
# setnames(epi_EGFR_L858R_t_test, old = "gene_s", new = "variant_B_s")
# epi_EGFR_L858R_t_test[,variant_B_s := factor(variant_B_s, levels = gene_order_KRAS_G12C, ordered = T)]

p_L858R <- ggplot(data= epi_EGFR_L858R, aes(x= data_type , y = ces_B_on_A))+ 
  geom_errorbar(aes(ymin = ci_low_95_ces_B_on_A_forplot, ymax = ci_high_95_ces_B_on_A, color = errorBarColor), width = 0,  na.rm = T) + 
  geom_point(aes(fill = data_type2, color = errorBarColor),shape = 21, size =2) + 
  scale_fill_manual(values =c("All" = "black", "Smo" = "red","nonSmo" = "blue","Gray" = "gray")) +
  scale_color_manual(values = c("Gray" = "gray","Black" = "black"))+
  labs(x = "Tumor suppressor gene")+
  ylab(bquote(bold(atop("Cancer effect",
                        "(KRAS G12D background)"))))+ 
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
#   data = epi_EGFR_L858R_t_test, 
#   aes(x=2.5,y=max(epi_EGFR_L858R$ci_high_95_ces_B_on_A)+1, label = p_value_anno),
#   hjust = 0.5, vjust = 0.5 ) 



### plot CES_B_on_A of TSG under a specific variant background in one figure  ####
jpeg(file = "ces_B_on_A_variant_TSG.jpeg", width = 13, height = 13, units = 'in', res = 300) 
p_G12C + p_G12D + p_V600E + p_L858R + plot_annotation(tag_levels = 'A') + plot_layout(nrow = 4, guides = "collect") 
dev.off()

### plot CES_B_on_A of TSG under a specific variant background; only show those points with significant difference between smokers and never-smokers ####
# ## use t_test p value as the standards 
# epi_EGFR_L858R_t_test_signifi <- epi_EGFR_L858R_t_test[p_value_anno %in% c("*","**","***")]
# gene_order_signifi <- epi_EGFR_L858R_t_test_signifi[order(smo_mean, decreasing = T),as.character(variant_B_s)]
# epi_EGFR_L858R_t_test_signifi[,p_value_anno_nonOverlap:= "*"]
# epi_EGFR_L858R_t_test_signifi[,variant_B_s := as.character(variant_B_s) ] 
# epi_EGFR_L858R_t_test_signifi[,variant_B_s := factor(variant_B_s, levels = gene_order_signifi, ordered = T)]

## use non-overlapping CIs as the standards
variant_tsg_signifi_gene <- c("TP53")
epi_EGFR_L858R_signfi <- epi_EGFR_L858R[variant_B_s %in% variant_tsg_signifi_gene]
variant_tsg_signifi <- data.frame(variant_tsg_signifi_gene = variant_tsg_signifi_gene,
                                  p_value_anno_nonOverlap = c("*"))

p_L858R_signifi <- ggplot(data= epi_EGFR_L858R_signfi, aes(x= data_type , y = ces_B_on_A))+ 
  geom_errorbar(  aes(ymin = ci_low_95_ces_B_on_A, ymax = ci_high_95_ces_B_on_A, color = errorBarColor), width = 0,  na.rm = T) + ##color = "black", #, linetype = data_type, position = position_jitterdodge(dodge.width = 0.75,jitter.width = 0.4)
  geom_point(aes(fill = data_type2, color = errorBarColor),shape = 21, size =6) + ## position = "jitter", size = 1## position=position_jitterdodge(dodge.width = 0.75,jitter.width = 0.4)
  scale_fill_manual(values =c("All" = "black", "Smo" = "red","nonSmo" = "blue","Gray" = "gray")) +
  scale_color_manual(values = c("Gray" = "gray","Black" = "black"))+
  labs(x = "Tumor suppressor gene")+
  ylab(bquote(bold(atop("Cancer effect", 
                        "(EGFR L858R background)"))))+ 
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = F))+
  scale_x_discrete(labels = c("All","Smoker","Never-smoker"))+ ## change the x-axis label
  theme_minimal()+
  facet_wrap(~variant_B_s, nrow=1)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size=16,face = "bold"),
        axis.text.y = element_text(size =16,face = "bold"),
        axis.title.y = element_text( size = 18,face= "bold"),
        axis.title.x = element_text( size = 18,face = "bold"),
        strip.text.x = element_text(angle = 45, hjust =0.5, vjust = 0.5,size =16, face = "bold"))+
  geom_text(
    data = variant_tsg_signifi, 
    aes(x=2.5,y=max(epi_EGFR_L858R_signfi$ci_high_95_ces_B_on_A)+1, label = p_value_anno_nonOverlap),
    hjust = 0.5, vjust = 0.5, size =6 ) ##size =4,


svg(file = "CES_VariantTSG_signifi_byNonOverlapping.svg", height =9.60, width = 2.93) ## 6inch = 1800 pixels. 1278 pixels, 1000
p_L858R_signifi
dev.off()






## CES value only in All group####
epi_All <- rbind(epi_KRAS_G12C[c("All"), on = "data_type"],epi_KRAS_G12D[c("All"), on = "data_type"],
                            epi_BRAF_V600E[c("All"), on = "data_type"], epi_EGFR_L858R[c("All"), on = "data_type"])
### CES value only in All group; plot CES_B_on_A in all samples####
gene_order_for_All <- epi_KRAS_G12C[c("All"), on = "data_type"][order(-ces_B_on_A)]$variant_B_s
epi_All[,variant_A_type := factor(variant_A_type, levels =c("KRAS_G12C","KRAS_G12D","BRAF_V600E","EGFR_L858R"), ordered = T) ]
epi_All[,variant_B_s := as.character(variant_B_s)]
epi_All[,variant_B_s := factor(variant_B_s, levels =gene_order_for_All, ordered = T) ]
str(epi_All)
ggplot(data= epi_All, aes(x= variant_B_s , y = ces_B_on_A))+ 
  geom_errorbar(  aes(ymin = ci_low_95_ces_B_on_A, ymax = ci_high_95_ces_B_on_A, color = errorBarColor), width = 0,  na.rm = T) + ##color = "black", #, linetype = data_type, position = position_jitterdodge(dodge.width = 0.75,jitter.width = 0.4)
  geom_point(aes(fill = data_type2, color = errorBarColor),shape = 21, size =2) + ## position = "jitter", size = 1## position=position_jitterdodge(dodge.width = 0.75,jitter.width = 0.4)
  scale_fill_manual(values =c("All" = "black", "Gray" = "gray")) + ## "Smo" = "red","nonSmo" = "blue",
  scale_color_manual(values = c("Gray" = "gray","Black" = "black"))+
  labs(x = "Tumor suppressor gene")+
  ylab(bquote(bold(atop("Cancer effect", 
                        "(Specific Oncogenic Variant background)"))))+ 
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = F))+
  #scale_x_discrete(labels = c("All","Smoker","Never-smoker"))+ ## change the x-axis label
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size= 6),
        axis.title.y = element_text(size = 8, face= "bold"),
        axis.title.x = element_text(size = 8, face = "bold"),
        legend.position = "none") +
  facet_wrap(~variant_A_type, nrow=4)+
  theme(strip.text.x = element_text(angle = 0, hjust =0.5, vjust = 0.5,size =7, face = "bold"))

### CES value only in All group; sample size forLRT ####
epi_All_forLRT<- epi_All[,.(variant_A_type, variant_B, ces_B_on_A, nA0, nB0, nAB, n00,data_type)]
View(epi_All_forLRT)
write.csv(epi_All_forLRT, file = "sampleSize_forLRT_TSGinOncovariant_allSamples.csv", row.names = F)
dim(cesa$gene_rates) ## 20091
head(cesa$gene_rates)

## gene-level neutral mutation rates
gene_rates_allSample <- cesa$gene_rates[c(tsg, "KRAS", "BRAF","EGFR"),on= "gene"] ## column 2 name 是 rate，而不是 rate_grp_1
setnames(gene_rates_allSample, old = "rate", new = "rate_grp_1")
gene_rates_smoking <- cesa_smoking_w_panel$gene_rates[c(tsg, "KRAS", "BRAF","EGFR"),on= "gene"]
gene_rates_nonsmoking <- cesa_nonsmoking_w_panel$gene_rates[c(tsg, "KRAS", "BRAF","EGFR"),on= "gene"]
gene_rates <- rbind(gene_rates_allSample[,data_type := "All"],gene_rates_smoking[, data_type := "Smo"],gene_rates_nonsmoking[, data_type := "nonSmo"])
View(gene_rates)
write.csv(gene_rates, file = "gene-level_neutral_mutation_rates.csv", row.names = F)

## baseline mutation rate: neutral mutation rates at specific sites based on gene mutation rates and the relative trinucleotide-context-specific SNV mutation rates of each sample
variants <- c("KRAS_G12C","KRAS_G12D","BRAF_V600E","EGFR_L858R")
## get variant_id
variant_id <- cesa$variants[stringr::str_detect(variant_id,"KRAS_G12C|KRAS_G12D|BRAF_V600E|EGFR_L858R"), variant_id]

baseline_rates_allSample <- baseline_mutation_rates(cesa = cesa, variant_ids = variant_id)
head(baseline_rates_allSample)
baseline_rates_smoking <- baseline_mutation_rates(cesa = cesa_smoking_w_panel, variant_ids = variant_id, 
                                                  samples = cesa_smoking_w_panel$samples[Unique_Patient_Identifier %in% c(smoking_samples, panel_smoking_samples), Unique_Patient_Identifier])
baseline_rates_nonsmoking <- baseline_mutation_rates(cesa = cesa_nonsmoking_w_panel, variant_ids = variant_id, 
                                                     samples = cesa_nonsmoking_w_panel$samples[Unique_Patient_Identifier %in% c(nonsmoking_samples, panel_nonsmoking_samples), Unique_Patient_Identifier])
write.csv(baseline_rates_allSample, file = "baseline_mutation_rates_allSamples.csv", row.names = F)
write.csv(baseline_rates_smoking, file = "baseline_mutation_rates_smoking.csv", row.names = F)
write.csv(baseline_rates_nonsmoking, file = "baseline_mutation_rates_nonsmoking.csv", row.names = F)

## sample size forLRT in smokers and never sokers ####
G12C_samples_forLRT <- epi_KRAS_G12C[,.(variant_A_type, variant_B, nA0, nB0, nAB, n00,data_type)]
View(G12C_samples_forLRT)
setnames(G12C_samples_forLRT,"variant_A_type","variant_A")
G12C_samples_forLRT <- G12C_samples_forLRT[!"All", on="data_type"]
write.csv(G12C_samples_forLRT, file = "G12C_samples_forLRT.csv", row.names = F)

G12D_samples_forLRT <- epi_KRAS_G12D_v4[,.(variant_A_type, variant_B, nA0, nB0, nAB, n00,data_type)]
View(G12D_samples_forLRT)
setnames(G12D_samples_forLRT,"variant_A_type","variant_A")
G12D_samples_forLRT <- G12D_samples_forLRT[!"All", on="data_type"]
write.csv(G12D_samples_forLRT, file = "G12D_samples_forLRT.csv", row.names = F)

V600E_samples_forLRT <- epi_BRAF_V600E[,.(variant_A_type, variant_B, nA0, nB0, nAB, n00,data_type)]
View(V600E_samples_forLRT)
setnames(V600E_samples_forLRT,"variant_A_type","variant_A")
V600E_samples_forLRT <- V600E_samples_forLRT[!"All", on="data_type"]
write.csv(V600E_samples_forLRT, file = "V600E_samples_forLRT.csv", row.names = F)

L858R_samples_forLRT <- epi_EGFR_L858R[,.(variant_A_type, variant_B, nA0, nB0, nAB, n00,data_type)]
View(L858R_samples_forLRT)
setnames(L858R_samples_forLRT,"variant_A_type","variant_A")
L858R_samples_forLRT <- L858R_samples_forLRT[!"All", on="data_type"]
write.csv(V600E_samples_forLRT, file = "L858R_samples_forLRT.csv", row.names = F)












