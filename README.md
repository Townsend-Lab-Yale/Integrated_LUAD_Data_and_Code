## Integrated_LUAD_Data_and Code
### Integrated LUAD Data
A comprehensive LUAD dataset of 9230 samples with whole-genome sequence data (WGS), whole-exome sequence data (WES), and targeted sequence data (TGS), including 1066 smokers and 656 never-smokers. 
### Code
All necessary codes for the manuscript are included here. 
#### System Requirements
All codes have been tested on MacOS with R version 4.3.0.
#### Installation Guide
Users should install the following packages to use the code, from an R terminal
```
install.packages(c('ggplot2', 'data.table', 'dplyr', 'rtracklayer', 'stringr', 'ggpubr', 'patchwork'))
remotes::install_github("Townsend-Lab-Yale/cancereffectsizeR@v2.10.1", dependencies = T, force = T)
remotes::install_github("Townsend-Lab-Yale/ces.refset.hg19@*release", dependencies = T, force = T) 
```
#### Demo 
Expected output files can be found in R_data/epistasis_Oncogene_TSG_output.Rdata and R_data/epistasis_variant_TSG_output.Rdata
#### Instructions for Use
First, run `01_CES_variants_calculation_and_plot.R` file to process all maf files, calculate the cancer effect size information for all samples, smokers, and never-smokers, and also do the visualization.
Then, run `02_epistasis_variant_TSG.R` file to perform epistasis analysis between a variant and a gene, and also do the visualization.
More detailed annotations can be found in the code files.
#### License
All codes are licensed under GNU General Public License v3.0





