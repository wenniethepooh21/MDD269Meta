library(magrittr)
library(tidyr)
library(here)
#library(biomaRt)
library(dplyr)
library(readxl)
library(readr)
                                                          #This script runs all meta-analyses on the Labonte ranscriptomic dataset 
# detach("package:here", unload=TRUE)
# setwd('../../../school/thesis/')
# library(here)

#Read in raw male expression data provided by Labonte, et al. 
Labonte_Male <- read_xlsx(here("Raw_Data/LabonteEtAl/Male DEG all.xlsx"))
Labonte_Male %<>% rename(gene_symbol  = "Gene name", brain_region = "Brain Region", Male.pvalue = "p-value", Male.logFC = logFC, Male.AveExpr = AveExpr)

#Read in raw female expression data provided by Labonte, et al. 
Labonte_Female <- read_xlsx(here("Raw_Data/LabonteEtAl/Female DEG all.xlsx"))
Labonte_Female %<>% rename(brain_region = `Brain region`, gene_symbol = `Gene Name`, Female.pvalue = `p-value`, Female.logFC = logFC, Female.AveExpr = AveExpr)

Labonte <- left_join(Labonte_Male %>% select(-Biotype, -Description), Labonte_Female %>% select(-Biotype, -Description), by = c("brain_region", "ENSG", "gene_symbol"))

########################################### update excel converted gene symbols from dates to gene symbols using ENSG number################
#convert ENSG to gene gene_symbol to remove dates from gene_symbols column
#ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

#This line was previously run to generate associated gene name, output was saved to save time
#ENSG_Genes <- getBM(attributes = c('ensembl_gene_id','hgnc_gene_symbol'), filters = 'ensembl_gene_id', values = unique(Labonte$ENSG), mart = ensembl) %>% as_tibble() 
#ENSG_Genes %>% write_csv(here("Processed_Data/LabonteEtAl/getBMGenes.csv"))

ENSG_Genes <- read_csv(here("Processed_Data/LabonteEtAl/getBMGenes.csv")) # read in ENSG to hgnc symbol mapping file previously generated for Labonte's genes
Labonte_genes <- Labonte %>% select(ENSG, gene_symbol)
Labonte_genes %<>% left_join(ENSG_Genes, by = c('ENSG' = 'ensembl_gene_id'))

#merge gene gene_symbols into one column 
Labonte_genes %<>% rowwise() %>% mutate(hgnc_symbol = if_else(is.na(hgnc_symbol), gene_symbol, hgnc_symbol))
Labonte_genes %<>% select(-gene_symbol)	
Labonte_genes %<>% distinct()
#Find if all ENSG's are distinct (manual)
Labonte_genes %>% group_by(ENSG) %>% filter(n() > 1)
#What should ENSG00000230417 map to?
Labonte %>% filter(ENSG == "ENSG00000230417") %>% select(gene_symbol)
#keep the hgnc symbol == "LINC00856", remove for "LINC00595":  ENSG00000230417 == LINC00856
Labonte_genes %<>% mutate(hgnc_symbol = if_else(ENSG == "ENSG00000230417" & hgnc_symbol == "LINC00595", "LINC00856", hgnc_symbol)) %>% distinct() 
Labonte %<>% left_join(Labonte_genes)
Labonte %<>% select(-gene_symbol)
Labonte %<>% rename(gene_symbol = "hgnc_symbol")
#re-order columns
Labonte <- Labonte[,c(1:2,9,3:8)]
#########################################################################################################################################################
#remove space between the words
Labonte %<>% mutate(brain_region = ifelse(brain_region == "Anterior Insula", "Anterior_Insula",brain_region))

#filter transcriptomic dataset for MAGMA genes tested by Howard, et al. 
magma_table <- read_csv(here("Raw_Data/HowardEtAl/FullMagmaGenes.csv")) %>% select(Labonte_genes) %>% distinct() %>% na.omit()
magma_labonte <- Labonte %>% right_join(magma_table , by = c('gene_symbol' = 'Labonte_genes'))

#Expand table row wise to hold male and female data 
Labonte_long <- bind_rows(magma_labonte %>% select(brain_region, gene_symbol, logFC = Male.logFC, pvalue = Male.pvalue) %>% mutate(sex = "male"),
                          magma_labonte %>% select(brain_region, gene_symbol, logFC = Female.logFC, pvalue = Female.pvalue) %>% mutate(sex = "female"))
Labonte_long %>% write_csv(here("Processed_Data/LabonteEtAl/CompleteLabonteTableMagma.csv"))
###########################################################
###### REGULAR META-ANALYSIS (FULL, FEMALE AND MALE) ######
###########################################################
#This file holds all functions used for meta-analyses - call this file
source(here("R/transcriptomic_meta/Labonte_Meta_Analysis.R"))
#Perform meta-analysis on all brain regions across both sexes
#number of p-values we need to filter = 12 (6 regions * 2 sexes) 
#number of brain regions = 6
Labonte_summary_results <- Labonte_long %>% LabonteMetaAnalysis(12,6)
Labonte_summary_results %>% write_csv(here("Processed_Data/LabonteEtAl/FullLabonteTableMagma.csv"))

#Perform the sex-specific meta-analyses
#Extract Female data 
Labonte_Female_Magma <- magma_labonte %>% select(brain_region, gene_symbol, logFC = Female.logFC, pvalue = Female.pvalue) %>% mutate(sex = "female")
#Perform meta-analysis on all brain regions in female data 
#number of p-values we need to filter = 6 (6 regions * 1 sex) 
#number of brain regions = 6
Labonte_Female_results <- Labonte_Female_Magma %>% LabonteMetaAnalysis(6,6)
Labonte_Female_results %>% write_csv(here("Processed_Data/LabonteEtAl/FemaleLabonteTableMagma.csv"))

#Extract Male data
Labonte_Male_Magma <- magma_labonte %>% select(brain_region, gene_symbol, logFC = Male.logFC, pvalue = Male.pvalue) %>% mutate(sex = "male")
#Perform meta-analysis on all brain regions in female data 
#number of p-values we need to filter = 6 (6 regions * 1 sex) 
#number of brain regions = 6
Labonte_Male_results <- Labonte_Male_Magma %>% LabonteMetaAnalysis(6,6)
Labonte_Male_results %>% write_csv(path = here("Processed_Data/LabonteEtAl/MaleLabonteTableMagma.csv"))


###########################################
###### CORTICAL ANALYSIS BOTH SEXES ######
###########################################
#Pick out only the cortical regions that Labonte studied
Labonte_Cortical <- Labonte_long %>% filter(brain_region != "Nac") %>% filter( brain_region != "Subic")
#Perform meta-analysis on all brain regions across both sexes
#number of p-values we need to filter = 8 (4 regions * 2 sexes) 
#number of brain regions = 4
cortical_Labonte_summary_results <- Labonte_Cortical %>% LabonteMetaAnalysis(8,4)
cortical_Labonte_summary_results %>% write_csv(path = here("Processed_Data/LabonteEtAl/CorticalLabonteTableMagma.csv"))


#######################################
####SEX-INTERACTION FULL ANALYSIS ####
#######################################
#Flip male logFC value
Labonte_long_flip <- Labonte_long %>% mutate(logFC = if_else(sex == "male",logFC *-1, logFC))
#Perform meta-analysis on all brain regions across both sexes
#number of p-values we need to filter = 12 (6 regions * 2 sexes) 
#number of brain regions = 6
Labonte_full_flip <- Labonte_long_flip %>% LabonteMetaAnalysis(12,6)
#change directions back to original
full_labonte_dir <- Labonte_summary_results %>% select(1:2)
Labonte_summary_full_flip <- full_labonte_dir %>% left_join(Labonte_full_flip %>% select(-2))
Labonte_summary_full_flip %>% write_csv(here("Processed_Data/LabonteEtAl/FullLabonteTableMagma_flipped.csv"))

############################################
#### SEX-INTERACTION CORTICAL ANALYSIS ####
############################################
Labonte_long_cortical_flip <- Labonte_Cortical %>% mutate(logFC = if_else(sex == "male",logFC *-1, logFC))
#Perform meta-analysis on all brain regions across both sexes
#number of p-values we need to filter = 8 (4 regions * 2 sexes) 
#number of brain regions = 4
Labonte_cortical_flip <-Labonte_long_cortical_flip %>% LabonteMetaAnalysis(8,4)
#change directions back to original
cortical_labonte_dir <- cortical_Labonte_summary_results %>% select(1:2)
Labonte_summary_cortical_flip <- cortical_labonte_dir %>% left_join(Labonte_cortical_flip %>% select(-2))
Labonte_summary_cortical_flip %>% write_csv(here("Processed_Data/LabonteEtAl/CorticalLabonteTableMagma_flipped.csv"))


###############################################################
###### GENOME PERCENTILE RANKING FOR ALL ABOVE ANALYSES ######
##############################################################
#load script that holds the genome percentile ranking function used by all transcriptomic studies
source(here("R/transcriptomic_meta/Percentile_Rank_Analysis.R"))

#Run genome percentile ranking analysis - calculates the percentage of genes that have a smaller meta p-value than the current gene on our full meta-analysis
Labonte_summary_results %<>% getRank()
Labonte_summary_results %>% write_csv(here("Processed_Data/LabonteEtAl/FullLabonteTableMagma.csv"))

#Run genome percentile ranking analysis - calculates the percentage of genes that have a smaller meta p-value than the current gene on our female meta-analysis
Labonte_Female_results %<>% getRank()
Labonte_Female_results %>% write_csv(path = here("Processed_Data/LabonteEtAl/FemaleLabonteTableMagma.csv"))

#Run genome percentile ranking analysis - calculates the percentage of genes that have a smaller meta p-value than the current gene on our male meta-analysis
Labonte_Male_results %<>% getRank()
Labonte_Male_results %>% write_csv(path = here("Processed_Data/LabonteEtAl/MaleLabonteTableMagma.csv"))

#Run genome percentile ranking analysis - calculates the percentage of genes that have a smaller meta p-value than the current gene on our cortical meta-analysis
cortical_Labonte_summary_results %<>% getRank()
cortical_Labonte_summary_results %>% write_csv(here("Processed_Data/LabonteEtAl/CorticalLabonteTableMagma.csv"))

#Run genome percentile ranking analysis - calculates the percentage of genes that have a smaller meta p-value than the current gene on our sex-interaction full meta-analysis
Labonte_summary_full_flip %<>% getRank()
Labonte_summary_full_flip %>% write_csv(here("Processed_Data/LabonteEtAl/FullLabonteTableMagma_flipped.csv"))

#Run genome percentile ranking analysis - calculates the percentage of genes that have a smaller meta p-value than the current gene on our sex-interaction cortical meta-analysis
Labonte_summary_cortical_flip %<>% getRank()
Labonte_summary_cortical_flip %>% write_csv(here("Processed_Data/LabonteEtAl/CorticalLabonteTableMagma_flipped.csv"))
