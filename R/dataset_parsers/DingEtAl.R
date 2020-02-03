library(here)
library(tidyr)
library(magrittr)
library(readr)
library(dplyr)


# Columns C-E tells you which genes were used and the reasons: Low mean meants too little expression for analysis, low variance meant invariant expression
# Columns G-N are p-values per gene in the individual datasets (See paper for methods)
# Columns P-W are effect sizes per genes in the individual datasets
# Columns Y-Z and AB-AC are the p- and q-values for two different methods of meta-analysis we used (See paper for methods)
# Columns AE-AF and AH-AI are the results of the meta-analysis in the male and female cohorts only
# Columns AK-AL are the results of the meta-regression for effects of sex.
# 
# For the initial analysis, I would use the REM-all p-values, but we should include the genes that were not analyzed due to low variance. These genes were detected but do not show any MDD-related effect (assume p-values ~0.9999)

#read in data provided by Dr. Sibille for this study
DingTable <- read_csv(here("Raw_Data/DingEtAl/MDD-metaAR_8cohorts_Final.csv"))

#remove the columns for the encoded subtable
colnamesDing <- colnames(DingTable)
goodCols <- colnamesDing[1:which(colnamesDing =="X39")-1]
DingTable %<>% select(one_of(goodCols))
DingTable %<>% select(-starts_with("X"))
colnames(DingTable) <- gsub("[.]1", "_effectsize", colnames(DingTable))
#use the 10680 genes that has p-values and effectsizes 
#in the selected column filter for genes that were selected for analysis in the study
DingTable %<>% filter(selected == 1)

#correct the dates to gene symbols
Ding_no_date_table <- DingTable %>% filter(grepl("^[[:digit:]]+-", SYMBOL)) %>% select(SYMBOL)
Ding_no_date_table %<>% mutate(month_val = toupper(gsub("^[[:digit:]]+-", "", SYMBOL)))
#only march and sept, add the missing letters to the name 
Ding_no_date_table %<>% mutate(month_val = ifelse(month_val == "MAR","MARCH", "SEPT"))
Ding_no_date_table %<>% mutate(date_val = gsub("-[[:alpha:]]+","", SYMBOL))
Ding_no_date_table %<>% unite(new_symbol, c("month_val", "date_val"), sep = "")
DingTable %<>% full_join(Ding_no_date_table, by = c('SYMBOL' = 'SYMBOL'))
DingTable %<>% mutate(new_symbol = ifelse(is.na(new_symbol), SYMBOL, new_symbol)) %>% select(-SYMBOL) 
DingTable %<>% rename(gene_symbol = new_symbol)

###########################################################
###### REGULAR META-ANALYSIS (FULL, FEMALE AND MALE) ######
###########################################################
#This file holds all functions used for meta-analyses and other preprocessing of data - call this file
source(here("R/transcriptomic_meta/Ding_Meta_Analysis.R"))
fullDingTable <- DingTable %>% ProcessDingTable() 
#Read in the associated list of MAGMA genes that howard tested (17,842)
magma_table <- read_csv(here("Raw_Data/HowardEtAl/FullMagmaGenes.csv")) %>% select(Ding_genes) %>% distinct() %>% na.omit()
#filter ding table for magma genes
fullDingMagma <- fullDingTable %>% right_join(magma_table, by = c("gene_symbol" = "Ding_genes")) 
fullDingMagma %>% write_csv(here("Processed_Data/DingEtAl/CompleteDingTableMagma.csv"))

#Perform meta-analysis on all brain regions across both sexes
#number of p-values we need to filter = 8 (4 regions * 2 sexes) 
#number of brain regions = 4
Ding_summary_results <- fullDingMagma %>% DingMetaAnalysis(8, 4)
Ding_summary_results %>% write_csv(path = here("Processed_Data/DingEtAl/FullDingTableMagma.csv"))

#Perform meta-analysis on all brain regions in males
#number of p-values we need to filter = 4 (4 regions * 1 sex) 
#number of brain regions = 4
Ding_Male <- fullDingMagma %>% filter(sex == "male") 
#combine 4 (4 studies * 1 sexes) 
Ding_Male_results <- Ding_Male %>% DingMetaAnalysis(4,4)
Ding_Male_results %>% write_csv(path = here("Processed_Data/DingEtAl/MaleDingTableMagma.csv"))

#Perform meta-analysis on all brain regions in females
#number of p-values we need to filter = 4 (4 regions * 1 sex) 
#number of brain regions = 4
Ding_Female <- fullDingMagma %>% filter(sex == "female") 
#combine 4 (4 studies * 1 sexes) 
Ding_Female_results <- Ding_Female %>% DingMetaAnalysis(4, 4)
Ding_Female_results %>% write_csv(path = here("Processed_Data/DingEtAl/FemaleDingTableMagma.csv"))

#########################################
###### CORTICAL META-ANALYSIS  ######
#########################################

#Remve the data for AMY brain region
Ding_cortical <- fullDingMagma %>% filter(brain_region != "AMY")
#Perform meta-analysis on cortical brain regions in both sexes
#number of p-values we need to filter = 6 (3 regions * 2 sex) 
#number of brain regions = 3
Ding_cortical_results <- Ding_cortical %>% DingMetaAnalysis(6,3)
Ding_cortical_results %>% write_csv(here("Processed_Data/DingEtAl/CorticalDingTableMagma.csv"))

#######################################
####SEX-INTERACTION FULL ANALYSIS ####
#######################################
#Flip male effectsize values
fullDingMagma_flip <- fullDingMagma %>% rowwise() %>% mutate(effectsize = if_else(sex == "male", effectsize*-1, effectsize))
full_Ding_flip_results <- fullDingMagma_flip %>% DingMetaAnalysis(8,4)
#get the original directions of expressions 
original_ding_dir <- Ding_summary_results %>% select(1:2)
full_Ding_flip_results_summary <- original_ding_dir %>% left_join(full_Ding_flip_results %>% select(-2))
full_Ding_flip_results_summary %>% write_csv(here("Processed_Data/DingEtAl/FullDingTableMagma_flipped.csv"))

#######################################
####SEX-INTERACTION CORTICAL ANALYSIS ####
#######################################
Ding_cortical_flip <- Ding_cortical %>% rowwise() %>% mutate(effectsize = if_else(sex == "male", effectsize*-1, effectsize))
cortical_Ding_flip_results <- Ding_cortical_flip %>% DingMetaAnalysis(6,3)
original_ding_cortical_dir <- Ding_cortical_results %>% select(1:2)
cortical_Ding_flip_results_summary <- original_ding_cortical_dir %>% left_join(cortical_Ding_flip_results %>% select(-2))
cortical_Ding_flip_results_summary %>% write_csv(here("Processed_Data/DingEtAl/CorticalDingTableMagma_flipped.csv"))


###############################################################
###### GENOME PERCENTILE RANKING FOR ALL ABOVE ANALYSES ######
##############################################################
#load script that holds the genome percentile ranking function used by all transcriptomic studies
source(here("R/transcriptomic_meta/Percentile_Rank_Analysis.R"))

#Run genome percentile ranking analysis - calculates the percentage of genes that have a smaller meta p-value than the current gene on our full meta-analysis
Ding_summary_results %<>% getRank()
Ding_summary_results %>% write_csv(here("Processed_Data/DingEtAl/FullDingTableMagma.csv"))

#Run genome percentile ranking analysis - calculates the percentage of genes that have a smaller meta p-value than the current gene on our female meta-analysis
Ding_Female_results %<>% getRank()
Ding_Female_results %>% write_csv(path = here("Processed_Data/DingEtAl/FemaleDingTableMagma.csv"))

#Run genome percentile ranking analysis - calculates the percentage of genes that have a smaller meta p-value than the current gene on our male meta-analysis
Ding_Male_results %<>% getRank()
Ding_Male_results %>% write_csv(path = here("Processed_Data/DingEtAl/MaleDingTableMagma.csv"))

#Run genome percentile ranking analysis - calculates the percentage of genes that have a smaller meta p-value than the current gene on our cortical meta-analysis
Ding_cortical_results %<>% getRank()
Ding_cortical_results %>% write_csv(here("Processed_Data/DingEtAl/CorticalDingTableMagma.csv"))

#Run genome percentile ranking analysis - calculates the percentage of genes that have a smaller meta p-value than the current gene on our sex-interaction full meta-analysis
full_Ding_flip_results_summary %<>% getRank()
full_Ding_flip_results_summary %>% write_csv(here("Processed_Data/DingEtAl/FullDingTableMagma_flipped.csv"))

#Run genome percentile ranking analysis - calculates the percentage of genes that have a smaller meta p-value than the current gene on our sex-interaction cortical meta-analysis
cortical_Ding_flip_results_summary %<>% getRank()
cortical_Ding_flip_results_summary %>% write_csv(here("Processed_Data/DingEtAl/CorticalDingTableMagma_flipped.csv"))
