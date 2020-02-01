library(here)
library(tidyr)
library(magrittr)
library(readr)
library(reshape)
library(dplyr)

source(here("R", "loading datasets", "DingMetaPAnalysis.R"))
#From email from Etienne - subject: MetaA-MDD results
# Columns C-E tells you which genes were used and the reasons: Low mean meants too little expression for analysis, low variance meant invariant expression
# Columns G-N are p-values per gene in the individual datasets (See paper for methods)
# Columns P-W are effect sizes per genes in the individual datasets
# Columns Y-Z and AB-AC are the p- and q-values for two different methods of meta-analysis we used (See paper for methods)
# Columns AE-AF and AH-AI are the results of the meta-analysis in the male and female cohorts only
# Columns AK-AL are the results of the meta-regression for effects of sex.
# 
# For the initial analysis, I would use the REM-all p-values, but we should include the genes that were not analyzed due to low variance. These genes were detected but do not show any MDD-related effect (assume p-values ~0.9999)
# We can discuss if any of this is not clear.


fullDingTable <- read_csv(here("data","DingEtAl","MDD-metaAR_8cohorts_Final.csv"))

#remove the columns for the encoded subtable
colnamesDing <- colnames(fullDingTable)
goodCols <- colnamesDing[1:which(colnamesDing =="X39")-1]
fullDingTable %<>% dplyr::select(one_of(goodCols))
fullDingTable %<>% dplyr::select(-starts_with("X"))

colnames(fullDingTable) <- gsub("[.]1", "_effectsize", colnames(fullDingTable))
#use the 10680
fullDingTable %<>% filter(selected == 1)

#correct the dates 
change_date <- fullDingTable %>% filter(grepl("^[[:digit:]]+-", SYMBOL)) %>% select(SYMBOL)
change_date %<>% mutate(month_val = toupper(gsub("^[[:digit:]]+-", "", SYMBOL)))
#only march and sept, add the missing letters to the name 
change_date %<>% mutate(month_val = ifelse(month_val == "MAR","MARCH", "SEPT"))
change_date %<>% mutate(date_val = gsub("-[[:alpha:]]+","", SYMBOL))
change_date %<>% unite(new_symbol, c("month_val", "date_val"), sep = "")

change_date %<>% full_join(fullDingTable, by = c('SYMBOL' = 'SYMBOL'))
change_date %<>% mutate(new_symbol = ifelse(is.na(new_symbol), SYMBOL, new_symbol)) %>% select(-SYMBOL) %>% dplyr::rename(SYMBOL = new_symbol)

fullDingTable <- change_date #table is updated with no dates

newFullDingTable <- fullDingTable %>% processDingTable() %>% write_csv(here("ProcessedData", "DingEtAl", "fullDingTableFisher.csv"))

#combine 8 (4 studies * 2 sexes) 8 pvalues for each gene
Ding_summary_results <- newDingMeta(newFullDingTable,8, 4)
Ding_summary_results %>% write_csv(path = here("ProcessedData", "DingEtAl", "DingTableFisher.slim.csv"))

#get the sex data separately
Ding_Male <- newFullDingTable %>% filter(sex == "male") %>% write_csv(here("ProcessedData", "DingEtAl", "maleDingTableFisher.csv"))
Ding_Female <- newFullDingTable %>% filter(sex == "female") %>% write_csv(here("ProcessedData", "DingEtAl", "femaleDingTableFisher.csv"))

#combine 4 (4 studies * 1 sexes) 
Ding_Female_results <- newDingMeta(Ding_Female,4, 4)
Ding_Female_results %>% write_csv(path = here("ProcessedData", "DingEtAl", "DingTableFisher.slim.Female.csv"))

#combine 4 (4 studies * 1 sexes) 
Ding_Male_results <- newDingMeta(Ding_Male,4,4)
Ding_Male_results %>% write_csv(path = here("ProcessedData", "DingEtAl", "DingTableFisher.slim.Male.csv"))

#--------- Genome Percentile Ranking MAGMA--------------------------------#
magma <- read_csv(here("data", "HowardEtAl", "FullMagmaGenes.csv"))

newFullDingTable %<>% mutate(SYMBOL = gsub("C([X0-9]+)ORF([0-9]+)", "C\\1orf\\2", SYMBOL)) 

Ding_magma <- left_join(magma %>% select(Ding_genes) %>% distinct(), newFullDingTable, by = c('Ding_genes' = 'SYMBOL'))
Ding_magma %<>% dplyr::rename(SYMBOL = Ding_genes) %>% na.omit() %>% write_csv(here("ProcessedData", "DingEtAl", "fullDingTableFisherMagma.csv"))

Ding_summary_magma<- Ding_magma %>% newDingMeta(8, 4)
Ding_summary_magma %>% write_csv(path = here("ProcessedData", "DingEtAl", "DingTableFisher.slim_magma.csv"))

Ding_summary_magma_Female<- Ding_magma %>% filter(sex == "female") %>% newDingMeta(4, 4)
Ding_summary_magma_Female %>% write_csv(path = here("ProcessedData", "DingEtAl", "DingTableFisher.slim.Female_magma.csv"))

Ding_summary_magma_Male<- Ding_magma %>% filter(sex == "male") %>% newDingMeta(4, 4)
Ding_summary_magma_Male %>% write_csv(path = here("ProcessedData", "DingEtAl", "DingTableFisher.slim.Male_magma.csv"))

fullDing_meta <- read_csv(here("ProcessedData", "DingEtAl", "DingTableFisher.slim_magma.csv"))
femaleDing_meta <- read_csv (here("ProcessedData", "DingEtAl", "DingTableFisher.slim.Female_magma.csv"))
maleDing_meta <- read_csv(here("ProcessedData", "DingEtAl", "DingTableFisher.slim.Male_magma.csv"))

num_genes_meta <- length(fullDing_meta$gene_symbol)
num_genes_female <- length(femaleDing_meta$gene_symbol)
num_genes_males <- length(maleDing_meta$gene_symbol)

fullDing <- getRankFishers(fullDing_meta, num_genes_meta)
merged_full <- left_join(fullDing_meta, fullDing %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_full %>% write_csv(path = here("ProcessedData", "DingEtAl", "DingTableFisher.slim_magma.csv"))

femaleDing <- getRankFishers(femaleDing_meta, num_genes_female)
merged_female <- left_join(femaleDing_meta, femaleDing %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_female %>% write_csv(path = here("ProcessedData", "DingEtAl", "DingTableFisher.slim.Female_magma.csv"))

maleDing <- getRankFishers(maleDing_meta, num_genes_males)
merged_male <- left_join(maleDing_meta, maleDing %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_male %>% write_csv(path = here("ProcessedData", "DingEtAl", "DingTableFisher.slim.Male_magma.csv"))
