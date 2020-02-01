library(magrittr)
library(tidyverse)
library(here)
library(dplyr)


source(here("R", "loading datasets", "LabonteMetaPAnalysis.R"))

#Import Labonte data 
Labonte <- read_csv(here("ProcessedData", "LabonteEtAl", "CompleteLabonteTableMagma.csv"))

#Flip male logFC value
Labonte %<>% mutate(Male.logFC = Male.logFC *-1)

#Write Data to .csv files
Labonte %>% write_csv(path = here("ProcessedData", "LabonteEtAl","CompleteLabonteTableMagma_flipped.csv"))

Labonte_long <- bind_rows(Labonte %>% dplyr::select(brain_region, symbol, logFC = Male.logFC, pvalue = Male.pvalue) %>% mutate(sex = "male"),
                          Labonte %>% dplyr::select(brain_region, symbol, logFC = Female.logFC, pvalue = Female.pvalue) %>% mutate(sex = "female"))

#combine 12 (6 regions * 2 sexes) + 12 pvalues for each gene
Labonte_summary_results <- LabonteMeta(Labonte_long,12,6)
#change directions back to original

#unflipped directions 
labonte_unfliped <- read_csv(here("ProcessedData", "LabonteEtAl", "fullLabonteTable.csv"))
labonte_unfliped %<>% select(1:2)
labonte_unfliped %<>% right_join(Labonte_summary_results %>% select(-2))

#needs writing out
labonte_unfliped %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "fullLabonteTableMagma_flipped.csv"))

#flip direction in male data
Labonte_Male <- Labonte %>% select(brain_region, symbol, logFC = Male.logFC, pvalue = Male.pvalue) %>% mutate(sex = "male")

#combine 6 (6 regions * 1 sexes) 
Labonte_Male_results <- LabonteMeta(Labonte_Male,6,6)
#change directions back to original
male_unfliped <- read_csv(here("ProcessedData", "LabonteEtAl", "MaleLabonteTable.csv"))
male_unfliped %<>% select(1:2)
male_unfliped %<>% right_join(Labonte_Male_results %>% select(-2))

#needs writing out
male_unfliped %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "MaleLabonteTableMagma_flipped.csv"))

#------- get rankings

#get the meta-p value from the analysis
fullLabonte_meta <- read_csv(here("ProcessedData", "LabonteEtAl", "fullLabonteTableMagma_flipped.csv"))
maleLabonte_meta <- read_csv(here("ProcessedData", "LabonteEtAl", "MaleLabonteTableMagma_flipped.csv"))

num_genes_meta <- length(fullLabonte_meta$gene_symbol)
num_genes_males <- length(maleLabonte_meta$gene_symbol)

fullLabonte <- getRank(fullLabonte_meta, num_genes_meta)
merged_full <- left_join(fullLabonte_meta, fullLabonte %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_full %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "fullLabonteTableMagma_flipped.csv"))

maleLabonte <- getRank(maleLabonte_meta, num_genes_males)
merged_male <- left_join(maleLabonte_meta, maleLabonte %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_male %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "MaleLabonteTableMagma_flipped.csv"))

















