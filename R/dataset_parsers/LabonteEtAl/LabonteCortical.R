library(magrittr)
library(dplyr)
library(tidyverse)
library(here)
library(metap)


source(here("R", "loading datasets", "LabonteMetaPAnalysis.R"))

Labonte <- read_csv(here("ProcessedData", "LabonteEtAl","CompleteLabonteTableMagma.csv"))

#Pick out only the cortical regions that Labonte studied
LabonteCortical <- Labonte %>% filter(brain_region != "Nac") %>% filter( brain_region != "Subic")

Labonte_long <- bind_rows(LabonteCortical %>% dplyr::select(brain_region, symbol, logFC = Male.logFC, pvalue = Male.pvalue) %>% mutate(sex = "male"),
LabonteCortical %>% dplyr::select(brain_region, symbol, logFC = Female.logFC, pvalue = Female.pvalue) %>% mutate(sex = "female"))

#combine 8 (4 regions * 2 sexes) 
Labonte_summary_results <- LabonteMeta(Labonte_long,8,4)

#needs writing out
Labonte_summary_results %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "CorticalLabonteTableMagma.csv"))

Labonte_Male <- LabonteCortical %>% dplyr::select(brain_region, symbol, logFC = Male.logFC, pvalue = Male.pvalue) %>% mutate(sex = "male")
Labonte_Female <- LabonteCortical %>% dplyr::select(brain_region, symbol, logFC = Female.logFC, pvalue = Female.pvalue) %>% mutate(sex = "female")

#combine 4 (4 regions * 1 sexes) 
Labonte_Female_results <- LabonteMeta(Labonte_Female, 4,4)

#combine 4 (4 regions * 1 sexes) 
Labonte_Male_results <- LabonteMeta(Labonte_Male,4,4)

Labonte_Female_results %>% write_csv(path= here("ProcessedData", "LabonteEtAl", "CorticalFemaleLabonteTableMagma.csv"))
Labonte_Male_results %>% write_csv(path= here("ProcessedData", "LabonteEtAl", "CorticalMaleLabonteTableMagma.csv"))

#---------------------------------- Genome Percentile Ranking

#Get Labonte cortical meta analysis data 
fullLabonte_cortical <- read_csv(here("ProcessedData", "LabonteEtAl", "CorticalLabonteTableMagma.csv"))
femaleLabonte_cortical <- read_csv(here("ProcessedData", "LabonteEtAl", "CorticalFemaleLabonteTableMagma.csv"))
maleLabonte_cortical <- read_csv(here("ProcessedData", "LabonteEtAl", "CorticalMaleLabonteTableMagma.csv"))

num_genes_cortical <- length(fullLabonte_cortical$gene_symbol)
num_genes_cortical_female <- length(femaleLabonte_cortical$gene_symbol)
num_genes_cortical_males <- length(maleLabonte_cortical$gene_symbol)


fullCortical <- getRank(fullLabonte_cortical, num_genes_cortical)
merged_full <- left_join(fullLabonte_cortical, fullCortical %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_full %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "CorticalLabonteTableMagma.csv"))

femaleCortical <- getRank(femaleLabonte_cortical, num_genes_cortical_female)
merged_female <- left_join(femaleLabonte_cortical, femaleCortical %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_female %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "CorticalFemaleLabonteTableMagma.csv"))

maleCortical <- getRank(maleLabonte_cortical, num_genes_cortical_males)
merged_male <- left_join(maleLabonte_cortical, maleCortical %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_male %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "CorticalMaleLabonteTableMagma.csv"))

