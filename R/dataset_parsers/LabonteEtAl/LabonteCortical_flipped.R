library(magrittr)
library(dplyr)
library(tidyverse)
library(here)
library(metap)


source(here("R", "loading datasets", "LabonteMetaPAnalysis.R"))

Labonte <- read_csv(here("ProcessedData", "LabonteEtAl","CompleteLabonteTableMagma.csv"))

#Pick out only the cortical regions that Labonte studied
LabonteCortical <- Labonte %>% filter(brain_region != "Nac") %>% filter( brain_region != "Subic")

#Flip male logFC value
LabonteCortical %<>% mutate(Male.logFC = Male.logFC *-1)

Labonte_long <- bind_rows(LabonteCortical %>% dplyr::select(brain_region, symbol, logFC = Male.logFC, pvalue = Male.pvalue) %>% mutate(sex = "male"),
                          LabonteCortical %>% dplyr::select(brain_region, symbol, logFC = Female.logFC, pvalue = Female.pvalue) %>% mutate(sex = "female"))

#combine 8 (4 regions * 2 sexes) 
Labonte_summary_results <- LabonteMeta(Labonte_long,8,4)

#change directions back to original

#unflipped directions 
labonte_unfliped <- read_csv(here("ProcessedData", "LabonteEtAl", "CorticalLabonteTableMagma.csv"))
labonte_unfliped %<>% select(1:2)
labonte_unfliped %<>% right_join(Labonte_summary_results %>% select(-2))

#needs writing out
labonte_unfliped %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "CorticalLabonteTableMagma_flipped.csv"))

#flip direction in male data
Labonte_Male <- LabonteCortical %>% select(brain_region, symbol, logFC = Male.logFC, pvalue = Male.pvalue) %>% mutate(sex = "male")

#combine 4 (4 regions * 1 sexes) 
Labonte_Male_results <- LabonteMeta(Labonte_Male,4,4)
#change directions back to original
male_unfliped <- read_csv(here("ProcessedData", "LabonteEtAl", "CorticalMaleLabonteTableMagma.csv"))
male_unfliped %<>% select(1:2)
male_unfliped %<>% right_join(Labonte_Male_results %>% select(-2))

#needs writing out
male_unfliped %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "CorticalMaleLabonteTableMagma_flipped.csv"))



#---------------------------------- Genome Percentile Ranking

#Get Labonte cortical meta analysis data 
fullLabonte_cortical <- read_csv(here("ProcessedData", "LabonteEtAl", "CorticalLabonteTableMagma_flipped.csv"))
maleLabonte_cortical <- read_csv(here("ProcessedData", "LabonteEtAl", "CorticalMaleLabonteTableMagma_flipped.csv"))

num_genes_cortical <- length(fullLabonte_cortical$gene_symbol)
num_genes_cortical_males <- length(maleLabonte_cortical$gene_symbol)


fullCortical <- getRank(fullLabonte_cortical, num_genes_cortical)
merged_full <- left_join(fullLabonte_cortical, fullCortical %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_full %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "CorticalLabonteTableMagma_flipped.csv"))

maleCortical <- getRank(maleLabonte_cortical, num_genes_cortical_males)
merged_male <- left_join(maleLabonte_cortical, maleCortical %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_male %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "CorticalMaleLabonteTableMagma_flipped.csv"))

