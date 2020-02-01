library(here)
library(tidyr)
library(magrittr)
library(dplyr)
library(readr)


source(here("R", "loading datasets", "DingMetaPAnalysis.R"))

Ding <- read_csv(here("ProcessedData", "DingEtAl", "fullDingTableFisherMagma.csv"))

#Remve the data for AMY brain region
Ding_cortical <- Ding %>% filter(brain_region != "AMY")
Ding_cortical %>% write_csv(path = here("ProcessedData", "DingEtAl", "CorticalDingTableFisherMagma.slim.csv"))

#Extract the female data separately from the male data
DingFemale <- Ding_cortical %>% filter(sex == "female")

DingMale <- Ding_cortical %>% filter(sex == "male")


Ding_cortical_results <- newDingMeta(Ding_cortical, 6,3)
Ding_cortical_results %>% write_csv(here("ProcessedData", "DingEtAl", "CorticalDingTableFisherMagma.slim.csv"))

Ding_Female_cortical_results <- newDingMeta(DingFemale,3,3)
Ding_Female_cortical_results %>% write_csv(here("ProcessedData", "DingEtAl", "CorticalDingTableFisherMagma.slim.Female.csv"))

Ding_Male_cortical_results <- newDingMeta(DingMale,3,3)
Ding_Male_cortical_results %>% write_csv(here("ProcessedData", "DingEtAl", "CorticalDingTableFisherMagma.slim.Male.csv"))


#-------------------- CORRECTING FOR INDEPENDENCE ---------------------#
fullDing_cortical <- read_csv(here("ProcessedData", "DingEtAl", "CorticalDingTableFisherMagma.slim.csv"))
femaleDing_cortical <- read_csv(here("ProcessedData", "DingEtAl", "CorticalDingTableFisherMagma.slim.Female.csv"))
maleDing_cortical <- read_csv(here("ProcessedData", "DingEtAl", "CorticalDingTableFisherMagma.slim.Male.csv"))

num_genes_cortical <- length(fullDing_cortical$gene_symbol)
num_genes_female_cortical <- length(femaleDing_cortical$gene_symbol)
num_genes_males_cortical <- length(maleDing_cortical$gene_symbol)

fullCortical <- getRankFishers(fullDing_cortical, num_genes_cortical)
merged_full <- left_join(fullDing_cortical, fullCortical %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_full %>% write_csv(path = here("ProcessedData", "DingEtAl", "CorticalDingTableFisherMagma.slim.csv"))

femaleCortical <- getRankFishers(femaleDing_cortical, num_genes_female_cortical)
merged_female <- left_join(femaleDing_cortical, femaleCortical %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_female %>% write_csv(path = here("ProcessedData", "DingEtAl", "CorticalDingTableFisherMagma.slim.Female.csv"))

maleCortical <- getRankFishers(maleDing_cortical, num_genes_males_cortical)
merged_male <- left_join(maleDing_cortical, maleCortical %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_male %>% write_csv(path = here("ProcessedData", "DingEtAl", "CorticalDingTableFisherMagma.slim.Male.csv"))

