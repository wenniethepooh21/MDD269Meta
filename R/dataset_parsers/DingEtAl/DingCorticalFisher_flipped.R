library(here)
library(tidyr)
library(magrittr)
library(dplyr)
library(readr)

source(here("R", "loading datasets", "DingMetaPAnalysis.R"))

Ding <- read_csv(here("ProcessedData", "DingEtAl","fullDingTableFisherMagma_flipped.csv"))

Ding_cortical <- Ding %>% filter(brain_region != "AMY")
#Extract the female data separately from the male data
DingMale <- Ding_cortical %>% filter(sex == "male")

Ding_cortical_results <- newDingMeta(Ding_cortical,6,3)

DingCortUnflip <- read_csv(here("ProcessedData", "DingEtAl", "CorticalDingTableFisherMagma.slim.csv"))
DingCortUnflip %<>% select(1:2)
DingCortUnflip %<>% right_join(Ding_cortical_results %>% select(-2))
DingCortUnflip %>% write_csv(here("ProcessedData", "DingEtAl", "CorticalDingTableFisherMagma.slim_flipped.csv"))


Ding_Male_cortical_results <- newDingMeta(DingMale,3,3)
DingMaleCortUnflip <- read_csv(here("ProcessedData", "DingEtAl", "CorticalDingTableFisherMagma.slim.Male.csv"))
DingMaleCortUnflip %<>% select(1:2)
DingMaleCortUnflip %<>% right_join(Ding_Male_cortical_results %>% select(-2))
DingMaleCortUnflip %>% write_csv(here("ProcessedData", "DingEtAl", "CorticalDingTableFisherMagma.slim.Male_flipped.csv"))



#-------------------- CORRECTING FOR INDEPENDENCE ---------------------#
fullDing_cortical <- read_csv(here("ProcessedData", "DingEtAl", "CorticalDingTableFisherMagma.slim_flipped.csv"))
maleDing_cortical <- read_csv(here("ProcessedData", "DingEtAl", "CorticalDingTableFisherMagma.slim.Male_flipped.csv"))

num_genes_cortical <- length(fullDing_cortical$gene_symbol)
num_genes_males_cortical <- length(maleDing_cortical$gene_symbol)

fullCortical <- getRankFishers(fullDing_cortical, num_genes_cortical)
merged_full <- left_join(fullDing_cortical, fullCortical %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_full %>% write_csv(path = here("ProcessedData", "DingEtAl", "CorticalDingTableFisherMagma.slim_flipped.csv"))

maleCortical <- getRankFishers(maleDing_cortical, num_genes_males_cortical)
merged_male <- left_join(maleDing_cortical, maleCortical %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_male %>% write_csv(path = here("ProcessedData", "DingEtAl", "CorticalDingTableFisherMagma.slim.Male_flipped.csv"))
