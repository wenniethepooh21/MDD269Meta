library(here)
library(tidyr)
library(metap)
library(magrittr)
library(dplyr)
library(readr)


source(here("R", "loading datasets", "DingMetaPAnalysis.R"))

fullDingTable <- read_csv(here("ProcessedData", "DingEtAl", "fullDingTableFisherMagma.csv"))

#flip the male effectsizes
fullDingTable %<>% mutate(effectsize = ifelse(sex == "male",effectsize*-1, effectsize))
fullDingTable %>% write_csv( path = here("ProcessedData", "DingEtAl","fullDingTableFisherMagma_flipped.csv")) #Keep as a csv file 

Ding_summary_results <- fullDingTable %>% newDingMeta(8,4)

DingMale <- fullDingTable %>% filter(sex == "male")
Ding_summary_results_Male<- DingMale %>% newDingMeta(4,4)


### writing of slimmed/updated version

#flip the directions back to original
Ding_unfliped <- read_csv(here("ProcessedData", "DingEtAl", "DingTableFisher.slim_magma.csv"))
Ding_unfliped %<>% select(1:2)
Ding_unfliped %<>% right_join(Ding_summary_results %>% select(-2))
Ding_unfliped %>% write_csv(here("ProcessedData", "DingEtAl", "DingTableFisherMagma.slim_flipped.csv"))


#flip the directions back to original
Male_unfliped <- read_csv(here("ProcessedData", "DingEtAl", "DingTableFisher.slim.Male_magma.csv"))
Male_unfliped %<>% select(1:2)
Male_unfliped %<>% right_join(Ding_summary_results_Male %>% select(-2))

Male_unfliped %>% write_csv(path = here("ProcessedData", "DingEtAl", "DingTableFisherMagma.slim.Male_flipped.csv"))


#--------- Ranking Analysis --------------------------------#
fullDing_meta <- read_csv(here("ProcessedData", "DingEtAl", "DingTableFisherMagma.slim_flipped.csv"))
maleDing_meta <- read_csv(here("ProcessedData", "DingEtAl", "DingTableFisherMagma.slim.Male_flipped.csv"))

num_genes_meta <- length(fullDing_meta$gene_symbol)
num_genes_males <- length(maleDing_meta$gene_symbol)

fullDing <- getRankFishers(fullDing_meta, num_genes_meta)
merged_full <- left_join(fullDing_meta, fullDing %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_full %>% write_csv(path = here("ProcessedData", "DingEtAl", "DingTableFisherMagma.slim_flipped.csv"))

maleDing <- getRankFishers(maleDing_meta, num_genes_males)
merged_male <- left_join(maleDing_meta, maleDing %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_male %>% write_csv(path = here("ProcessedData", "DingEtAl", "DingTableFisherMagma.slim.Male_flipped.csv"))


















