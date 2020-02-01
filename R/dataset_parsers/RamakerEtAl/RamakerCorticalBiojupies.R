library(here)
library(tidyr)
library(magrittr)
library(dplyr)
library(readr)

source(here("R", "loading datasets", "RamakerFullMetaAnalysis.R"))

metadata <- read_tsv(here("data", "RamakerEtAl", "GSE80655 from biojupies", "GSE80655-metadata.txt.txt"))
regions <- unique(metadata$`brain region`)
regions <- regions[which(regions != "nAcc")]

ramaker_cortical<- read_csv(here("ProcessedData", "RamakerEtAl", "CompleteRamakerTable.csv"))
female_ramaker_cortical <- read_csv(here("ProcessedData", "RamakerEtAl", "CompleteFemaleRamakerTable.csv"))
male_ramaker_cortical<- read_csv(here("ProcessedData", "RamakerEtAl", "CompleteMaleRamakerTable.csv"))

#Extract the cortical region data & run analysis
ramaker_cortical %<>% filter(target_region != "nAcc")
write_csv(ramaker_cortical, path = here("ProcessedData", "RamakerEtAl", "CorticalRamakerTable.csv"))
ramaker_cortical %<>% RamakerAnalysis(regions)

#Extract the cortical region data & run analysis
female_ramaker_cortical %<>% filter(target_region != "nAcc")
write_csv(female_ramaker_cortical, path = here("ProcessedData", "RamakerEtAl", "CorticalFemaleRamakerTable.csv"))
female_ramaker_cortical %<>% RamakerAnalysis(regions)
female_ramaker_cortical %<>% rename(AnCg_DLPFC_Female_directions = AnCg_DLPFC_directions)

#Extract the cortical region data & run analysis
male_ramaker_cortical %<>% filter(target_region != "nAcc")
write_csv(male_ramaker_cortical, path = here("ProcessedData", "RamakerEtAl", "CorticalMaleRamakerTable.csv"))
male_ramaker_cortical %<>% RamakerAnalysis(regions)
male_ramaker_cortical %<>% rename(AnCg_DLPFC_Male_directions = AnCg_DLPFC_directions)

#merge all gender directions into summary_results table for better visualization
#summary <- summary_results %>% select(gene_symbol, AnCg_nAcc_DLPFC_directions)

summary <- left_join(female_ramaker_cortical %>% select(gene_symbol,AnCg_DLPFC_Female_directions), male_ramaker_cortical %>% select(gene_symbol,AnCg_DLPFC_Male_directions ))
summary %<>% unite("AnCg.F_DLPFC.F_AnCg.M_DLPFC.M", AnCg_DLPFC_Female_directions, AnCg_DLPFC_Male_directions, sep = "")
summary %<>% left_join(ramaker_cortical) %>% select(-AnCg_DLPFC_directions) %>% distinct()

#Save data in .csv files
summary %>% write_csv(here("ProcessedData", "RamakerEtAl", "fullCorticalRamakerTable.csv"))
female_ramaker_cortical %>% write_csv(here("ProcessedData", "RamakerEtAl", "fullCorticalFemaleRamakerTable.csv"))
male_ramaker_cortical %>% write_csv(here("ProcessedData", "RamakerEtAl", "fullCorticalMaleRamakerTable.csv"))


#---------------Genome Ranking -------------------#

fullRamaker_cortical <- read_csv(here("ProcessedData", "RamakerEtAl", "fullCorticalRamakerTable.csv"))
femaleRamaker_cortical <- read_csv(here("ProcessedData", "RamakerEtAl", "fullCorticalFemaleRamakerTable.csv"))
maleRamaker_cortical <- read_csv(here("ProcessedData", "RamakerEtAl", "fullCorticalMaleRamakerTable.csv"))

num_genes_cortical <- length(fullRamaker_cortical$gene_symbol)
num_genes_female_cortical <- length(femaleRamaker_cortical$gene_symbol)
num_genes_males_cortical <- length(maleRamaker_cortical$gene_symbol)

fullCortical <- getRank(fullRamaker_cortical, num_genes_cortical)
merged_full <- left_join(fullRamaker_cortical, fullCortical %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_full %>% write_csv(path = here("ProcessedData", "RamakerEtAl", "fullCorticalRamakerTable.csv"))

femaleCortical <- getRank(femaleRamaker_cortical, num_genes_female_cortical)
merged_female <- left_join(femaleRamaker_cortical, femaleCortical %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_female %>% write_csv(path = here("ProcessedData", "RamakerEtAl", "fullCorticalFemaleRamakerTable.csv"))

maleCortical <- getRank(maleRamaker_cortical, num_genes_males_cortical)
merged_male <- left_join(maleRamaker_cortical, maleCortical %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_male %>% write_csv(path = here("ProcessedData", "RamakerEtAl", "fullCorticalMaleRamakerTable.csv"))

