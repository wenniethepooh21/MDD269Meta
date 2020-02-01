library(here)
library(tidyr)
library(magrittr)
library(dplyr)

source(here("R", "loading datasets", "RamakerFullMetaAnalysis.R"))

metadata <- read_tsv(here("data", "RamakerEtAl", "GSE80655 from biojupies", "GSE80655-metadata.txt.txt"))
regions <- unique(metadata$`brain region`)

#get female directions across all brain regions
full_female_results <- read_csv(here("ProcessedData", "RamakerEtAl", "FemaleRamakerTable_magma.csv"))
female_summary_results <- read_csv(here("ProcessedData", "RamakerEtAl", "CompleteFemaleRamakerTableMagma.csv"))

#add sex for unique identifier 
female_summary_results %<>% mutate(sex = "female")

male_summary_results <- read_csv(here("ProcessedData", "RamakerEtAl", "CompleteMaleRamakerTableMagma.csv"))
male_summary_results_flip <- male_summary_results %>% mutate(t = t*-1)
write_csv(male_summary_results_flip, path = here("ProcessedData", "RamakerEtAl", "CompleteMaleRamakerTableMagma_flipped.csv"))

#add sex for unique identifier
male_summary_results_flip %<>% mutate(sex = "male")

#merge male summary results with female summary results into one table to perform calculations
summary <- rbind(male_summary_results_flip, female_summary_results)
write_csv(summary, path = here("ProcessedData", "RamakerEtAl", "CompleteRamakerTableMagma_flipped.csv"))

#Calculations from female data and male flipped data combined
summary %<>% RamakerAnalysis(regions)
summary %<>% flipDirections("male")

#Calculate the new male data flipped
male_summary_results_flip %<>% RamakerAnalysis(regions)
male_summary_results_flip %<>% flipDirections("male") #flip the male directions 

male_summary_results_flip %<>% rename(AnCg_nAcc_DLPFC_Male_directions = AnCg_nAcc_DLPFC_directions)


#combine male directions and female directions 
full_summary <- left_join(full_female_results %>% select(gene_symbol, AnCg_nAcc_DLPFC_Female_directions), male_summary_results_flip %>% select(gene_symbol, AnCg_nAcc_DLPFC_Male_directions))
#concatenate the two directions columns
full_summary %<>% unite("AnCg.F_nAcc.F_DLPFC.F_AnCg.M_nAcc.M_DLPFC.M",AnCg_nAcc_DLPFC_Female_directions, AnCg_nAcc_DLPFC_Male_directions, sep = "")
#Combine directions with meta values calculated 
summary %<>% select(-sex, -AnCg_nAcc_DLPFC_directions)
#join calculations for each gene
full_summary %<>% left_join(summary) %>% distinct()

write_csv(full_summary, path = here("ProcessedData", "RamakerEtAl", "fullRamakerTableMagma_flipped.csv"))

#for excel and google sheet viewing, add to beginning of directions string
male_summary_results_flip %<>% select(-sex)
male_summary_results_flip %>% write_csv(path = here("ProcessedData", "RamakerEtAl", "MaleRamakerTableMagma_flipped.csv"))



#---------------------- Genome Ranking
# 
# #get the meta-p value from the analysis
fullRamaker_meta <- read_csv(here("ProcessedData", "RamakerEtAl", "fullRamakerTableMagma_flipped.csv"))
maleRamaker_meta <- read_csv(here("ProcessedData", "RamakerEtAl", "MaleRamakerTableMagma_flipped.csv"))


num_genes_meta <- length(fullRamaker_meta$gene_symbol)
num_genes_males <- length(maleRamaker_meta$gene_symbol)


fullRamaker <- getRank(fullRamaker_meta, num_genes_meta)
merged_full <- left_join(fullRamaker_meta, fullRamaker %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_full %>% write_csv(path = here("ProcessedData", "RamakerEtAl", "fullRamakerTableMagma_flipped.csv"))


maleRamaker <- getRank(maleRamaker_meta, num_genes_males)
merged_male <- left_join(maleRamaker_meta, maleRamaker %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_male %>% write_csv(path = here("ProcessedData", "RamakerEtAl", "MaleRamakerTableMagma_flipped.csv"))

