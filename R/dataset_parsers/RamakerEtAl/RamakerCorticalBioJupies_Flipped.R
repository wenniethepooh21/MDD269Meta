library(here)
library(tidyr)
library(magrittr)
library(dplyr)
library(readr)

#Flip the male directions for the cortical data 
source(here("R", "loading datasets", "RamakerFullMetaAnalysis.R"))

metadata <- read_tsv(here("data", "RamakerEtAl", "GSE80655 from biojupies", "GSE80655-metadata.txt.txt"))
regions <- unique(metadata$`brain region`)
regions <- regions[which(regions != "nAcc")]

#get female directions across all brain regions
full_female_results <- read_csv(here("ProcessedData", "RamakerEtAl", "fullCorticalFemaleRamakerTable.csv"))
#get female data for calculations
female_summary_results <- read_csv(here("ProcessedData", "RamakerEtAl", "CorticalFemaleRamakerTable.csv"))
#add sex for unique identifier 
female_summary_results %<>% mutate(sex = "female")

#get male data for calculations
male_summary_results <- read_csv(here("ProcessedData", "RamakerEtAl", "CorticalMaleRamakerTable.csv"))
#flip the male direction
male_summary_results_flip <- male_summary_results %>% mutate(t = t*-1)
write_csv(male_summary_results_flip, path = here("ProcessedData", "RamakerEtAl", "CorticalMaleRamakerTable_flipped.csv"))
#add sex for unique identifier
male_summary_results_flip %<>% mutate(sex = "male")

#merge male summary results with female summary results into one table to perform calculations
summary <- rbind(male_summary_results_flip, female_summary_results)
write_csv(summary, path = here("ProcessedData", "RamakerEtAl", "CorticalRamakerTable_flipped.csv"))

#Calculations from female data and male flipped data combined
summary %<>% RamakerAnalysis(regions)

#Calculate the new male cortical data flipped
male_summary_results_flip %<>% RamakerAnalysis(regions)
male_summary_results_flip %<>% rename(AnCg_DLPFC_Male_directions = AnCg_DLPFC_directions)


#combine male directions and female directions 
full_summary <- left_join(full_female_results %>% select(gene_symbol, AnCg_DLPFC_Female_directions), male_summary_results_flip %>% select(gene_symbol, AnCg_DLPFC_Male_directions))
#concatenate the two directions columns
full_summary %<>% unite("AnCg.F_DLPFC.F_AnCg.M_DLPFC.M", AnCg_DLPFC_Female_directions, AnCg_DLPFC_Male_directions, sep = "")
#Combine directions with meta values calculated 
#join calculations for each gene
full_summary %<>% left_join(summary) %>% select(-sex, -AnCg_DLPFC_directions) %>% distinct()

write_csv(full_summary, path = here("ProcessedData", "RamakerEtAl", "fullCorticalRamakerTable_flipped.csv"))

#for excel and google sheet viewing, add to beginning of directions string
male_summary_results_flip %<>% select(-sex)
male_summary_results_flip %>% write_csv(path = here("ProcessedData", "RamakerEtAl", "fullCorticalMaleRamakerTable_flipped.csv"))



#---------------CORRECTING FOR INDEPENDENCE -------------------#

fullRamaker_cortical <- read_csv(here("ProcessedData", "RamakerEtAl", "fullCorticalRamakerTable_flipped.csv"))
maleRamaker_cortical <- read_csv(here("ProcessedData", "RamakerEtAl", "fullCorticalMaleRamakerTable_flipped.csv"))

num_genes_cortical <- length(fullRamaker_cortical$gene_symbol)
num_genes_males_cortical <- length(maleRamaker_cortical$gene_symbol)

fullCortical <- getRank(fullRamaker_cortical, num_genes_cortical)
merged_full <- left_join(fullRamaker_cortical, fullCortical %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_full %>% write_csv(path = here("ProcessedData", "RamakerEtAl", "fullCorticalRamakerTable_flipped.csv"))


maleCortical <- getRank(maleRamaker_cortical, num_genes_males_cortical)
merged_male <- left_join(maleRamaker_cortical, maleCortical %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_male %>% write_csv(path = here("ProcessedData", "RamakerEtAl", "fullCorticalMaleRamakerTable_flipped.csv"))

