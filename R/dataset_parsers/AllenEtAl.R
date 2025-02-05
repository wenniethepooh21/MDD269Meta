library(here)
library(readr)
library(magrittr)
library(dplyr)

#This script aims to identify the brain region that maximally expresses the 269 genes 

#Read in the gene expression data in the brain 
four_donors_reannotated <- read_csv(here("Processed_Data/AllenEtAl/6_brains_reannotated_aggregated_4_donors.csv")) # file provided by Derek Howard
# identify the maximum cell-type taxon that expresses each human gene
# function is in this file - source it
source(here("R/transcriptomic_meta/Max_Expression_Matrix.R"))
max_four_donors_reannotated <- get_max_expression(four_donors_reannotated, "structure_name", "gene_symbol")

#remove the comma's in the region name
max_four_donors_reannotated %<>% mutate(structure_name = gsub(",","",structure_name))
#remove additional spaces 
max_four_donors_reannotated %<>% mutate(structure_name = gsub("\\s+", " ", structure_name))


howard <- read_csv(here("Processed_Data/HowardEtAl/fullHowardTable.csv")) 
howard_genes <- howard %>% dplyr::select(gene_symbol) %>% dplyr::pull()
#filter the brain expression data for the 269
brain_expression <- max_four_donors_reannotated %>% filter(gene_symbol %in% howard_genes)
#howard %<>% dplyr::select(gene_symbol, Updated_Gene_Names)
howard %<>% dplyr::select(gene_symbol, gene_name)

Howard_Table <- howard %>% left_join(brain_expression %>% dplyr::select(gene_symbol, structure_name), by = c('gene_symbol' = 'gene_symbol'))
Howard_Table %<>% dplyr::rename(brain_region = structure_name)

Howard_Table %>% filter(is.na(brain_region))
#found 3 missing genes = PQLC2L, PRR34, C7orf72 update manually, C7orf72 is not found
#PRR34 has an old gene symbol = C22orf26 found in the expression matrix, add it in
Howard_Table %<>% mutate(brain_region = if_else(gene_symbol == "PRR34", max_four_donors_reannotated %>% filter(gene_symbol == "C22orf26") %>% dplyr::select(structure_name) %>% as.character(), brain_region))
#PQLC2L old gene symbol is C3orf55
Howard_Table %<>% mutate(brain_region = if_else(gene_symbol == "PQLC2L", max_four_donors_reannotated %>% filter(gene_symbol == "C3orf55") %>% dplyr::select(structure_name) %>% as.character(), brain_region))

#Map brain structures to ambiguous brain regions 
# Full list of enclosing brain regions
#enclosing_regions_2 <- read_csv(here("Processed_Data/AllenEtAl/full_brain_region_hierarchy.csv"))
enclosing_regions <- read_csv(here("Processed_Data/AllenEtAl/full_brain_region_hierarchy.csv"))

Howard_Table %<>% left_join(enclosing_regions, by = c('brain_region' = 'brain_region'))
#slimmed list of enclosing brain regions
brain_slim <- read_csv(here("Processed_Data/AllenEtAl/brain_regions_slim.csv"))
Howard_Table %<>% left_join(brain_slim, by = c('brain_region' = 'BrainRegion'))
Howard_Table %<>% mutate(location = if_else(is.na(location), region_location, location)) %>% dplyr::rename(slim_region_location = location)

Howard_Table %>% write_csv(here("Processed_Data/HowardEtAl/HowardRegions_four.csv"))

# ##############################calculate the probability for the brain regions
brain_pop <- nrow(max_four_donors_reannotated)
structure_count <- max_four_donors_reannotated %>% group_by(structure_name) %>% summarize(genome_tissue_count = n())
#
howard_brain <- Howard_Table %>% dplyr::select(gene_symbol, gene_name,brain_region,slim_region_location)
howard_count <- howard_brain %>% group_by(brain_region) %>% summarise(sample_tissue_count = n()) %>% na.omit()

full_count <- structure_count %>% full_join(howard_count, by = c('structure_name' = 'brain_region'))
full_count %<>% rowwise() %>% mutate(genome_tissue_count = if_else(is.na(genome_tissue_count) & !is.na(sample_tissue_count), 0, as.numeric(genome_tissue_count))) 
# perform hypergeometric test - read in file for hyper_test function
source(here("R/transcriptomic_meta/hyper_test.R"))
tissue_expected_probs <- full_count %>% rowwise() %>% mutate(expected_tissue_count = genome_tissue_count/brain_pop * colSums(na.omit(howard_count)[,2]))
tissue_expected_probs %<>% mutate(hypergeometric_p = hyper_test(sample_tissue_count, genome_tissue_count, brain_pop, colSums(na.omit(howard_count)[,2])))

#correct by number of possible brain structures to choose from
tissue_expected_probs %<>% mutate(corrected_hypergeometric_p = p.adjust(hypergeometric_p, method = "bonferroni", n = nrow(structure_count)))
tissue_expected_probs %<>% arrange(hypergeometric_p,-sample_tissue_count)
tissue_expected_probs %<>% left_join(brain_slim, by = c('structure_name'='BrainRegion')) %>% ungroup()
tissue_expected_probs %<>% left_join(enclosing_regions, by = c('structure_name' = 'brain_region'))%>% ungroup()
tissue_expected_probs %<>% mutate(location = if_else(is.na(location), region_location, location)) %>% dplyr::rename(enclosing_regions = location)

tissue_expected_probs$hypergeometric_p <- signif(as.numeric(tissue_expected_probs$hypergeometric_p),digits=3)
tissue_expected_probs$corrected_hypergeometric_p <- signif(as.numeric(tissue_expected_probs$corrected_hypergeometric_p),digits=3)

write_csv(tissue_expected_probs, path = here('Results', 'supplementary_tables', 'hypergeometric_brain_regions_tissue_hyper_expected_four_full.csv'))
write_csv(tissue_expected_probs %>% dplyr::select(-region_location) , path = here('Results', 'supplementary_tables', 'hypergeometric_brain_regions_tissue_hyper_expected_four.csv'))

## fix probe mappings to multiple protocadherin genes 
howard_count_fix <- howard_brain %>% filter(gene_symbol != 'PCDHA4' & gene_symbol != 'PCDHA7') %>% group_by(brain_region) %>% summarise(sample_tissue_count = n()) %>% na.omit()
full_count_fix <- structure_count %>% full_join(howard_count_fix, by = c('structure_name' = 'brain_region'))
full_count_fix %<>% rowwise() %>% mutate(genome_tissue_count = if_else(is.na(genome_tissue_count) & !is.na(sample_tissue_count), 0, as.numeric(genome_tissue_count))) 
tissue_expected_probs_corrected <- full_count_fix %>% rowwise() %>% mutate(expected_tissue_count = genome_tissue_count/brain_pop * colSums(na.omit(howard_count)[,2]))
tissue_expected_probs_corrected %<>% mutate(hypergeometric_p = hyper_test(sample_tissue_count, genome_tissue_count, brain_pop, colSums(na.omit(howard_count_fix)[,2])))
#
#correct by number of possible brain structures to choose from
tissue_expected_probs_corrected %<>% mutate(corrected_hypergeometric_p = p.adjust(hypergeometric_p, method = "bonferroni", n = nrow(structure_count)))
tissue_expected_probs_corrected %<>% arrange(hypergeometric_p,-sample_tissue_count)
tissue_expected_probs_corrected %<>% left_join(brain_slim, by = c('structure_name'='BrainRegion')) %>% ungroup()
tissue_expected_probs_corrected %<>% left_join(enclosing_regions, by = c('structure_name' = 'brain_region'))%>% ungroup()
tissue_expected_probs_corrected %<>% mutate(location = if_else(is.na(location), region_location, location)) %>% dplyr::rename(enclosing_regions = location)

tissue_expected_probs_corrected$hypergeometric_p <- signif(as.numeric(tissue_expected_probs_corrected$hypergeometric_p),digits=3)
tissue_expected_probs_corrected$corrected_hypergeometric_p <- signif(as.numeric(tissue_expected_probs_corrected$corrected_hypergeometric_p),digits=3)

write_csv(tissue_expected_probs_corrected, path = here('Results', 'supplementary_tables', 'hypergeometric_brain_regions_tissue_hyper_expected_four_full_corr.csv'))
write_csv(tissue_expected_probs_corrected %>% dplyr::select(-region_location) , path = here('Results', 'supplementary_tables', 'hypergeometric_brain_regions_tissue_hyper_expected_four_corr.csv'))

