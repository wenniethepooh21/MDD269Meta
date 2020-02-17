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
# max_four_donors_reannotated %>% write_csv(here("Processed_Data/AllenEtAl/max_expression_in_four_donors.csv"))

howard <- read_csv(here("Processed_Data/HowardEtAl/fullHowardTable.csv")) 
howard_genes <- howard %>% select(gene_symbol) %>% pull()
#filter the brain expression data for the 269
brain_expression <- max_four_donors_reannotated %>% filter(gene_symbol %in% howard_genes)
# brain_expression %<>% gather("brain_region", "expression", 2:ncol(brain_expression))

#Howard's updated gene names don't map to any brain regions -- all regions can be found using gene_symbol (don't need to do extra comparison)
howard %<>% select(gene_symbol, Updated_Gene_Names)

Howard_Table <- howard %>% left_join(brain_expression %>% select(gene_symbol, structure_name), by = c('gene_symbol' = 'gene_symbol'))
Howard_Table %<>% rename(brain_region = structure_name)

Howard_Table %>% filter(is.na(brain_region))
#found 3 missing genes = PQLC2L, PRR34, C7orf72 update manually, C7orf72 is not found
#PRR34 has an old gene symbol = C22orf26 found in the expression matrix, add it in
Howard_Table %<>% mutate(brain_region = if_else(gene_symbol == "PRR34", brain_expression %>% filter(gene_symbol == "C22orf26") %>% select(structure_name) %>% as.character(), brain_region))
#PQLC2L old gene symbol is C3orf55
Howard_Table %<>% mutate(brain_region = if_else(gene_symbol == "PQLC2L", brain_expression %>% filter(gene_symbol == "C3orf55") %>% select(structure_name) %>% as.character(), brain_region))

#Map brain structures to ambiguous brain regions 
enclosing_regions <- read_csv(here("Processed_Data/AllenEtAl/full_brain_region_hierarchy.csv"))
# brain_structures<-brain_structures[!duplicated(brain_structures$brain_region),] #remove duplicates 

#remove the comma's in the region name
Howard_Table %<>% mutate(brain_region = gsub(",","",brain_region))
#remove additional spaces 
Howard_Table %<>% mutate(brain_region = gsub("\\s+", " ", brain_region))
Howard_Table %<>% left_join(enclosing_regions, by = c('brain_region' = 'brain_region'))

brain_slim <- read_csv(here("Processed_Data/AllenEtAl/brain_regions_slim.csv"))
Howard_Table %<>% left_join(brain_slim, by = c('brain_region' = 'BrainRegion'))
Howard_Table %<>% mutate(location = if_else(is.na(location), region_location, location)) %>% rename(slim_region_location = location)

Howard_Table %>% write_csv(here("Processed_Data/HowardEtAl/HowardRegions_four.csv"))

# ##############################calculate the probability for the brain regions
brain_pop <- nrow(max_four_donors_reannotated)
structure_count <- max_four_donors_reannotated %>% group_by(structure_name) %>% summarize(genome_tissue_count = n())
#
howard_brain <- Howard_Table %>% select(gene_symbol, Updated_Gene_Names,brain_region,slim_region_location)
howard_brain %<>% mutate(brain_region = ifelse(brain_region == "character(0)", NA, brain_region))
howard_count <- howard_brain %>% group_by(brain_region) %>% summarise(sample_tissue_count = n())

full_count <- structure_count %>% left_join(howard_count, by = c('structure_name' = 'brain_region'))
# perform hypergeometric test - read in file for hyper_test function
source(here("R/transcriptomic_meta/hyper_test.R"))
tissue_expected_probs <- full_count %>% rowwise() %>% mutate(hypergeometric_p = hyper_test(sample_tissue_count, genome_tissue_count, brain_pop, colSums(na.omit(howard_count)[,2])))
#
#correct by number of possible brain structures to choose from
tissue_expected_probs %<>% mutate(corrected_hypergeometric_p = p.adjust(hypergeometric_p, method = "bonferroni", n = nrow(structure_count)))
tissue_expected_probs %<>% arrange(hypergeometric_p)


#upload to google drive
sheets_auth(token = drive_token())

region <- drive_get("~/Thesis/Manuscript/Supplement_Tables/tissue_hyper_expected_four")
if(nrow(region) != 0) {
  drive_rm(region)
}
#create the google worksheet
region <- sheets_create("tissue_hyper_expected_four",sheets = c('hypergeometric_brain_regions'))
sheets_write(tissue_expected_probs, region,  sheet = "hypergeometric_brain_regions")

drive_mv(file = "tissue_hyper_expected_four", path = "~/Thesis/Manuscript/Supplement_Tables/")  # move Sheets file


