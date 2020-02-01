library(mygene)
library(tidyverse)
library(magrittr)
library(dplyr)
library(homologene)
#Using Howard 102 gene list 

#Import the Howard data 
source(here("R", "loading datasets", "fileManipulation.R"))
source(here("R", "loading datasets", "geneDirection.R"))
source(here("R", "loading datasets","ZeiselPolygenic.R"))

GeneReferences <-read_excel(here("Gene_Name_Reference.xlsx"))

HowardGenes <- read_excel(here("data", "HowardEtAl", "433367-22.xlsx"), skip = 1)

HowardGenes %<>% dplyr::select(`Gene Name`, `P-value`)
HowardGenes %<>% dplyr::rename(Howard.Genes = `Gene Name`, Howard.pvalue = `P-value`)
HowardGenes %<>% left_join(GeneReferences %>% dplyr::select(Howard.Gene,Updated_Gene_Names), by = c('Howard.Genes' = 'Howard.Gene'))

#Get the full gene name from the gene symbol
getGeneName <- queryMany(HowardGenes$Updated_Gene_Names, scope = "symbol", species = "human") %>% as_tibble() %>% distinct(query,name)

HowardGenes %<>% left_join(getGeneName, by = c('Updated_Gene_Names' = 'query'))

HowardGenes %<>% dplyr::rename(gene_name = name)

HowardGenes %<>%left_join(GeneReferences %>% dplyr:: select(Howard.Gene, Ramaker.Gene,Labonte.Gene,Ding.Gene, `Match Across Studies?`), by = c('Howard.Genes' = 'Howard.Gene'))

HowardGenes %>% write_csv( path = here("ProcessedData", "HowardEtAl","fullHowardTable.csv"))


#----------------------------------------------
#these files created in R/data_processing/max_expression_matrix.R
HowardGenes <- read_csv(here("ProcessedData", "HowardEtAl", "fullHowardTable.csv"))
#six_brains <- read_csv(here("ProcessedData", "AllenEtAl", "max_expression_in_six_donors.csv"))
six_brains <- read_csv(here("ProcessedData", "AllenEtAl", "max_expression_in_four_donors.csv"))

cell_types <- read_csv(here("ProcessedData", "ZeiselEtAl", "max_cell_type_expression.csv"))
cns_cell_types <- read_csv(here("ProcessedData", "ZeiselEtAl", "max_cns_cell_type_expression.csv"))

#Howard's updated gene names don't map to any brain regions -- all regions can be found using Howard.Genes (don't need to do extra comparison)
HowardGenes %<>% select(Howard.Genes, Updated_Gene_Names)
#HowardGenes %<>% group_by(Howard.Genes) %>% mutate(BrainRegion = six_brains %>% filter(gene_symbol == Howard.Genes) %>% select(structure_name) %>% as.character())
HowardGenes %<>% left_join(six_brains %>% select(gene_symbol, structure_name), by = c('Howard.Genes' = 'gene_symbol'))
HowardGenes %<>% rename(BrainRegion = structure_name)

#found 3 missing genes = PQLC2L, PRR34, C7orf72
#PRR34 has an old gene symbol = C22orf26 found in the expression matrix, add it in
HowardGenes %<>% mutate(BrainRegion = ifelse(Howard.Genes == "PRR34", six_brains %>% filter(gene_symbol == "C22orf26") %>% select(structure_name) %>% as.character(), BrainRegion))
#PQLC2L old gene symbol is C3orf55
HowardGenes %<>% mutate(BrainRegion = ifelse(Howard.Genes == "PQLC2L", six_brains %>% filter(gene_symbol == "C3orf55") %>% select(structure_name) %>% as.character(), BrainRegion))


#Howard genes must be converted to mouse 
#Howard.Genes and Updated_Gene_Names mapped to the same mouse genes, no need for double comparison
howard_mouse_genes <- human2mouse(HowardGenes$Howard.Genes) %>% as_tibble()
howard_mouse_genes %<>% mutate(humanGene = toupper(humanGene))
HowardGenes %<>% mutate(Howard.Genes_upper = toupper(Howard.Genes))
HowardGenes %<>% left_join(howard_mouse_genes %>% select(humanGene, mouseGene), by = c('Howard.Genes_upper' = 'humanGene'))
HowardGenes %<>% left_join(cns_cell_types %>% select(Gene,cell_type_taxon), by = c('mouseGene' = 'Gene'))
HowardGenes %<>% rename(cns_cell_type_taxon = cell_type_taxon)
HowardGenes %<>% left_join(cell_types %>% select(Gene, cell_type_taxon), by = c('mouseGene' = 'Gene'))


#HowardGenes %<>% group_by(Howard.Genes) %>% mutate(CellType = cell_types %>% filter(Gene == mouseGene) %>% select(cell_types) %>% as.character())

#for the genes that mapped to more than 1 mouse gene - remove them from cell type analysis - only one gene -- OR2B2
#HowardGenes %<>% mutate(CellType = ifelse(Howard.Genes == "OR2B2", "1+ mouse homolog - removed", CellType), mouseGene = ifelse(Howard.Genes == "OR2B2","1+ mouse homolog - removed", mouseGene)) %>% distinct() %>% select(-Updated_Gene_Names)

#remove the mouseGene for Olfr1360 and human gene OR2B2
HowardGenes %<>% filter(!(Howard.Genes == "OR2B2" &mouseGene == "Olfr1360"))

#Map brain structures to ambiguous brain regions 
brain_structures <- read_csv(here("ProcessedData", "AllenEtAl", "full_brain_region_hierarchy.csv"))
# brain_structures<-brain_structures[!duplicated(brain_structures$brain_region),] #remove duplicates 

HowardGenes %<>% mutate(BrainRegion = gsub(",","",BrainRegion))
HowardGenes %<>% mutate(BrainRegion = gsub("\\s+", " ", BrainRegion))
HowardGenes %<>% left_join(brain_structures, by = c('BrainRegion' = 'brain_region'))

brain_slim <- read_csv(here("ProcessedData", "AllenEtAl", "brain_regions_slim.csv"))
HowardGenes %<>% left_join(brain_slim)
HowardGenes %<>% mutate(location = ifelse(is.na(location), region_location, location)) %>% rename(slim_region_location = location)


HowardGenes %>% write_csv(path = here("ProcessedData", "HowardEtAl", "HowardRegionsPolygenicCellTypes_four.csv"))
#HowardGenes %>% write_csv(path = here("ProcessedData", "HowardEtAl", "HowardRegionsPolygenicCellTypes.csv"))

