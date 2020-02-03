library(homologene)
library(here)
library(magrittr)
library(dplyr)

#This script matches the howard 269 genes to the max brain region and cell type taxon that expressed it
source(here("R", "loading datasets", "fileManipulation.R"))
source(here("R", "loading datasets", "geneDirection.R"))
source(here("R", "loading datasets","ZeiselPolygenic.R"))


#-------------------------------------------------------------- Brain expression data 
four_donors_reannotated <- read_csv(here("Processed_Data/AllenEtAl/6_brains_reannotated_aggregated_4_donors.csv")) # file provided by Derek Howard
# identify the maximum cell-type taxon that expresses each human gene
# function is in this file - source it
source(here("R/transcriptomic_meta/Max_Expression_Matrix.R"))
max_four_donors_reannotated <- get_max_expression(four_donors_reannotated, "structure_name", "gene_symbol")
max_four_donors_reannotated %>% write_csv(path = here("Processed_Data/AllenEtAl/max_expression_in_four_donors.csv"))

howard<- read_csv(here("Processed_Data/HowardEtAl/fullHowardTable.csv")) 
howard_genes <- howard %>% select(gene_symbol) %>% pull()
#filter the brain expression data for the 269
brain_expression <- max_four_donors_reannotated %>% filter(gene_symbol %in% howard_genes)
brain_expression %<>% gather("brain_region", "expression", 2:ncol(brain_expression))

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


#-------- read cell-type matrix data 
loom_file <- connect(filename = here("Raw_Data/ZeiselEtAl/l6_r4.agg.loom"))
cell_types <- as_tibble(t(as.matrix(loom_file[["matrix"]][,])))
colnames(cell_types) <- loom_file$col.attrs$TaxonomyRank4[]
cell_types %<>% mutate(mouse_gene = loom_file$row.attrs$mouse_gene[]) %>% select(mouse_gene, everything())
#deal with duplicate genes - avg the expression
full_cell_type <- cell_types %<>% group_by(mouse_gene) %>% mutate_each(list(mean)) %>% distinct() %>% ungroup()
# identify the maximum cell-type taxon that expresses each mouse gene
max_cell_types <- get_max_expression(full_cell_type,"cell_type_taxon","mouse_gene")
#for those genes that had a max expression level of 0, updated the cell_type description
max_cell_types %<>% mutate(cell_type_taxon = if_else(expression_levels == 0,"Gene detected; No expression measured", cell_type_taxon))
max_cell_types %>% write_csv(here("Processed_Data/ZeiselEtAl/max_cell_type_expression.csv"))

#---------
#max cell type with PNS neurons taxon removed 
taxonomy_group_two <- loom_file$col.attrs$TaxonomyRank2[] %>% as_tibble()
taxonomy_group_two %<>% rename(Taxon_2 = value)
taxonomy_group_four <- loom_file$col.attrs$TaxonomyRank4[] %>% as_tibble()
taxonomy_group_four %<>% rename(Taxon_4 = value)

Taxonomy_assignment <- cbind(taxonomy_group_two, taxonomy_group_four)
cns_colnames<- Taxonomy_assignment %>% filter(Taxon_2 != "PNS neurons") %>% select(Taxon_4) %>% pull()

exp_matrix_slim <- cell_types %>% select(mouse_gene, cns_colnames)
# exp_matrix_slim %>% write_csv(here("ProcessedData", "ZeiselEtAl", "cns_cell_types_zeisel.R"))

#deal with duplicate genes - avg the expression
exp_matrix_slim %<>% group_by(mouse_gene) %>% mutate_each(list(mean))  %>% distinct() %>% ungroup()

cns_max_cell_types <- get_max_expression(exp_matrix_slim, "cell_type_taxon","mouse_gene")

#for those genes that had a max expression level of 0, updated the cell_type description
cns_max_cell_types %<>% mutate(cell_type_taxon = if_else(expression_levels == 0,"Gene detected; No expression measured", cell_type_taxon))
cns_max_cell_types %>% write_csv(path = here("ProcessedData", "ZeiselEtAl", "max_cns_cell_type_expression.csv"))


####################################################################
## Find the brain region that maximally expresses the 269 genes ####
####################################################################
howard<- read_csv(here("Processed_Data/HowardEtAl/fullHowardTable.csv")) 
howard_genes <- howard %>% select(gene_symbol) %>% pull()
howard %<>% select(gene_symbol, Updated_Gene_Names)

#Howard genes must be converted to mouse 
#gene_symbol and Updated_Gene_Names mapped to the same mouse genes, no need for double comparison
howard_mouse_genes <- human2mouse(howard$gene_symbol) %>% as_tibble()
howard_mouse_genes %<>% mutate(humanGene = toupper(humanGene))
Howard_Table <- howard %>% mutate(gene_symbol_upper = toupper(gene_symbol))
Howard_Table %<>% left_join(howard_mouse_genes %>% select(humanGene, mouseGene), by = c('gene_symbol_upper' = 'humanGene'))
Howard_Table %<>% left_join(cns_max_cell_types %>% select(mouse_gene,cell_type_taxon), by = c('mouseGene' = 'mouse_gene'))
Howard_Table %<>% rename(cns_cell_type_taxon = cell_type_taxon)
Howard_Table %<>% left_join(max_cell_types %>% select(mouse_gene, cell_type_taxon), by = c('mouseGene' = 'mouse_gene'))


#HowardGenes %>% group_by(gene_symbol) %>% mutate(cell_type = max_cell_types %>% filter(mouse_gene == mouseGene) %>% select(max_cell_types) %>% as.character())

#for the genes that mapped to more than 1 mouse gene - remove them from cell type analysis - only one gene -- OR2B2
HowardGenes %<>% mutate(cell_type = if_else(gene_symbol == "OR2B2", "1+ mouse homolog - removed", cell_type), mouseGene = if_else(gene_symbol == "OR2B2","1+ mouse homolog - removed", mouseGene)) %>% distinct() %>% select(-Updated_Gene_Names)

#remove the mouseGene for Olfr1360 and human gene OR2B2
Howard_Table %<>% filter(gene_symbol != "OR2B2") %>% filter(mouseGene != "Olfr1360")

#Map brain structures to ambiguous brain regions 
brain_structures <- read_csv(here("Processed_Data/AllenEtAl/full_brain_region_hierarchy.csv"))
# brain_structures<-brain_structures[!duplicated(brain_structures$brain_region),] #remove duplicates 

Howard_Table %<>% mutate(brain_region = gsub(",","",brain_region))
Howard_Table %<>% mutate(brain_region = gsub("\\s+", " ", brain_region))
Howard_Table %<>% left_join(brain_structures, by = c('brain_region' = 'brain_region'))

brain_slim <- read_csv(here("Processed_Data/AllenEtAl/brain_regions_slim.csv"))
Howard_Table %<>% left_join(brain_slim)
Howard_Table %<>% mutate(location = if_else(is.na(location), region_location, location)) %>% rename(slim_region_location = location)

Howard_Table %>% write_csv(path = here("Processed_Data/HowardEtAl/HowardRegionsPolygenicCellTypes_four.csv"))



