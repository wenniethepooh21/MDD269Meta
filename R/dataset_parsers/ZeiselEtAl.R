library(loomR)
library(dplyr)
library(homologene)

#-------- read cell-type matrix data 
loom_file <- connect(filename = here("Raw_Data/ZeiselEtAl/l6_r4.agg.loom"))
cell_types <- as_tibble(t(as.matrix(loom_file[["matrix"]][,])))
colnames(cell_types) <- loom_file$col.attrs$TaxonomyRank4[]
cell_types %<>% mutate(mouse_gene = loom_file$row.attrs$Gene[]) %>% select(mouse_gene, everything())
#deal with duplicate genes - avg the expression
full_cell_type <- cell_types %>% group_by(mouse_gene) %>% mutate_each(list(mean)) %>% distinct() %>% ungroup()
# identify the maximum cell-type taxon that expresses each mouse gene - read in file for get_max_expression function
source(here("R/transcriptomic_meta/Max_Expression_Matrix.R"))
max_cell_types <- get_max_expression(full_cell_type,"cell_type_taxon","mouse_gene")
#for those genes that had a max expression level of 0, updated the cell_type description
max_cell_types %<>% mutate(cell_type_taxon = if_else(expression_levels == 0,"Gene detected; No expression measured", cell_type_taxon))
# max_cell_types %>% write_csv(here("Processed_Data/ZeiselEtAl/max_cell_type_expression.csv"))

#---------
#max cell type with PNS neurons taxon removed 
taxonomy_group_two <- tibble(taxon_2 = loom_file$col.attrs$TaxonomyRank2[])
taxonomy_group_four <- tibble(taxon_4 = loom_file$col.attrs$TaxonomyRank4[])

taxonomy_assignment <- cbind(taxonomy_group_two, taxonomy_group_four)
#remove PNS cell-type taxons
cns_cell_type_taxa<- taxonomy_assignment %>% filter(taxon_2 != "PNS neurons") %>% select(taxon_4) %>% pull()
cns_cell_types <- cell_types %>% select(mouse_gene, cns_cell_type_taxa)
#deal with duplicate genes - avg the expression
cns_cell_types %<>% group_by(mouse_gene) %>% mutate_each(list(mean))  %>% distinct() %>% ungroup()
#re-run max cell type 
cns_max_cell_types <- get_max_expression(cns_cell_types,"cell_type_taxon","mouse_gene")
#for those genes that had a max expression level of 0, updated the cell_type description
cns_max_cell_types %<>% mutate(cell_type_taxon = if_else(expression_levels == 0,"Gene detected; No expression measured", cell_type_taxon))

# cns_max_cell_types %>% write_csv(path = here("ProcessedData", "ZeiselEtAl", "max_cns_cell_type_expression.csv"))


####################################################################
## Find the cell-type taxon that maximally expresses the 269 genes ####
####################################################################
howard<- read_csv(here("Processed_Data/HowardEtAl/fullHowardTable.csv")) 
howard_genes <- howard %>% select(gene_symbol) %>% pull()
howard %<>% select(gene_symbol, Updated_Gene_Names)

#Howard genes must be converted to mouse 
#gene_symbol and Updated_Gene_Names mapped to the same mouse genes, no need for double comparison
howard_mouse_genes <- human2mouse(howard_genes) %>% as_tibble()
howard_mouse_genes %<>% mutate(gene_symbol_upper = toupper(humanGene))
#Find which human gene mapped to more than 1 mouse homolog
dub_mouse_gene <- howard_mouse_genes %>% group_by(humanGene) %>% filter(n() > 1)
#check if OR2B2 is in howard_genes
howard %>% filter(gene_symbol == "OR2B2")
#check to see what cell-type was assigned to these two mouse genes
dub_mouse_gene %<>% ungroup() %>% select(mouseGene) %>% pull()
max_cell_types %>% filter(mouse_gene %in% dub_mouse_gene)
# assign OR2B2 to Olfr1359
howard_mouse_genes %<>% mutate(mouseGene = if_else(gene_symbol_upper == "OR2B2", "Olfr1359", mouseGene))

Howard_Table <- howard %>% mutate(gene_symbol_upper = toupper(gene_symbol))
Howard_Table %<>% left_join(howard_mouse_genes %>% select(gene_symbol_upper, mouseGene), by = c('gene_symbol_upper' = 'gene_symbol_upper'))
Howard_Table %<>% left_join(cns_max_cell_types %>% select(mouse_gene,cell_type_taxon), by = c('mouseGene' = 'mouse_gene'))
Howard_Table %<>% rename(cns_cell_type_taxon = cell_type_taxon)
Howard_Table %<>% left_join(max_cell_types %>% select(mouse_gene, cell_type_taxon), by = c('mouseGene' = 'mouse_gene'))
Howard_Table %<>% distinct() %>% select(-gene_symbol_upper)
#Check if changes were applies correctly
Howard_Table %>% filter(gene_symbol == "OR2B2") #should map to one cell-type taxon
Howard_Table %>% write_csv(here("Processed_Data/ZeiselEtAl/HowardCellTypes.csv"))
