library(loomR)
library(dplyr)
library(homologene)
library(here)

#-------- read cell-type matrix data 
loom_file <- connect(filename = here("Raw_Data/ZeiselEtAl/l6_r4.agg.loom"))
cell_types <- as_tibble(t(as.matrix(loom_file[["matrix"]][,])))
colnames(cell_types) <- loom_file$col.attrs$TaxonomyRank4[]
cell_types %<>% mutate(mouse_gene = loom_file$row.attrs$Gene[]) %>% dplyr::select(mouse_gene, everything())
#deal with duplicate genes - avg the expression
full_cell_type <- cell_types %>% group_by(mouse_gene) %>% mutate_each(list(mean)) %>% distinct() %>% ungroup()

# identify the maximum cell-type taxon that expresses each mouse gene - read in file for get_max_expression function
source(here("R/transcriptomic_meta/Max_Expression_Matrix.R"))
# max_cell_types <- get_max_expression(full_cell_type,"cell_type_taxon","mouse_gene")
#for those genes that had a max expression level of 0, updated the cell_type description
# max_cell_types %<>% mutate(cell_type_taxon = if_else(expression_levels == 0,"Gene detected; No expression measured", cell_type_taxon))

cell_types_z <- get_z_score(full_cell_type)
cell_types_z %<>% rename(cell_type_taxon_zscore = zscore)
cell_types_z %>% write_csv(here("Processed_Data/ZeiselEtAl/full_cell_taxon_zscore.csv"))
max_cell_types <- get_z_score_max(cell_types_z)
#for those genes that had a max expression level of 0, updated the cell_type description
max_cell_types %<>% mutate(cell_type_taxon = if_else(expression_levels == 0,"Gene detected; No expression measured", cell_type_taxon))

#---------
#max cell type with PNS neurons taxon removed 
taxonomy_group_two <- tibble(taxon_2 = loom_file$col.attrs$TaxonomyRank2[])
taxonomy_group_four <- tibble(taxon_4 = loom_file$col.attrs$TaxonomyRank4[])

taxonomy_assignment <- cbind(taxonomy_group_two, taxonomy_group_four)
#remove PNS cell-type taxons
cns_cell_type_taxa<- taxonomy_assignment %>% filter(taxon_2 != "PNS neurons") %>% dplyr::select(taxon_4) %>% pull()
cns_cell_types <- cell_types %>% dplyr::select(mouse_gene, cns_cell_type_taxa)
#deal with duplicate genes - avg the expression
cns_cell_types %<>% group_by(mouse_gene) %>% mutate_each(list(mean))  %>% distinct() %>% ungroup()
#re-run max cell type 
# cns_max_cell_types <- get_max_expression(cns_cell_types,"cns_cell_type_taxon","mouse_gene")

cns_max_cell_types <- get_z_score(cns_cell_types)
cns_max_cell_types %<>% rename(cns_cell_type_taxon = cell_type_taxon, cns_cell_type_taxon_zscore = zscore)
#for those genes that had a max expression level of 0, updated the cell_type description
cns_max_cell_types %<>% mutate(cns_cell_type_taxon = if_else(expression_levels == 0,"Gene detected; No expression measured", cns_cell_type_taxon))
cns_max_cell_types %<>% get_z_score_max()

####################################################################
## Find the cell-type taxon that maximally expresses the 269 genes ####
####################################################################
howard<- read_csv(here("Processed_Data/HowardEtAl/fullHowardTable.csv")) 
howard_genes <- howard %>% dplyr::select(gene_symbol) %>% pull()
howard %<>% dplyr::select(gene_symbol, Updated_Gene_Names)


#Howard genes must be converted to mouse 
#gene_symbol and Updated_Gene_Names mapped to the same mouse genes, no need for double comparison
howard_mouse_genes <- human2mouse(howard_genes) %>% as_tibble()
howard_mouse_genes %<>% mutate(gene_symbol_upper = toupper(humanGene))
#Find which human gene mapped to more than 1 mouse homolog
dub_mouse_gene <- howard_mouse_genes %>% group_by(humanGene) %>% filter(n() > 1)
#check if OR2B2 is in howard_genes
howard %>% filter(gene_symbol == "OR2B2") #should map to one cell-type taxon
#check to see what cell-type was assigned to these two mouse genes
dub_mouse_gene %<>% ungroup() %>% dplyr::select(mouseGene) %>% pull()
max_cell_types %>% filter(mouse_gene %in% dub_mouse_gene)
# assign OR2B2 to Olfr1359
howard_mouse_genes %<>% mutate(mouseGene = if_else(gene_symbol_upper == "OR2B2", "Olfr1359", mouseGene))


Howard_Table <- howard %>% mutate(gene_symbol_upper = toupper(gene_symbol))
Howard_Table %<>% left_join(howard_mouse_genes %>% dplyr::select(gene_symbol_upper, mouseGene), by = c('gene_symbol_upper' = 'gene_symbol_upper'))
Howard_Table %<>% left_join(cns_max_cell_types %>% dplyr::select(mouse_gene,cns_cell_type_taxon,cns_cell_type_taxon_zscore), by = c('mouseGene' = 'mouse_gene')) %>% distinct()
Howard_Table %<>% left_join(max_cell_types %>% dplyr::select(mouse_gene, cell_type_taxon, cell_type_taxon_zscore), by = c('mouseGene' = 'mouse_gene'))
Howard_Table %<>% dplyr::select(-gene_symbol_upper)

Howard_Table %>% filter(is.na(cell_type_taxon_zscore))
Howard_Table %>% write_csv(here("Processed_Data/ZeiselEtAl/HowardCellTypes_zscore.csv"))


# #############################calculate the probability for the cell types
cell_pop <- nrow(max_cell_types)
cell_count <- max_cell_types %>% group_by(cell_type_taxon) %>% summarize(genome_cell_count = n())

howard_cell <- Howard_Table
howard_cell %<>% mutate(cell_type_taxon = ifelse(cell_type_taxon == "character(0)", NA, cell_type_taxon))
howard_count <- howard_cell %>% group_by(cell_type_taxon) %>% summarise(sample_cell_count = n())

full_cell_count <- cell_count %>% left_join(howard_count, by=c('cell_type_taxon' = 'cell_type_taxon') )
# perform hypergeometric test - read in file for hyper_test function
source(here("R/transcriptomic_meta/hyper_test.R"))
#subtract one for the gene that wasn't assigned a cell type but was observed
cell_expected_probs  <- full_cell_count %>% rowwise() %>% mutate(hypergeometric_p = hyper_test(sample_cell_count, genome_cell_count, cell_pop, colSums(na.omit(howard_count)[,2]) - 1))
#corrected by the number of possible cell types to choose from subtract one because of the 'gene not expressed option'
cell_expected_probs %<>% mutate(corrected_hypergeometric_p = p.adjust(hypergeometric_p, method = "bonferroni", n= (nrow(cell_count) - 1)))
cell_expected_probs %<>% mutate(corrected_hypergeometric_p = ifelse(cell_type_taxon == "Gene detected; No expression measured", NA, corrected_hypergeometric_p))
cell_expected_probs %<>% arrange(hypergeometric_p)

cell_expected_probs$hypergeometric_p <- signif(as.numeric(cell_expected_probs$hypergeometric_p),digits=3)
cell_expected_probs$corrected_hypergeometric_p <- signif(as.numeric(cell_expected_probs$corrected_hypergeometric_p),digits=3)

write_csv(cell_expected_probs, path = here('Results', 'supplementary_tables', 'hypergeometric_cell_type_taxons_cell_expected_probs.csv'))


# #################################
cns_cell_pop <- nrow(cns_max_cell_types)
cns_cell_count <- cns_max_cell_types %>% group_by(cns_cell_type_taxon) %>% summarize(genome_cell_count = n())

cns_howard_cell <- Howard_Table
cns_howard_cell %<>% mutate(cns_cell_type_taxon = ifelse(cns_cell_type_taxon == "character(0)", NA, cns_cell_type_taxon))
cns_howard_count <- cns_howard_cell %>% group_by(cns_cell_type_taxon) %>% summarise(sample_cell_count = n())

#full_count <- howard_count %>% left_join(cell_count, by=c('cell_types' = 'cell_types') )
cns_full_count <- cns_cell_count %>% left_join(cns_howard_count, by=c('cns_cell_type_taxon' = 'cns_cell_type_taxon') )

#subtract one for the gene that wasn't assigned a cell type but was observed
cns_cell_expected_probs  <- cns_full_count %>% rowwise() %>% mutate(hypergeometric_p = hyper_test(sample_cell_count, genome_cell_count, cell_pop, colSums(na.omit(howard_count)[,2]) - 1))
#corrected by the number of possible cell types to choose from subtract one because of the 'gene not expressed option'
cns_cell_expected_probs %<>% mutate(corrected_hypergeometric_p = p.adjust(hypergeometric_p, method = "bonferroni", n = nrow(cell_count)-1))
cns_cell_expected_probs %<>% mutate(corrected_hypergeometric_p = ifelse(cns_cell_type_taxon == "Gene detected; No expression measured", NA, corrected_hypergeometric_p))
cns_cell_expected_probs %<>% arrange(hypergeometric_p)

cns_cell_expected_probs$hypergeometric_p <- signif(as.numeric(cns_cell_expected_probs$hypergeometric_p),digits=3)
cns_cell_expected_probs$corrected_hypergeometric_p <- signif(as.numeric(cns_cell_expected_probs$corrected_hypergeometric_p),digits=3)
write_csv(cns_cell_expected_probs, path = here('Results', 'supplementary_tables', 'hypergeometric_cns_cell_type_taxons_cell_expected_probs.csv'))

