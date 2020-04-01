# visulalize the 39 and 265 classification for each gene
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(magrittr)
library(homologene)

table <- read_tsv(here('./Raw_Data/ZeiselEtAl/l5_all.agg.tab'), col_names = FALSE)
descriptions <- table[29, 9:273]

table %<>% filter(!is.na(X1) | X8 == 'ClusterName')
colnames(table) <- c(unlist(table[2,1:7]), unlist(table[1, 8:ncol(table)]))
table <- table[-c(1,2),]

# columns to drop
table <- table[!is.na(names(table))]
dropcols <- c('Accession', '_LogCV', '_LogMean', '_Selected', '_Total', '_Valid', 'ClusterName')
table %<>% select(-one_of(dropcols))

#to change all but the 'Gene' column name into numeric type
table[,2:length(colnames(table))] %<>% lapply(function(x) as.numeric(as.character(x)))
#deal with duplicate genes - avg the expression
table %<>% group_by(Gene) %>% mutate_each(list(mean)) %>% distinct() %>% ungroup()

long_table <- gather(table, cluster_id, expression_level, 2:ncol(table))

chol_monamin_neurons_id <- c('DECHO1', 'HBADR', 'HBCHO3', 'HBCHO4', 'HBNOR', 'HBSER1', 'HBSER2', 'HBSER3', 'HBSER4', 'HBSER5', 'HYPEP6', 'HYPEP7', 'MBDOP1', 'MBDOP2', 'MEGLU14', 'TECHO')
significant_cluster <- tibble(cluster_id = chol_monamin_neurons_id, cell_type_taxon = 'cholinergic and monoaminergic neurons')
enteric_neurons_id <- c('ENT1', 'ENT2', 'ENT3', 'ENT4', 'ENT5', 'ENT6', 'ENT7', 'ENT8', 'ENT9')
significant_cluster %<>% rbind(tibble(cluster_id = enteric_neurons_id, cell_type_taxon = 'Enteric neurons'))

#both cell-type taxa
full_long_table <- long_table %>% filter(cluster_id %in% significant_cluster$cluster_id)
#cholinergic taxon
chol_long_table <- long_table %>% filter(cluster_id %in% chol_monamin_neurons_id)
#enteric taxon
ent_long_table <- long_table %>% filter(cluster_id %in% enteric_neurons_id)


full_slim_table <- spread(full_long_table,cluster_id, expression_level)
chol_slim_table <- spread(chol_long_table,cluster_id, expression_level)
ent_slim_table <- spread(ent_long_table,cluster_id, expression_level)

source(here("R/transcriptomic_meta/Max_Expression_Matrix.R"))
full_max_cell_types <- get_max_expression(full_slim_table,"cluster_id","Gene")
#for those genes that had a max expression level of 0, updated the cell_type description
full_max_cell_types %<>% mutate(cluster_id = if_else(expression_levels == 0,"Gene detected; No expression measured", cluster_id))

chol_max_cell_types <- get_max_expression(chol_slim_table,"cluster_id","Gene")
#for those genes that had a max expression level of 0, updated the cell_type description
chol_max_cell_types %<>% mutate(cluster_id = if_else(expression_levels == 0,"Gene detected; No expression measured", cluster_id))

ent_max_cell_types <- get_max_expression(ent_slim_table,"cluster_id","Gene")
#for those genes that had a max expression level of 0, updated the cell_type description
ent_max_cell_types %<>% mutate(cluster_id = if_else(expression_levels == 0,"Gene detected; No expression measured", cluster_id))

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
full_max_cell_types %>% filter(Gene %in% dub_mouse_gene)
chol_max_cell_types %>% filter(Gene %in% dub_mouse_gene)
ent_max_cell_types %>% filter(Gene %in% dub_mouse_gene)

# assign OR2B2 to Olfr1359 mouse homolog
howard_mouse_genes %<>% mutate(mouseGene = if_else(gene_symbol_upper == "OR2B2", "Olfr1359", mouseGene))
Howard_Table <- howard %>% mutate(gene_symbol_upper = toupper(gene_symbol))
Howard_Table %<>% left_join(howard_mouse_genes %>% select(gene_symbol_upper, mouseGene), by = c('gene_symbol_upper' = 'gene_symbol_upper'))
full_Howard_Table <- Howard_Table %>% left_join(full_max_cell_types %>% select(Gene, cluster_id), by = c('mouseGene' = 'Gene')) %>% distinct() %>% select(-gene_symbol_upper)
chol_Howard_Table <- Howard_Table %>% left_join(chol_max_cell_types %>% select(Gene, cluster_id), by = c('mouseGene' = 'Gene')) %>% distinct() %>% select(-gene_symbol_upper)
ent_Howard_Table <- Howard_Table %>% left_join(ent_max_cell_types %>% select(Gene, cluster_id), by = c('mouseGene' = 'Gene')) %>% distinct() %>% select(-gene_symbol_upper)

# #############################calculate the probability for the cell types

get_hyper <- function(max_table, howard_table) {
  
  cell_pop <- nrow(max_table)
  cell_count <- max_table %>% group_by(cluster_id) %>% summarize(genome_cell_count = n())
  
  howard_cell <- howard_table
  howard_cell %<>% mutate(cluster_id = ifelse(cluster_id == "character(0)", NA, cluster_id))
  howard_count <- howard_cell %>% group_by(cluster_id) %>% summarise(sample_cell_count = n())
  
  full_cell_count <- cell_count %>% left_join(howard_count, by=c('cluster_id' = 'cluster_id') )
  #
  # perform hypergeometric test - read in file for hyper_test function
  source(here("R/transcriptomic_meta/hyper_test.R"))
  #subtract one for the gene that wasn't assigned a cell type but was observed
  cell_expected_probs  <- full_cell_count %>% rowwise() %>% mutate(hypergeometric_p = hyper_test(sample_cell_count, genome_cell_count, cell_pop, colSums(na.omit(howard_count)[,2]) - 1))
  #corrected by the number of possible cell types to choose from subtract one because of the 'gene not expressed option'
  cell_expected_probs %<>% ungroup() %>% mutate(corrected_hypergeometric_p = p.adjust(hypergeometric_p, method = "bonferroni"))
  cell_expected_probs %<>% mutate(corrected_hypergeometric_p = ifelse(cluster_id == "Gene detected; No expression measured", NA, corrected_hypergeometric_p))
  cell_expected_probs %<>% arrange(hypergeometric_p)
  
  cell_expected_probs$hypergeometric_p <- signif(as.numeric(cell_expected_probs$hypergeometric_p),digits=3)
  cell_expected_probs$corrected_hypergeometric_p <- signif(as.numeric(cell_expected_probs$corrected_hypergeometric_p),digits=3)
  
  return(cell_expected_probs)
}
#Hyper test results 
full_hyper <- get_hyper(full_max_cell_types, full_Howard_Table)
chol_hyper <- get_hyper(chol_max_cell_types, chol_Howard_Table)
ent_hyper <- get_hyper(ent_max_cell_types, ent_Howard_Table)
