library(homologene)
library(readr)
library(dplyr)
library(tidyr)
library(reshape2)
library(magrittr)
library(here)
library(googledrive)
library(googlesheets4)


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

# identify the maximum cell-type taxon that expresses each mouse gene - read in file for get_max_expression function
source(here("R/transcriptomic_meta/Max_Expression_Matrix.R"))
max_cell_types <- get_max_expression(table,"cluster_id","Gene")
#for those genes that had a max expression level of 0, updated the cell_type description
max_cell_types %<>% mutate(cluster_id = if_else(expression_levels == 0,"Gene detected; No expression measured", cluster_id))

#change cluster id to more descriptive name
cell_type_info <- rbind(colnames(table[2:266]), descriptions)
cell_type_info <- t(cell_type_info)
colnames(cell_type_info) <- c('cluster_id', 'cluster_id')
cell_type_info <- as_tibble(cell_type_info)
cell_type_info %<>% rename(cluster_description = V2)

# max_cell_types %<>% left_join(cell_type_info, by = c('cluster_id' = 'cluster_id'))
# max_cell_types %<>% mutate(cluster_description = if_else(cluster_id == "Gene detected; No expression measured","Gene detected; No expression measured", cluster_description))
max_cell_types %<>% rename(mouse_gene = Gene)

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
Howard_Table %<>% left_join(max_cell_types %>% select(mouse_gene, cluster_id), by = c('mouseGene' = 'mouse_gene'))
Howard_Table %<>% distinct() %>% select(-gene_symbol_upper)
#Check if changes were applies correctly
Howard_Table %>% filter(gene_symbol == "OR2B2") #should map to one cell-type taxon
Howard_Table %>% left_join(cell_type_info, by = c('cluster_id' = 'cluster_id')) %>% write_csv(here("Processed_Data/ZeiselEtAl/HowardCellTypes_265.csv"))


# #############################calculate the probability for the cell types
cell_pop <- nrow(max_cell_types)
cell_count <- max_cell_types %>% group_by(cluster_id) %>% summarize(genome_cell_count = n())

howard_cell <- Howard_Table
howard_cell %<>% mutate(cluster_id = ifelse(cluster_id == "character(0)", NA, cluster_id))
howard_count <- howard_cell %>% group_by(cluster_id) %>% summarise(sample_cell_count = n())

full_cell_count <- cell_count %>% left_join(howard_count, by=c('cluster_id' = 'cluster_id') )
# perform hypergeometric test - read in file for hyper_test function
source(here("R/transcriptomic_meta/hyper_test.R"))
#subtract one for the gene that wasn't assigned a cell type but was observed
cell_expected_probs  <- full_cell_count %>% rowwise() %>% mutate(hypergeometric_p = hyper_test(sample_cell_count, genome_cell_count, cell_pop, colSums(na.omit(howard_count)[,2]) - 1))
#corrected by the number of possible cell types to choose from subtract one because of the 'gene not expressed option'
cell_expected_probs %<>% mutate(corrected_hypergeometric_p = p.adjust(hypergeometric_p, method = "bonferroni", n= (nrow(cell_count) - 1)))
cell_expected_probs %<>% mutate(corrected_hypergeometric_p = ifelse(cluster_id == "Gene detected; No expression measured", NA, corrected_hypergeometric_p))
cell_expected_probs %<>% arrange(hypergeometric_p)

cell_expected_probs$hypergeometric_p <- signif(as.numeric(cell_expected_probs$hypergeometric_p),digits=3)
cell_expected_probs$corrected_hypergeometric_p <- signif(as.numeric(cell_expected_probs$corrected_hypergeometric_p),digits=3)
cell_expected_probs<-cell_type_info %>% left_join(cell_expected_probs, by = c('cluster_id' = 'cluster_id'))  %>% arrange(hypergeometric_p)
#upload to google drive
sheets_auth(token = drive_token())

cells <- drive_get("~/Thesis/Manuscript/Supplement_Tables/cell_hyper_expected_265")
if(nrow(cells) != 0) {
  drive_rm(cells)
}
#create the google worksheet
cells <- sheets_create("cell_hyper_expected_265",sheets = c('hypergeometric_cell_type_taxons_265'))
sheets_write(cell_expected_probs, cells,  sheet = "hypergeometric_cell_type_taxons_265")

drive_mv(file = "cell_hyper_expected_265", path = "~/Thesis/Manuscript/Supplement_Tables/")  # move Sheets file

