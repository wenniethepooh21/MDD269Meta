library(dplyr)
library(loomR)
library(here)
library(stringr)
library(readr)
library(tidyr)

#script to generate the max expression matrix for cell types and brain regions 


#Function that transforms the expression matrix data to find the cell types and brain regions with the highest gene level expression
get_max_expression <-function(full_expression_matrix, key.col, group.col) {
  #expand matrix to have the genes as rows instead of columns 
  exp_matrix_long <- gather_(full_expression_matrix, colnames(full_expression_matrix)[-1], key = key.col, value = "expression_levels")
  max_expression <- exp_matrix_long %>% group_by_(group.col) %>% slice(which.max(expression_levels)) %>% ungroup()
  return(max_expression)
}

#-------- Cell_Types expression matrix
loom_file <- connect(filename = here("data", "ZeiselEtAl", "l6_r4.agg.loom"))
cell_types <- as_tibble(t(as.matrix(loom_file[["matrix"]][,])))
colnames(cell_types) <- loom_file$col.attrs$TaxonomyRank4[]
cell_types %<>% mutate(Gene = loom_file$row.attrs$Gene[]) %>% select(Gene, everything())

#cell_types %>% write_csv(here("ProcessedData", "ZeiselEtAl", "39_types_zeisel.R"))

#deal with duplicate genes - avg the expression
full_cell_type <- cell_types %<>% group_by(Gene) %>% mutate_each(list(mean))  %>% distinct() %>% ungroup()

max_cell_types <- get_max_expression(full_cell_type, "cell_type_taxon","Gene")

#for those genes that had a max expression level of 0, updated the cell_type description
max_cell_types %<>% mutate(cell_type_taxon = ifelse(expression_levels == 0,"Gene detected; No expression measured", cell_type_taxon))
max_cell_types %>% write_csv(path = here("ProcessedData", "ZeiselEtAl", "max_cell_type_expression.csv"))

#---------
      #max cell type with PNS neurons taxon removed 
      
      taxonomy_group_two <- loom_file$col.attrs$TaxonomyRank2[] %>% as_tibble() 
      taxonomy_group_two %<>% rename(Taxon_2 = value)
      taxonomy_group_four <- loom_file$col.attrs$TaxonomyRank4[] %>% as_tibble()
      taxonomy_group_four %<>% rename(Taxon_4 = value)
      
      Taxonomy_assignment <- cbind(taxonomy_group_two, taxonomy_group_four)
      cns_colnames<- Taxonomy_assignment %>% filter(Taxon_2 != "PNS neurons") %>% select(Taxon_4) %>% pull()
      
      exp_matrix_slim <- cell_types %>% select(Gene, cns_colnames)
      # exp_matrix_slim %>% write_csv(here("ProcessedData", "ZeiselEtAl", "cns_cell_types_zeisel.R"))
      
      #deal with duplicate genes - avg the expression
      exp_matrix_slim %<>% group_by(Gene) %>% mutate_each(list(mean))  %>% distinct() %>% ungroup()
      
      cns_max_cell_types <- get_max_expression(exp_matrix_slim, "cell_type_taxon","Gene")
      
      #for those genes that had a max expression level of 0, updated the cell_type description
      cns_max_cell_types %<>% mutate(cell_type_taxon = ifelse(expression_levels == 0,"Gene detected; No expression measured", cell_type_taxon))
      cns_max_cell_types %>% write_csv(path = here("ProcessedData", "ZeiselEtAl", "max_cns_cell_type_expression.csv"))
      
      

#--------------------------------------------------------------
#---------- 6 donors re-annotated expression matrix
six_donors_reannotated <- read_csv(here("ProcessedData", "AllenEtAl", "6_donors_reannotated_aggregated.csv"))

four_donors_reannotated <- read_csv(here("ProcessedData", "AllenEtAl", "6_brains_reannotated_aggregated_4_donors.csv"))
max_four_donors_reannotated <- get_max_expression(four_donors_reannotated, "structure_name", "gene_symbol")
max_six_donors_reannotated <- get_max_expression(six_donors_reannotated, "gene_symbol", "gene_symbol")
max_four_donors_reannotated %>% write_csv(path = here("ProcessedData", "AllenEtAl", "max_expression_in_four_donors.csv"))
max_six_donors_reannotated %>%  write_csv(path = here("ProcessedData", "AllenEtAl", "max_expression_in_six_donors.csv"))


#--- assess expression levels of 269 genes in the sampled cortical and subcortical brain regions
howard_genes <- read_csv(here("ProcessedData", "HowardEtAl", "fullHowardTable.csv")) %>% select(gene_symbol) %>% pull()
gwas_expression <- four_donors_reannotated %>% filter(gene_symbol %in% howard_genes)
gwas_expression %<>% gather("brain_region", "expression", 2:ncol(gwas_expression))

#find the matching brain regions to transcriptomic studies


