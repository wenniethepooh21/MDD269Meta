library(dplyr)
library(tidyr)

#script to generate the max expression matrix for cell types and brain regions 


#Function that transforms the expression matrix data to find the cell types and brain regions with the highest gene level expression
get_max_expression <-function(full_expression_matrix, key.col, group.col) {
  #expand matrix to have the genes as rows instead of columns 
  exp_matrix_long <- gather_(full_expression_matrix, colnames(full_expression_matrix)[-1], key = key.col, value = "expression_levels")
  max_expression <- exp_matrix_long %>% group_by_(group.col) %>% dplyr::slice(which.max(expression_levels)) %>% ungroup()
  return(max_expression)
}

get_z_score <- function(full_expression_matrix) {
  data <- get_mean_sd(full_expression_matrix)
  expression <- gather_(full_expression_matrix, colnames(full_expression_matrix)[-1], key = "cell_type_taxon", value = "expression_levels")
  expression %<>% left_join(data)
  z_score_table <- expression%>% calc_z_score()
  return(z_score_table)
}

get_z_score_max <- function(full_zscore_matrix) {
  max_zscore_table <- full_zscore_matrix %>% group_by(mouse_gene) %>% dplyr::slice(which.max(expression_levels)) %>% ungroup()
  return(max_zscore_table)
}

get_mean_sd <- function(df) {
  gene_means <- df %>% mutate(avg = rowMeans(dplyr::select(., -mouse_gene))) %>% dplyr::select(mouse_gene,avg)
  gene_sd <- df %>% mutate(sd = apply(dplyr::select(.,-mouse_gene), 1, sd)) %>% dplyr::select(mouse_gene, sd)
  data <- inner_join(gene_means,gene_sd)
  return(data)
  
}

calc_z_score <- function(df) {
  df %<>% rowwise() %>% mutate(zscore = (expression_levels - avg)/sd)
  return(df)
}


