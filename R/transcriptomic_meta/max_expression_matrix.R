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

