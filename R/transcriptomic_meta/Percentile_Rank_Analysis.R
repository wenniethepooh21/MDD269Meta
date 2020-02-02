library(magrittr)
library(dplyr)


#Function that calculates the proportion of genes within the genome that has a smaller meta p_value than the current gene 
#Which gene had the smallest meta p-value compared to the entire genome
getRank <- function(dataset){
  print("Running rank calculations.. ")
  #get the number of genes
  num_genes <- dataset %>% select(gene_symbol) %>% nrow()
  dataset %<>% mutate(meta_Up = NA, meta_Down = NA)
  higher_expression <- dataset %>% filter(meta_direction == +1)
  lower_expression <- dataset %>% filter(meta_direction == -1) 
  for (i in 1:num_genes) {
    if(i == num_genes/2){
      print("50% complete")
    }
    
    gene_meta <- dataset$meta_p[i]
    direction_val <- dataset$meta_direction[i]
    if(direction_val == +1) {
      dataset$meta_Up[i] <- higher_expression %>% filter(meta_p <= gene_meta) %>% select(gene_symbol) %>% nrow() / num_genes
      dataset$meta_Down[i] <- (lower_expression %>% filter(meta_p <= gene_meta) %>% select(gene_symbol) %>% nrow() + 1) / (num_genes + 1)
      
      #this current gene will be in the down regulation list, offset it
    }else if (direction_val == -1) {
      dataset$meta_Up[i] <- (higher_expression %>% filter(meta_p <= gene_meta) %>% select(gene_symbol) %>% nrow() + 1) / (num_genes + 1)
      dataset$meta_Down[i] <- lower_expression %>% filter(meta_p <= gene_meta) %>% select(gene_symbol) %>% nrow() / num_genes
    }
    
  }
  print("Complete!")
  dataset %<>% rowwise() %>% mutate(empirical_meta_p = if_else(2 * min(meta_Up, meta_Down) > 1, 1, 2 * min(meta_Up, meta_Down)))
  return(dataset)
}


