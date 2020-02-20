library(magrittr)
library(dplyr)


#Function that calculates the proportion of genes within the genome that has a smaller meta p_value than the current gene 
#Which gene had the smallest meta p-value compared to the entire genome
getRank <- function(dataset){
  print("Running rank calculations.. ")
  #get the number of genes
  num_genes <- (dataset %>% select(gene_symbol) %>% nrow())
  dataset %<>% mutate(meta_Up = NA, meta_Down = NA)
  #filter for all genes that have higher expression levels
  higher_expression <- dataset %>% filter(meta_direction == +1)
  #filter all genes that have lower expression levels
  lower_expression <- dataset %>% filter(meta_direction == -1) 
  for (i in 1:num_genes) {
    if(i == ceiling(num_genes/2)){
      print("50% complete")
    }
    #get current genes calculated meta_p value
    gene_meta <- dataset$meta_p[i]
    #get current genes calculated meta direction
    direction_val <- dataset$meta_direction[i]
    #compare the current genes direction with the rest of the genome 
    if(direction_val == +1) {      
      #this current gene will be in the down regulation list, offset it
      dataset$meta_Up[i] <- higher_expression %>% filter(meta_p <= gene_meta) %>% select(gene_symbol) %>% nrow() / num_genes
      dataset$meta_Down[i] <- (lower_expression %>% filter(meta_p <= gene_meta) %>% select(gene_symbol) %>% nrow() + 1) / (num_genes + 1)
    
     }else if (direction_val == -1) {
      dataset$meta_Up[i] <- (higher_expression %>% filter(meta_p <= gene_meta) %>% select(gene_symbol) %>% nrow() + 1) / (num_genes + 1)
      dataset$meta_Down[i] <- lower_expression %>% filter(meta_p <= gene_meta) %>% select(gene_symbol) %>% nrow() / num_genes
    }
    
  }
  print("Complete!")
  dataset %<>% rowwise() %>% mutate(empirical_meta_p = if_else(2 * min(meta_Up, meta_Down) > 1, 1, 2 * min(meta_Up, meta_Down))) %>% ungroup()
  return(dataset)
}


