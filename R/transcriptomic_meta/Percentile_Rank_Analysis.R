library(magrittr)
library(dplyr)


#Function that calculates the proportion of genes within the genome that has a smaller meta p_value than the current gene 
#Which gene had the smallest meta p-value compared to the entire genome
getRank <- function(dataset){
  print("Running rank calculations.. ")
  #get the number of genes
  num_genes <- (dataset %>% select(gene_symbol) %>% nrow())
  dataset %<>% mutate(meta_Up = rank(meta_higher_in_MDD_pvalue)/num_genes)
  dataset %<>% mutate(meta_Down = rank(meta_lower_in_MDD_pvalue)/num_genes)
  
  dataset %<>% rowwise() %>% mutate(empirical_meta_p = if_else(2 * min(meta_Up, meta_Down) > 1, 1, 2 * min(meta_Up, meta_Down))) %>% ungroup()
  return(dataset)
}

