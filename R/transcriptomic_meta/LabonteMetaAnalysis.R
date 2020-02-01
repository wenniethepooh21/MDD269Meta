library(magrittr)
library(dplyr)
library(tidyverse)
library(devtools)
library(metap)

LabonteMeta <- function(tableName, numP, numRegions) {

#tempary drop of missing log fold changes - due to excel date issues that will fixed later or due to missing genes between male/female
tableName %<>% filter(!is.na(logFC))

#Meta-p analysis here
tableName %<>% mutate(higher_in_MDD_pvalue = two2one(pvalue, invert= 1 == sign(-1*logFC)))
tableName %<>% mutate(lower_in_MDD_pvalue = two2one(pvalue, invert=1 == sign(1*logFC)))
tableName %<>% rename(gene_symbol = symbol, P.Value = pvalue)

#maybe deal with these later
genes_with_pvalues <- unique(tableName %>% group_by(gene_symbol) %>% summarize(n=n(), unique_regions = length(unique(brain_region))) %>% filter(n==numP, unique_regions == numRegions) %>% .$gene_symbol)
tableName %<>% filter(gene_symbol %in% genes_with_pvalues)


Labonte_summary_results <- tableName %>% group_by(gene_symbol) %>% summarize(min_p_across_regions = min(P.Value), 
                                                                        meta_higher_in_MDD_pvalue = sumlog(c(higher_in_MDD_pvalue))$p,
                                                                        meta_lower_in_MDD_pvalue = sumlog(c(lower_in_MDD_pvalue))$p)
#convert from two one-sided meta pvalues to one
Labonte_summary_results %<>% rowwise() %>% mutate(meta_direction = if_else(meta_higher_in_MDD_pvalue < meta_lower_in_MDD_pvalue, 1, -1), meta_p = if_else((2 * min(meta_higher_in_MDD_pvalue, meta_lower_in_MDD_pvalue) )> 1,1, 2*min(meta_higher_in_MDD_pvalue,meta_lower_in_MDD_pvalue)))

tableName %<>% mutate(combined_region_sex = paste(sex, brain_region, sep="_"))

combined_order = sort(unique(tableName$combined_region_sex))
combined_order_string = paste(combined_order, collapse="_")

#add in individual directions for visualization
directions <- tableName %>% dplyr::select(gene_symbol, combined_region_sex, logFC)
#spread seems to do alphabetical sorting here too
directions %<>% spread(combined_region_sex, logFC)
directions %<>% mutate_at(combined_order, list(~ if_else(. > 0, "+", "-")))
directions %<>% unite(col = !!combined_order_string, combined_order, sep="", remove = TRUE)

Labonte_summary_results <- inner_join(directions, Labonte_summary_results)
Labonte_summary_results %<>% arrange(meta_p)

return(Labonte_summary_results)
}


getRank <- function(dataset, num_genes){
  dataset %<>% mutate(meta_Up = NA, meta_Down = NA)
  Labonte_Up <- dataset %>% filter(meta_direction == +1)
  Labonte_Down <- dataset %>% filter(meta_direction == -1)
  for (i in 1:num_genes) {
    gene_meta <- dataset$meta_p[i]
    direction_val <- dataset$meta_direction[i]
    
    #this current gene will be in the up regulated list, offset it
    if(direction_val == +1) {
      dataset$meta_Up[i] <- Labonte_Up %>% filter(meta_p <= gene_meta) %>% dplyr::select(gene_symbol) %>% nrow() / num_genes
      dataset$meta_Down[i] <- (Labonte_Down %>% filter(meta_p <= gene_meta) %>% dplyr::select(gene_symbol) %>% nrow() + 1) / (num_genes + 1)
      
      #this current gene will be in the down regulation list, offset it
    }else if (direction_val == -1) {
      dataset$meta_Up[i] <- (Labonte_Up %>% filter(meta_p <= gene_meta) %>% dplyr::select(gene_symbol) %>% nrow() + 1) / (num_genes + 1)
      dataset$meta_Down[i] <- Labonte_Down %>% filter(meta_p <= gene_meta) %>% dplyr::select(gene_symbol) %>% nrow() / num_genes
    }
  }
  dataset %<>% rowwise() %>% mutate(genome_percentile_rank = if_else(2 * min(meta_Up, meta_Down) > 1, 1, 2 * min(meta_Up, meta_Down))) 
  return(dataset)
}