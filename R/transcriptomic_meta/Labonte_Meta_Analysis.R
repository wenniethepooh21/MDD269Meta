library(magrittr)
library(dplyr)
library(tidyr)
library(metap)

#Function that calculates the meta p values using Fisher's method
#args: table_name == Labonte dataset, num_p_val_val == number of p-values that we expect to have, num_regions == 
LabonteMetaAnalysis <- function(full_results, num_p_val_val, num_regions) {

  #tempary drop of missing log fold changes - due to excel date issues that will fixed later or due to missing genes between male/female
  full_results %<>% filter(!is.na(logFC))
  
  #p-values are two tailed for negative and positive expression levels 
  #create six one-sided pvalues per gene based on the direction of expression (get a p-value for both directions separately to see which is more significant)
  full_results %<>% mutate(higher_in_MDD_pvalue = two2one(pvalue, invert= 1 == sign(-1*logFC)))#invert p-value if expression is lower (1-p_value) otherwise (pvalue/2)
  full_results %<>% mutate(lower_in_MDD_pvalue = two2one(pvalue, invert=1 == sign(1*logFC)))#invert p-value if expression is higher (1-p_value) otherwise (p_value/2)
  full_results %<>% rename(P.Value = pvalue)
  
  #Find the genes that have p-values in all brain regions
  genes_with_pvalues <- unique(full_results %>% group_by(gene_symbol) %>% summarize(n=n(), unique_regions = length(unique(brain_region))) %>% filter(n==num_p_val_val, unique_regions == num_regions) %>% .$gene_symbol)
  #perform meta-analysis 
  #Run Fisher's method grouping genes across brain regions (if there's more than one brain region)
  #combining the information in the p-values from different statistical tests to form a single overall test
  full_results %<>% filter(gene_symbol %in% genes_with_pvalues)
  summary_results <- full_results %>% group_by(gene_symbol) %>% summarize(min_p_across_regions = min(P.Value), 
                                                                          meta_higher_in_MDD_pvalue = sumlog(c(higher_in_MDD_pvalue))$p,
                                                                          meta_lower_in_MDD_pvalue = sumlog(c(lower_in_MDD_pvalue))$p)
  
  #meta p-value is the direction with the more significant p-value, multiply by 2 to change back to 2 sided p-value
  summary_results %<>% rowwise() %>% mutate(meta_direction = if_else(meta_higher_in_MDD_pvalue < meta_lower_in_MDD_pvalue, 1, -1), meta_p = 2 * min(meta_higher_in_MDD_pvalue, meta_lower_in_MDD_pvalue))
  
  #append different brain regions together for column header
  full_results %<>% mutate(combined_region_sex = paste(sex, brain_region, sep="_"))
  combined_order <- sort(unique(full_results$combined_region_sex))
  combined_order_string = paste(combined_order, collapse="_")
  
  #summarize the direction of expression with '+' and '-' 
  directions <- full_results %>% select(gene_symbol, combined_region_sex, logFC)
  directions %<>% spread(combined_region_sex, logFC)
  directions %<>% mutate_at(combined_order, list(~ if_else(. > 0, "+", "-")))
  directions %<>% unite(col = !!combined_order_string, combined_order, sep="", remove = TRUE)
  
  Labonte_summary_results <- inner_join(directions, summary_results)
  Labonte_summary_results %<>% arrange(meta_p)
  
  return(Labonte_summary_results)
}


