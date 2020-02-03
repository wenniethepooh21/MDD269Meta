library(magrittr)
library(dplyr)
library(tidyr)
library(metap)

#Function that calculates the meta p values using Fisher's Method 
DingMetaAnalysis <- function(full_results, num_p, num_regions) {
  #remove data where there is no valid effectsize value - it will break sumlog calculations
  full_results %<>% filter(!is.na(effectsize))
  
  #p-values are two tailed for negative and positive expression levels 
  #create six one-sided pvalues per gene based on the direction of expression (get a p-value for both directions separately to see which is more significant)
  full_results %<>% mutate(higher_in_MDD_pvalue = two2one(P.Value, invert=1 == sign(-1*effectsize)))
  full_results %<>% mutate(lower_in_MDD_pvalue = two2one(P.Value, invert=1 == sign(1*effectsize)))
  
  full_results %<>% rename(gene_symbol = gene_symbol) ##put in DIngEtAl script
  
  #Find the genes that have p-values in all brain regions
  genes_with_pvalues <- unique(full_results %>% group_by(gene_symbol) %>% summarize(n=n(), unique_regions = length(unique(brain_region))) %>% filter(n==num_p, unique_regions == num_regions) %>% .$gene_symbol)
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
  combined_order = sort(unique(full_results$combined_region_sex))
  combined_order_string = paste(combined_order, collapse="_")
  
  #summarize the direction of expression with '+' and '-' 
  directions <- full_results %>% select(gene_symbol, combined_region_sex, effectsize)
  directions %<>% spread(combined_region_sex, effectsize)
  directions %<>% mutate_at(combined_order, list(~ if_else(. > 0, "+", "-")))
  directions %<>% unite(col = !!combined_order_string, combined_order, sep="", remove = TRUE)
  
  summary_results <- inner_join(directions, summary_results)
  summary_results %<>% arrange(meta_p)
  
  return(summary_results)
}

#This function prepares the data in the correct formatting before performing the meta-analysis
ProcessDingTable<- function(DingTable) {
  Ding_Slim <- DingTable %<>% select(gene_symbol, GeneTitle, contains("MD"))
  Ding_Slim_p <- Ding_Slim %>% select(-contains("_effectsize")) #get rid of the effect size columns to keep only the pvalue columns 

  Ding_Slim_p_female <-Ding_Slim_p %>% select(-contains("_M"))
  #gather the data 
  Ding_Slim_p_female %<>% gather( "Study", "P.Value", 3:ncol(Ding_Slim_p_female))
  Ding_Slim_p_female %<>% mutate(brain_region = ifelse(grepl("ACC", Study), ifelse(grepl("MD3", Study),gsub("MD3_ACC", "ACC.1", Study), gsub("MD2_ACC", "ACC.2", Study)),gsub("^[^_]*_","" ,Study)) ,sex = "female")
  Ding_Slim_p_female %<>% mutate(brain_region = gsub("_F", "", brain_region))
  
  Ding_Slim_p_male <- Ding_Slim_p %>% select(-contains("_F"))
  #gather the data 
  Ding_Slim_p_male %<>% gather( "Study", "P.Value", 3:ncol(Ding_Slim_p_male))
  Ding_Slim_p_male %<>% mutate(brain_region = ifelse(grepl("ACC", Study), ifelse(grepl("MD1", Study),gsub("MD1_ACC", "ACC.1", Study), gsub("MD2_ACC", "ACC.2", Study)),gsub("^[^_]*_","" ,Study)) ,sex = "male")
  Ding_Slim_p_male %<>% mutate(brain_region = gsub("_M", "", brain_region))
  
  #combine male and female p_values together 
  Ding_long <- rbind(Ding_Slim_p_female, Ding_Slim_p_male)
  
  # add in effectsizes 
  Ding_Slim_effectsize <- Ding_Slim %>% select(gene_symbol, GeneTitle, contains("_effectsize"))
  #remove _effectsize from name to match with p.values 
  Ding_Slim_effectsize_female <- Ding_Slim_effectsize %>% select(gene_symbol, GeneTitle,contains("_F_"))
  Ding_Slim_effectsize_female %<>% gather("Study", "effectsize", 3:ncol(Ding_Slim_effectsize_female))
  Ding_Slim_effectsize_female %<>% mutate(Study = gsub("_effectsize", "", Study))
  
  Ding_Slim_effectsize_male <- Ding_Slim_effectsize %>% select(gene_symbol, GeneTitle,contains("_M_"))
  Ding_Slim_effectsize_male %<>% gather("Study", "effectsize", 3:ncol(Ding_Slim_effectsize_male))
  Ding_Slim_effectsize_male %<>% mutate(Study = gsub("_effectsize", "", Study))
  
  #combine male and female effectsize values together  
  Ding_effectsize_long <- rbind(Ding_Slim_effectsize_female, Ding_Slim_effectsize_male)
  
  #combine effectsize with pvalue 
  Ding_long %<>% left_join(Ding_effectsize_long)
  return(Ding_long)
  
}
