library(reshape2)
library(magrittr)
library(dplyr)
library(tidyverse)
library(devtools)
library(metap)

getRankFishers <- function(dataset, num_genes){
  dataset %<>% mutate(meta_Up = NA, meta_Down = NA)
  Ding_Up <- dataset %>% filter(meta_direction == +1)
  Ding_Down <- dataset %>% filter(meta_direction == -1)
  
  for (i in 1:length(dataset$gene_symbol)) {
    gene_meta <- dataset$meta_p[i]
    direction_val <- dataset$meta_direction[i]
    
    #this current gene will be in the up regulated list, offset it
    if(direction_val == +1) {
      dataset$meta_Up[i] <- Ding_Up %>% filter(meta_p <= gene_meta) %>% dplyr::select(gene_symbol) %>% nrow() / num_genes
      dataset$meta_Down[i] <- (Ding_Down %>% filter(meta_p <= gene_meta) %>% dplyr::select(gene_symbol) %>% nrow() + 1) / (num_genes + 1)
      
      #this current gene will be in the down regulation list, offset it
    }else if (direction_val == -1) {
      dataset$meta_Up[i] <- (Ding_Up %>% filter(meta_p <= gene_meta) %>% dplyr::select(gene_symbol) %>% nrow() + 1) / (num_genes + 1)
      dataset$meta_Down[i] <- Ding_Down %>% filter(meta_p <= gene_meta) %>% dplyr::select(gene_symbol) %>% nrow() / num_genes
    }
    
  }
  dataset %<>% rowwise() %>% mutate(genome_percentile_rank = if_else(2 * min(meta_Up, meta_Down) > 1, 1, 2 * min(meta_Up, meta_Down)))
  return(dataset)
}

#Used in Fisher Analysis
newDingMeta <- function(full_results, numP, numRegions) {
  #remove data where there is no valid effectsize value - it will break sumlog calculations
  full_results %<>% filter(!is.na(effectsize))
  
  #create six one-sided pvalues per gene
  full_results %<>% mutate(higher_in_MDD_pvalue = two2one(P.Value, invert=1 == sign(-1*effectsize)))
  full_results %<>% mutate(lower_in_MDD_pvalue = two2one(P.Value, invert=1 == sign(1*effectsize)))
  
  full_results %<>% dplyr::rename(gene_symbol = SYMBOL)
  
  #maybe deal with these later
  genes_with_pvalues <- unique(full_results %>% group_by(gene_symbol) %>% summarize(n=n(), unique_regions = length(unique(brain_region))) %>% filter(n==numP, unique_regions == numRegions) %>% .$gene_symbol)
  
  full_results %<>% filter(gene_symbol %in% genes_with_pvalues)
  
  summary_results <- full_results %>% group_by(gene_symbol) %>% summarize(min_p_across_regions = min(P.Value), 
                                                                            meta_higher_in_MDD_pvalue = sumlog(c(higher_in_MDD_pvalue))$p,
                                                                            meta_lower_in_MDD_pvalue = sumlog(c(lower_in_MDD_pvalue))$p)

  #convert from two one-sided meta pvalues to one
  summary_results %<>% rowwise() %>% mutate(meta_direction = if_else(meta_higher_in_MDD_pvalue < meta_lower_in_MDD_pvalue, 1, -1), meta_p = 2 * min(meta_higher_in_MDD_pvalue, meta_lower_in_MDD_pvalue))
  
  full_results %<>% mutate(combined_region_sex = paste(sex, brain_region, sep="_"))
  
  combined_order = sort(unique(full_results$combined_region_sex))
  combined_order_string = paste(combined_order, collapse="_")
  
  #add in individual directions for visualization
  directions <- full_results %>% dplyr::select(gene_symbol, combined_region_sex, effectsize)
  #spread seems to do alphabetical sorting here too
  directions %<>% spread(combined_region_sex, effectsize)
  directions %<>% mutate_at(combined_order, list(~ if_else(. > 0, "+", "-")))
  directions %<>% unite(col = !!combined_order_string, combined_order, sep="", remove = TRUE)
  
  summary_results <- inner_join(directions, summary_results)
  summary_results %<>% arrange(meta_p)
  
  
  #code to check direction FAM101B is higher in depression
  #mean(unlist(rawcount_dataframe["FAM101B", metadata %>% filter(clinical_diagnosis != "Major_Depression") %>% .$Sample_geo_accession]))
  
  return(summary_results)
}
processDingTable<- function(DingTable) {
  Ding_Slim <- DingTable %<>% select(SYMBOL, GeneTitle, contains("MD"))
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
  Ding_Slim_effectsize <- Ding_Slim %>% select(SYMBOL, GeneTitle, contains("_effectsize"))
  #remove _effectsize from name to match with p.values 
  Ding_Slim_effectsize_female <- Ding_Slim_effectsize %>% select(SYMBOL, GeneTitle,contains("_F_"))
  Ding_Slim_effectsize_female %<>% gather("Study", "effectsize", 3:ncol(Ding_Slim_effectsize_female))
  Ding_Slim_effectsize_female %<>% mutate(Study = gsub("_effectsize", "", Study))
  
  Ding_Slim_effectsize_male <- Ding_Slim_effectsize %>% select(SYMBOL, GeneTitle,contains("_M_"))
  Ding_Slim_effectsize_male %<>% gather("Study", "effectsize", 3:ncol(Ding_Slim_effectsize_male))
  Ding_Slim_effectsize_male %<>% mutate(Study = gsub("_effectsize", "", Study))
  
  #combine male and female effectsize values together  
  Ding_effectsize_long <- rbind(Ding_Slim_effectsize_female, Ding_Slim_effectsize_male)
  
  #combine effectsize with pvalue 
  Ding_long %<>% left_join(Ding_effectsize_long)
  # Ding_long %>% left_join(Ding_effectsize_long,by = c('SYMBOL' = 'SYMBOL', 'GeneTitle' = 'GeneTitle', 'Study' = 'Study'))
  
  return(Ding_long)
  
}

DingMeta <- function(tableName) {
	for_visualization <- tableName %>% dplyr::select(SYMBOL, contains("_effectsize"))
	for_visualization %<>% mutate_at(vars(contains("effectsize")), list(~ if_else(. > 0, "+", "-")))
	#switch to base R
	for_visualization <- as.data.frame(for_visualization)
	colnames(for_visualization) <- gsub("^[^_]*_", "", colnames(for_visualization))
	f_indx <- for_visualization[which(grepl('_F_', colnames(for_visualization)))]
	f_indx <- f_indx[,sort(colnames(f_indx))]
	m_indx <- for_visualization[which(grepl('_M_', colnames(for_visualization)))]
	m_indx <- m_indx[,sort(colnames(m_indx))]

	sym_indx <- for_visualization[which(grepl('SYMBOL', colnames(for_visualization)))]
	for_visualization <- cbind(sym_indx,f_indx,m_indx)	

	study_string <- paste(colnames(for_visualization), collapse="_")
	study_string <- gsub("SYMBOL_", "", study_string)
	study_string <- gsub("_effectsize", "", study_string)
	for_visualization <- as_tibble(for_visualization)
	# for_visualization %<>% select(gene_symbol = SYMBOL, everything())

	for_visualization %<>% unite(col = !!study_string, contains("effectsize"), sep="", remove = TRUE)
	for_visualization %<>% dplyr::rename(gene_symbol = SYMBOL)

	#get minimum p-value
	min_pvalue <- tableName %>% dplyr::select(gene_symbol = SYMBOL, matches("_[MF]$"))
	min_pvalue %<>% as.data.frame() %>% melt() %>% as_tibble()
	min_pvalue %<>% group_by(gene_symbol) %>% summarise(min_p_across_regions = min(value))

	#get direction via effect size to avoid ties
	general_direction <- tableName %>% dplyr::select(gene_symbol = SYMBOL, contains("effectsize"))
	general_direction %<>% as.data.frame() %>% melt() %>% as_tibble()
	general_direction %<>% group_by(gene_symbol) %>% summarise(direction_from_effectsizes = sign(sum(value,na.rm=T)))
		

	### writing of slimmed/updated version
	Ding_summary_results <- tableName %>% dplyr::select(gene_symbol = SYMBOL, meta_p)
	Ding_summary_results %<>% filter(!is.na(meta_p))
	#direction via effect size sum to deal with ties
	Ding_summary_results <- inner_join(Ding_summary_results, general_direction)

	Ding_summary_results %<>% mutate(meta_higher_in_MDD_pvalue = two2one(meta_p, invert = -1 == direction_from_effectsizes))
	Ding_summary_results %<>% mutate(meta_lower_in_MDD_pvalue = two2one(meta_p, invert = 1 == direction_from_effectsizes))

	#add in minp and visualization
	Ding_summary_results <- inner_join(min_pvalue, Ding_summary_results)
	Ding_summary_results <- inner_join(for_visualization, Ding_summary_results)

	Ding_summary_results %>% arrange(meta_p)

	Ding_summary_results %<>% mutate(!!study_string := paste0("\'", !!rlang::sym(study_string)))
	

	return(Ding_summary_results)
}





getRank <- function(dataset, num_genes){
  dataset %<>% mutate(meta_Up = NA, meta_Down = NA)
  Ding_Up <- dataset %>% filter(direction_from_effectsizes == +1)
  Ding_Down <- dataset %>% filter(direction_from_effectsizes == -1)
  
  for (i in 1:length(dataset$gene_symbol)) {
    gene_meta <- dataset$meta_p[i]
    direction_val <- dataset$direction_from_effectsizes[i]
    
    #this current gene will be in the up regulated list, offset it
    if(direction_val == +1) {
      dataset$meta_Up[i] <- Ding_Up %>% filter(meta_p <= gene_meta) %>% dplyr::select(gene_symbol) %>% nrow() / num_genes
      dataset$meta_Down[i] <- (Ding_Down %>% filter(meta_p <= gene_meta) %>% dplyr::select(gene_symbol) %>% nrow() + 1) / (num_genes + 1)
      
      #this current gene will be in the down regulation list, offset it
    }else if (direction_val == -1) {
      dataset$meta_Up[i] <- (Ding_Up %>% filter(meta_p <= gene_meta) %>% dplyr::select(gene_symbol) %>% nrow() + 1) / (num_genes + 1)
      dataset$meta_Down[i] <- Ding_Down %>% filter(meta_p <= gene_meta) %>% dplyr::select(gene_symbol) %>% nrow() / num_genes
    }
    
  }
  dataset %<>% rowwise() %>% mutate(genome_percentile_rank = if_else(2 * min(meta_Up, meta_Down) > 1, 1, 2 * min(meta_Up, meta_Down)))
  return(dataset)
}
