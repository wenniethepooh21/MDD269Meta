library(metap)
library(readr)
library(magrittr)
library(dplyr)
library(tidyr)
library(here)

#This function identifies the list of p-values for each 269 genes that exist in the given study (used for min_p study analysis)
get_p <- function(howard_df, study_df, study_name){
  study_col <- paste0(study_name,'.Gene')
  output_df <- howard_df %>% dplyr::select(gene_symbol, study_col) %>% 
    left_join(study_df %>% dplyr::select(gene_symbol, min_p_across_regions), by = setNames('gene_symbol', study_col)) %>% 
    mutate(min_p_study = study_name) %>% dplyr::select(-study_col) %>% na.omit() 
  
  return (output_df)
}

###########ToO USE? REFERENCE IN MERGEMETAONSLIMS
##########3
#This function determines if there are any gene:min_p combinations that are the same between studies and updates the min_p_study column if so 
get_min_study <- function(p_table) {
  studies_list <- p_table %>% dplyr::select(min_p_study) %>% distinct()
  num_studies <- studies_list %>% nrow() 
  study <- studies_list %>% pull()
  num_studies_df <- p_table %>% group_by(gene_symbol, min_p_across_regions) %>% summarize(count = n()) %>% arrange(-count)
  unique_studies <- num_studies_df %>% filter(count == 1) %>% left_join(p_table, by = c('gene_symbol' = 'gene_symbol', 'min_p_across_regions' = 'min_p_across_regions'))
  num_duplicates <- nrow(num_studies_df %>% filter(count > 1))
  if(num_duplicates > 0) {
    num_studies_df %<>% left_join(unique_studies)
    num_studies_df %<>% mutate(min_p_study = if_else(count == num_studies, paste(study, collapse = ' '), if_else(count == 1, min_p_study, "multiple studies") ))
  }

  return(num_studies_df)
}


#This function merges the study-specific meta-analysis results into one table filtered for the 269 genes
mergeMetaStudies <- function (H, L, Ldir, D, Ddir, R, Rdir) {
  #All studies were used
  merged_directions <- H %>% left_join(R %>% dplyr::select(gene_symbol, RamakerDir = Rdir), by = c('Ramaker.Gene' = 'gene_symbol'))
  merged_directions %<>% left_join(L %>% dplyr::select(gene_symbol, LabonteDir= Ldir),by = c('Labonte.Gene' = 'gene_symbol'))
  merged_directions %<>% left_join( D%>% dplyr::select(gene_symbol, DingDir= Ddir),by = c('Ding.Gene' = 'gene_symbol'))
  
  #Label which genes Ding filtered out due to low expression and variance
  merged_directions %<>% mutate(DingDir = ifelse(Ding.Gene == "unused", "Filtered_Out", DingDir))
  #Label which genes that we filtered out due to missing p values
  merged_directions %<>% mutate(LabonteDir = ifelse(Labonte.Gene == "UNUSED", "Missing_P_Values", LabonteDir))
  
  #update Labonte data to make gene_symbols match to Howard gene symbols
  index<- which(H$gene_symbol != H$Labonte.Gene & H$Labonte.Gene != "UNUSED" & H$Labonte.Gene != "NA"& H$gene_symbol != "DCDC5")
  oldGeneName <- H$Labonte.Gene[index]
  newGeneName <- H$gene_symbol[index]
  
  for(i in 1:length(index)) {
    L$gene_symbol[which(L$gene_symbol == oldGeneName[i])] <- newGeneName[i]
  }
  
  #update Ding data to make gene_symbols match to Howard gene symbols
  index<- which(H$gene_symbol != H$Ding.Gene & H$Ding.Gene != "UNUSED" & H$Ding.Gene != "NA")
  oldGeneName <- H$Ding.Gene[index]
  newGeneName <- H$gene_symbol[index]
  
  for(i in 1:length(index)) {
    D$gene_symbol[which(D$gene_symbol == oldGeneName[i])] <- newGeneName[i]
  }
  
  #update Ramaker data to make gene_symbols match to Howard gene symbols
  index<- which(H$gene_symbol != H$Ramaker.Gene & H$Ramaker.Gene != "UNUSED" & H$Ramaker.Gene != "NA"& H$gene_symbol != "DCDC5")
  oldGeneName <- H$Ramaker.Gene[index]
  newGeneName <- H$gene_symbol[index]
  
  for(i in 1:length(index)) {
    R$gene_symbol[which(R$gene_symbol == oldGeneName[i])] <- newGeneName[i]
  }
  merged_p <- bind_rows(L, D, R)

  #merge two datasets together
  merged_directions %<>% dplyr::select(-`Match Across Studies?`,-Ramaker.Gene, -Labonte.Gene, - Ding.Gene)
  
  merged_p %<>% mutate(gene_symbol = gsub("C([X0-9]+)ORF([0-9]+)", "C\\1orf\\2", gene_symbol))
  
  merged_table <- left_join(merged_directions, merged_p)
  return(merged_table)
}

#This function takes the merged study-specific meta-analysis results and performes our meta-analysis combined using Fisher's method 
MetaAnalysis <- function(merged_table,minp_table=FALSE){
  
  #get just minp and the two one-sided meta p-values for each row
  merged_p <- merged_table %>% dplyr::select(gene_symbol, min_p_across_regions, meta_lower_in_MDD_pvalue, meta_higher_in_MDD_pvalue)
  merged_directions <- merged_table %>% dplyr:: select(gene_symbol, Howard_pvalue, gene_name, RamakerDir, LabonteDir, DingDir)
  
  meta_merged_p <- merged_p %>% group_by(gene_symbol) %>% summarize(list_of_meta_higher_in_MDD_pvalue = list(meta_higher_in_MDD_pvalue), 
                                                                    list_of_meta_lower_in_MDD_pvalue = list(meta_lower_in_MDD_pvalue), count_of_pvalues=n())
  min_p <- merged_p %>% group_by(gene_symbol) %>% summarize(minp = min(min_p_across_regions))
  
  if(minp_table!=FALSE) {
    min_p %<>% left_join(minp_table, by = c('gene_symbol' = 'gene_symbol', 'minp' = 'min_p_across_regions'))
  }
  #had trouble doing this with if_else
  twoOrMore <- meta_merged_p %>% filter(count_of_pvalues > 1) %>% rowwise() %>% 
    mutate(meta_meta_higher_in_MDD_pvalue = sumlog(c(list_of_meta_higher_in_MDD_pvalue))$p,
           meta_meta_lower_in_MDD_pvalue=  sumlog(c(list_of_meta_lower_in_MDD_pvalue))$p
    ) %>% dplyr::select(gene_symbol, meta_meta_higher_in_MDD_pvalue, meta_meta_lower_in_MDD_pvalue)
  
  singles <- meta_merged_p %>% filter(count_of_pvalues == 1) %>% 
    mutate(meta_meta_higher_in_MDD_pvalue = unlist(list_of_meta_higher_in_MDD_pvalue),
           meta_meta_lower_in_MDD_pvalue =  unlist(list_of_meta_lower_in_MDD_pvalue)
    ) %>% dplyr::select(gene_symbol, meta_meta_higher_in_MDD_pvalue, meta_meta_lower_in_MDD_pvalue)
  meta_merged_p <- bind_rows(twoOrMore, singles)
  
  #the doubling causes some meta_p values be above 1, check
  meta_merged_p %<>% rowwise() %>% mutate(meta_direction = if_else(meta_meta_higher_in_MDD_pvalue < meta_meta_lower_in_MDD_pvalue, '+', '-'), meta_p = if_else(2 * min(meta_meta_higher_in_MDD_pvalue, meta_meta_lower_in_MDD_pvalue) > 1, 1, 2 * min(meta_meta_higher_in_MDD_pvalue, meta_meta_lower_in_MDD_pvalue)))
  
  meta_merged_p %<>% mutate(meta_direction_val= if_else(meta_meta_higher_in_MDD_pvalue < meta_meta_lower_in_MDD_pvalue, meta_meta_higher_in_MDD_pvalue, meta_meta_lower_in_MDD_pvalue))
  meta_merged_p %<>% rowwise() %>% mutate(signed_log_meta = sign(meta_direction_val)*log10(as.numeric(meta_p))*-1)
  meta_merged_p %<>% dplyr::select(-meta_direction_val)
  
  full_table <- inner_join(meta_merged_p,min_p)
  full_table <- left_join(merged_directions,full_table) %>% distinct()
  
  #copy the values from DCDC1 to DCDC5 (same gene.. but have different p values for Howard)
  keep <- full_table %>% filter(gene_symbol == "DCDC1")
  keep_index <- keep[-c(1:6)]
  dcdc5 <- full_table %>% filter(gene_symbol == "DCDC5") %>% dplyr::select(1:6) %>% cbind(keep_index)
  full_table %<>% filter(gene_symbol != "DCDC5")
  full_table %<>% rbind(dcdc5) 
  
  full_table %<>% mutate(Corrected_p = meta_p*269)
  
  full_table %<>% mutate(Bonferroni_meta_p = p.adjust(meta_p, method = "bonferroni") ) 
  full_table %<>% arrange(Corrected_p)
  # three significant digits 
  full_table$Corrected_p <- signif(as.numeric(full_table$Corrected_p),digits=3)
  full_table$Bonferroni_meta_p <- signif(as.numeric(full_table$Bonferroni_meta_p),digits=3)
  
  if(length(is.na(full_table$Howard_pvalue)) == 269){
    full_table %<>% mutate(Spearman_corr = cor.test(full_table$Howard_pvalue, full_table$meta_p, use = 'pairwise.complete.obs', method = "spearman")$estimate)
    full_table %<>% mutate(Spearman_p = cor.test(full_table$Howard_pvalue, full_table$meta_p, use = 'pairwise.complete.obs', method = "spearman")$p.value)
    full_table$Spearman_corr <- signif(as.numeric(full_table$Spearman_corr),digits=3)
    full_table$Spearman_p <- signif(as.numeric(full_table$Spearman_p),digits=3)
    
  }
  full_table$meta_p <- signif(as.numeric(full_table$meta_p),digits=3)
  return(full_table)
}

#This function calculates combines the calculated empirical p-values of each gene across the studies to determine the meta-empirical p-value or the meta percentile rank
GenomeRank <- function(merged_table) {
  
  merged_p <- merged_table %>% dplyr::select(gene_symbol, min_p_across_regions, meta_Down, meta_Up)
  merged_directions <- merged_table %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name, RamakerDir, LabonteDir, DingDir)
  
  #get just minp and the two one-sided meta p-values for each row
  meta_merged_p <- merged_p %>% group_by(gene_symbol) %>% summarize(list_of_meta_higher_in_MDD_genome_percentile_rank = list(meta_Up), 
                                                                    list_of_meta_lower_in_MDD_genome_percentile_rank = list(meta_Down), count_of_pvalues=n())
  min_p <- merged_p %>% group_by(gene_symbol) %>% summarize(minp = min(min_p_across_regions))
  
  #had trouble doing this with if_else
  twoOrMore <- meta_merged_p %>% filter(count_of_pvalues > 1) %>% rowwise() %>% 
    mutate(meta_meta_higher_in_MDD_genome_percentile_rank = sumlog(c(list_of_meta_higher_in_MDD_genome_percentile_rank))$p,
           meta_meta_lower_in_MDD_genome_percentile_rank=  sumlog(c(list_of_meta_lower_in_MDD_genome_percentile_rank))$p
    ) %>% dplyr::select(gene_symbol, meta_meta_higher_in_MDD_genome_percentile_rank, meta_meta_lower_in_MDD_genome_percentile_rank)
  
  singles <- meta_merged_p %>% filter(count_of_pvalues == 1) %>% 
    mutate(meta_meta_higher_in_MDD_genome_percentile_rank = unlist(list_of_meta_higher_in_MDD_genome_percentile_rank),
           meta_meta_lower_in_MDD_genome_percentile_rank =  unlist(list_of_meta_lower_in_MDD_genome_percentile_rank)
    ) %>% dplyr::select(gene_symbol, meta_meta_higher_in_MDD_genome_percentile_rank, meta_meta_lower_in_MDD_genome_percentile_rank)
  meta_merged_p <- bind_rows(twoOrMore, singles)
  
  #the doubling causes some meta_p values be above 1, check
  meta_merged_p %<>% rowwise() %>% mutate(meta_direction = if_else(meta_meta_higher_in_MDD_genome_percentile_rank < meta_meta_lower_in_MDD_genome_percentile_rank, '+', '-'), meta_empirical_p = if_else(2 * min(meta_meta_higher_in_MDD_genome_percentile_rank, meta_meta_lower_in_MDD_genome_percentile_rank) > 1, 1, 2 * min(meta_meta_higher_in_MDD_genome_percentile_rank, meta_meta_lower_in_MDD_genome_percentile_rank)))
  
  
  full_table <- inner_join(meta_merged_p,min_p)
  full_table <- left_join(merged_directions,full_table) %>% distinct()
  
  #copy the values from DCDC1 to DCDC5 (same gene.. but have different p values for Howard)
  keep <- full_table %>% filter(gene_symbol == "DCDC1")
  keep_index <- keep[-c(1:6)]
  full_table[which(full_table$gene_symbol == "DCDC5"),c(7:ncol(full_table))]<- keep_index
  
  if(length(is.na(full_table$Howard_pvalue)) == 269){
    full_table %<>% mutate(Spearman_corr = cor.test(full_table$Howard_pvalue, full_table$meta_empirical_p, use = 'pairwise.complete.obs', method = "spearman")$estimate)
    full_table %<>% mutate(Spearman_p = cor.test(full_table$Howard_pvalue, full_table$meta_empirical_p, use = 'pairwise.complete.obs', method = "spearman")$p.value)
    full_table$Spearman_corr <- signif(as.numeric(full_table$Spearman_corr),digits=3)
    full_table$Spearman_p <- signif(as.numeric(full_table$Spearman_p),digits=3)
  }  
  full_table %<>% mutate(Corrected_meta_empirical_p = meta_empirical_p*269) #Bonferroni correction
  
  full_table %<>% mutate(Bonferroni_meta_empirical_p = p.adjust(meta_empirical_p, method = "bonferroni"))
  
  full_table %<>% arrange(Corrected_meta_empirical_p)
  full_table$Corrected_meta_empirical_p <- signif(as.numeric(full_table$Corrected_meta_empirical_p),digits=3)
  full_table$Bonferroni_meta_empirical_p <- signif(as.numeric(full_table$Bonferroni_meta_empirical_p),digits=3)
  full_table$meta_empirical_p <- signif(as.numeric(full_table$meta_empirical_p),digits=3)
  return(full_table)
}

#performs the meta-analysis or genome ranking on all magma genes (not filtered for the 269)
mergeMagmaMetaRank<- function(m, L, Ldir, D, Ddir, R, Rdir, analysis){
  
  #All studies were used
  merged_directions <- left_join(m, R %>% dplyr::select(gene_symbol, RamakerDir = Rdir), by = c('Ramaker_genes' = 'gene_symbol'))
  merged_directions %<>% left_join(L %>% dplyr::select(gene_symbol, LabonteDir= Ldir),by = c('Labonte_genes' = 'gene_symbol'))
  merged_directions %<>% left_join(D %>% dplyr::select(gene_symbol, DingDir= Ddir),by = c('Ding_genes' = 'gene_symbol')) %>% dplyr::select(-Previous_Symbol, - Updated_Symbol, - Labonte_genes, - Ding_genes, -Ramaker_genes) %>% distinct()
  
  Labonte <- left_join(m %>% dplyr::select(gene_symbol,Labonte_genes) %>% na.omit() , L, by = c('Labonte_genes' = 'gene_symbol')) %>% dplyr::select(-Labonte_genes) %>% distinct() %>% filter(!is.na(meta_p))
  Ding <- left_join(m %>% dplyr::select(gene_symbol,Ding_genes), D, by = c('Ding_genes' = 'gene_symbol')) %>% dplyr::select(-Ding_genes) %>% distinct() %>% filter(!is.na(meta_p))
  Ramaker <- left_join(m %>% dplyr::select(gene_symbol,Ramaker_genes), R, by = c('Ramaker_genes' = 'gene_symbol')) %>% dplyr::select(-Ramaker_genes) %>% distinct() %>% filter(!is.na(meta_p))
  
  merged_p <- bind_rows(Labonte, Ding, Ramaker)
  
  merged_table <- left_join(merged_directions, merged_p) %>% distinct()
  merged_table %<>% mutate(Howard_pvalue = NA, gene_name = NA)
  
  if(analysis == "meta") {
    merged_table %<>% MetaAnalysis() %>% dplyr::select(-Howard_pvalue, -gene_name)
  } else {
    merged_table %<>% GenomeRank() %>% dplyr::select(-Howard_pvalue, -gene_name)
  }
  return(merged_table)
  
}


























