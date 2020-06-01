library(here)
library(tidyr)
library(edgeR) #loads limma
library(metap)
library(magrittr)
library(dplyr)

#This script contains the meta-analysis functions used to calculate the meta-p and meta directions for the Ramaker dataset

#This function replicates the model used by Ramaker, et al. in their paper: https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-017-0458-5
RamakerDEModel <- function(metadata, read_counts, rawcount_dataframe, regions, full_results) {
	for(target_region in regions) {
	  print(target_region)
	  metadatacp <- metadata

	  #filter for brain region and diagnosis 
	  metadatacp %<>% filter(clinical_diagnosis %in% c("Major_Depression", "Control"), `brain region` == target_region)

	  #filter expression data the same way (the GSM id's that were MDD or Control )
	  rawcount_dataframecp <- rawcount_dataframe[,metadatacp$Sample_geo_accession]
	  
	  # print(paste("Samples remaining: ", nrow(metadatacp)))
	  ## check if rows are in the right order and if there are the same number of people
	  # print(paste("rows line up:", sum(metadatacp$Sample_geo_accession == colnames(rawcount_dataframecp)) == nrow(metadatacp)))
	  
	  #design expression model based on Ramaker et al. of how differential expression calculations
	  #ethnicity + gender  not used in their model
	  design <- model.matrix(~ clinical_diagnosis + `age at death` + `post-mortem interval` + brain_ph + read_counts, metadatacp)
	  
	  # Create DGEList object
	  dge <- DGEList(counts=rawcount_dataframecp)
	  
	  # Calculate normalization factors
	  dge <- calcNormFactors(dge)
	  
	  # Run VOOM
	  v <- voom(dge, plot=FALSE)
	  
	  # Fit linear model
	  fit <- lmFit(v, design)
	  
	  # Run DE
	  fit <- eBayes(fit)
	  
	  # Get results - how well the generated model fits to the data supplied 
	  #Extract a table of the top-ranked genes from a linear model fit
	  #fit is list containing a linear model fit
	  #coef is the column name specifying which coefficient or contrast of the linear model is of interest
	  #fdr method used to adjust the p-values for multiple testing (designed to control the expected proportion of "discoveries" that are false)
	  limma_dataframe <- topTable(fit, coef= "clinical_diagnosisMajor_Depression", adjust="fdr", number=nrow(rawcount_dataframecp))
	  limma_dataframe$gene_symbol <- rownames(limma_dataframe)
	  
	  #t-value is relative to control - postive means increased expression in depression
	  limma_dataframe <- as_tibble(limma_dataframe) %>% dplyr::select(gene_symbol, everything())
	  limma_dataframe %<>% select(gene_symbol, t, P.Value) %>% mutate(target_region = target_region)
	  #add data to full results table
	  full_results %<>% bind_rows(limma_dataframe)  
	}
	
	full_results %<>% arrange(P.Value) #p-value represents the probability of getting the observed expression for this gene within this model (
	return(full_results)
}

#Function that calculates the meta p values using Fisher's method
RamakerMetaAnalysis <- function(full_results, regions){
    #p-values are two tailed for negative and positive expression levels 
  	#create six one-sided pvalues per gene based on the direction of expression (get a p-value for both directions separately to see which is more significant)
  	full_results %<>% mutate(higher_in_MDD_pvalue = two2one(P.Value, invert=1 == sign(-1*t))) #invert p-value if expression is lower (1-p_value) otherwise (pvalue/2)
  	full_results %<>% mutate(lower_in_MDD_pvalue = two2one(P.Value, invert=1 == sign(1*t))) #invert p-value if expression is higher (1-p_value) otherwise (p_value/2)
  	
  	#count the number of occurrence each gene has
  	num_distinct <- full_results %>% select(gene_symbol) %>% group_by(gene_symbol) %>% summarize(total_num = n()) %>% select(total_num) %>% distinct() %>% pull()
  	#perform meta-analysis 
  	#Run Fisher's method grouping genes across brain regions (if there's more than one brain region)
  	  #combining the information in the p-values from different statistical tests to form a single overall test
  	if(num_distinct > 1) {
  		summary_results <- full_results %>% group_by(gene_symbol) %>% summarize(min_p_across_regions = min(P.Value), 
  	                                                      meta_higher_in_MDD_pvalue = sumlog(c(higher_in_MDD_pvalue))$p,
  	                                                      meta_lower_in_MDD_pvalue = sumlog(c(lower_in_MDD_pvalue))$p)
  		} else {
  		  #there's only data for one brain region not need to perform Fisher's method  
  			summary_results <- full_results %>% group_by(gene_symbol) %>% summarize(min_p_across_regions = min(P.Value), 
  	                                                      meta_higher_in_MDD_pvalue = higher_in_MDD_pvalue,
  	                                                      meta_lower_in_MDD_pvalue = lower_in_MDD_pvalue)
  		}
  	#meta p-value is the direction with the more significant p-value, multiply by 2 to change back to 2 sided p-value
  	summary_results %<>% rowwise() %>% mutate(meta_direction = if_else(meta_higher_in_MDD_pvalue < meta_lower_in_MDD_pvalue, 1, -1), meta_p = 2 * min(meta_higher_in_MDD_pvalue, meta_lower_in_MDD_pvalue))
  	
  	#add in individual directions for visualization
  	#handle flipping of male direction for sex-interaction meta-analysis 
  	if ("sex" %in% colnames(full_results)){
    	directions <- full_results %>% select(gene_symbol, target_region, t, sex)
    }else {
      directions <- full_results %>% select(gene_symbol, target_region, t)
    }
  	
  	#summarize the direction of expression with '+' and '-' 
  	directions %<>% spread(target_region, t)
  	directions %<>% mutate_at(regions, list(~ if_else(. > 0, "+", "-")))  
  	colname <- as.character(paste(c(regions, "directions"), collapse = "_"))
  	directions %<>% unite(col = !!colname, regions, sep="", remove = TRUE)
  
  	summary_results <- inner_join(directions, summary_results)
  	summary_results %>% arrange(meta_p)
  
	return(summary_results)
}

