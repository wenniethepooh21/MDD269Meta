library(here)
library(tidyr)
library(metap)
library(GEOquery)
library(magrittr)
library(dplyr)
library(readr)
library(limma)
library(edgeR)
library(stringr)

#based on signature.R by the biojupies team
#https://github.com/MaayanLab/biojupies-plugins/blob/1024a6ed702ad8b0958d4ccdd2afe89cbe493a51/library/core_scripts/signature/signature.R
#from https://amp.pharm.mssm.edu/biojupies/notebook/tlNmgbMxY (just load biojupies with GSE80655)

#Method for getting sequence counts
#"Raw RNA-seq data for GEO dataset GSE80655 was downloaded from the SRA database (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80655) and quantified to gene-level counts using the ARCHS4 pipeline (Lachmann et al., 2017). Gene counts were downloaded from the ARCHS4 gene expression matrix v6. For more information about ARCHS4, as well as free access to the quantified gene expression matrix, visit the project home page at the following URL: http://amp.pharm.mssm.edu/archs4/download.html."
#(Biojupies)


RamakerMeta <- function(metadata, read_counts, rawcount_dataframe, regions, full_results) {
	for(target_region in regions) {
	  print(target_region)

	  metadatacp <- metadata

	  #filter for brain region and diagnosis
	  metadatacp %<>% filter(clinical_diagnosis %in% c("Major_Depression", "Control"), `brain region` == target_region)

	  #filter expression data the same way
	  rawcount_dataframecp <- rawcount_dataframe[,metadatacp$Sample_geo_accession]
	  
	  print(paste("Samples remaining: ", nrow(metadatacp)))
	  #check if rows are in the right order
	  print(paste("rows line up:", sum(metadatacp$Sample_geo_accession == colnames(rawcount_dataframecp)) == nrow(metadatacp)))
	  
	  
	  #design matrix
	  #model based on Ramaker et al. 
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
	  
	  # Get results
	  limma_dataframe <- topTable(fit, coef= "clinical_diagnosisMajor_Depression", adjust="fdr", number=nrow(rawcount_dataframecp))
	  limma_dataframe$gene_symbol <- rownames(limma_dataframe)
	  
	  #t-value is relative to control - postive means increased expression in depression
	  limma_dataframe <- as_tibble(limma_dataframe) %>% dplyr::select(gene_symbol, everything())
	  limma_dataframe %<>% dplyr::select(gene_symbol, t, P.Value) %>% mutate(target_region = target_region)
	  full_results <- bind_rows(full_results, limma_dataframe)  
	}
	
	full_results %<>% arrange(P.Value)
	return(full_results)
}

RamakerAnalysis <- function(full_results, regions){
  	#create six one-sided pvalues per gene
  	full_results %<>% mutate(higher_in_MDD_pvalue = two2one(P.Value, invert=1 == sign(-1*t)))
  	full_results %<>% mutate(lower_in_MDD_pvalue = two2one(P.Value, invert=1 == sign(1*t)))
  	
  	if(length(regions) > 1) {
  		summary_results <- full_results %>% group_by(gene_symbol) %>% summarize(min_p_across_regions = min(P.Value), 
  	                                                      meta_higher_in_MDD_pvalue = sumlog(c(higher_in_MDD_pvalue))$p,
  	                                                      meta_lower_in_MDD_pvalue = sumlog(c(lower_in_MDD_pvalue))$p)
  		} else {
  			summary_results <- full_results %>% group_by(gene_symbol) %>% summarize(min_p_across_regions = min(P.Value), 
  	                                                      meta_higher_in_MDD_pvalue = higher_in_MDD_pvalue,
  	                                                      meta_lower_in_MDD_pvalue = lower_in_MDD_pvalue)
  		}
  	#convert from two one-sided meta pvalues to one
  	summary_results %<>% rowwise() %>% mutate(meta_direction = if_else(meta_higher_in_MDD_pvalue < meta_lower_in_MDD_pvalue, 1, -1), meta_p = 2 * min(meta_higher_in_MDD_pvalue, meta_lower_in_MDD_pvalue))
  	
  
  	#add in individual directions for visualization
  	#handle flipping of male direction analysis
  	if ("sex" %in% colnames(full_results)){
    	directions <- full_results %>% dplyr::select(gene_symbol, target_region, t, sex)
    }else {
      directions <- full_results %>% dplyr::select(gene_symbol, target_region, t)
    }
  	
  	directions %<>% spread(target_region, t)
  	directions %<>% mutate_at(regions, list(~ if_else(. > 0, "+", "-")))  
  	colname = as.character(paste(c(regions, "directions"), collapse = "_"))
  	directions %<>% unite(col = !!colname, regions, sep="", remove = TRUE)
  
  	summary_results <- inner_join(directions, summary_results)
  	summary_results %>% arrange(meta_p)
  
  	#code to check direction FAM101B is higher in depression
  	#mean(unlist(rawcount_dataframe["FAM101B", metadata %>% filter(clinical_diagnosis != "Major_Depression") %>% .$Sample_geo_accession]))

	return(summary_results)
}

flipDirections <- function(dataset, gender) {
  flipped_dataset <- dataset %>% filter(sex == gender)
  unflipped_dataset <- dataset %>% filter(sex != gender)
  
  get_col <- flipped_dataset %>% select(contains("directions"))
  flipped_dataset[[names(get_col)]]<- gsub("+", "_",flipped_dataset[[names(get_col)]], fixed = TRUE)
  flipped_dataset[[names(get_col)]]<- gsub("-", "+",flipped_dataset[[names(get_col)]], fixed = TRUE)
  flipped_dataset[[names(get_col)]]<- gsub("_", "-",flipped_dataset[[names(get_col)]], fixed = TRUE)
  
  new_dataset <- rbind(flipped_dataset, unflipped_dataset)
  return(new_dataset)
}
getRank <- function(dataset, num_genes){
  dataset %<>% mutate(meta_Up = NA, meta_Down = NA)
  Ramaker_Up <- dataset %>% filter(meta_direction == +1)
  Ramaker_Down <- dataset %>% filter(meta_direction == -1) 
  for (i in 1:length(dataset$gene_symbol)) {
    gene_meta <- dataset$meta_p[i]
    direction_val <- dataset$meta_direction[i]
    if(direction_val == +1) {
      dataset$meta_Up[i] <- Ramaker_Up %>% filter(meta_p <= gene_meta) %>% dplyr::select(gene_symbol) %>% nrow() / num_genes
      dataset$meta_Down[i] <- (Ramaker_Down %>% filter(meta_p <= gene_meta) %>% dplyr::select(gene_symbol) %>% nrow() + 1) / (num_genes + 1)
      
      #this current gene will be in the down regulation list, offset it
    }else if (direction_val == -1) {
      dataset$meta_Up[i] <- (Ramaker_Up %>% filter(meta_p <= gene_meta) %>% dplyr::select(gene_symbol) %>% nrow() + 1) / (num_genes + 1)
      dataset$meta_Down[i] <- Ramaker_Down %>% filter(meta_p <= gene_meta) %>% dplyr::select(gene_symbol) %>% nrow() / num_genes
    }
  }
  dataset %<>% rowwise() %>% mutate(genome_percentile_rank = if_else(2 * min(meta_Up, meta_Down) > 1, 1, 2 * min(meta_Up, meta_Down)))
  return(dataset)
}