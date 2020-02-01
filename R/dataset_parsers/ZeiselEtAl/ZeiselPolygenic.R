library(homologene)
library(stringr)
library(doMC)
library(readxl)

#Has the linnarssonMatrixMouse and linnarssonMatrixHumanReachable
load('./ProcessedData/ZeiselEtAl/processed_zeisel.Rdata', verbose=TRUE)
descriptions <- read_csv(here::here("ProcessedData", "ZeiselEtAl", "celltype_descriptions.csv"))
cell_types <- read_excel(here::here("data", "ZeiselEtAl", "mmc3.xlsx"), sheet= "7571_0_celltypes_summary_leafor")
cell_types %<>% select(`Cluster name`, Taxonomy_group)

#calculate the area under the roc
auroc_analytic <- function(scores, labels) {
  
  negatives <- which(labels == 0, arr.ind = TRUE)
  scores[negatives] <- 0
  
  p <- sum(scores, na.rm = TRUE)
  nL <- length(labels)
  np <- sum(labels, na.rm = TRUE)
  nn <- nL - np
  
  auroc <- (p/np - (np + 1)/2)/nn
  
  return(auroc)
} 

#function to convert human gene to mouse gene for comparison
convert_genes <- function(input_genes) {
  mouse_genes <- human2mouse(input_genes)
  return(unique(mouse_genes$mouseGene))
}

get_polygenic <- function(cleaned_gene_list) {
    total = length(cleaned_gene_list)
    registerDoMC(cores=2)
    
    #colnames(linnarssonMatrixHumanReachable) these are the cluster id's 
    #linnarssonMatrixHumanReachable is a table that ranks the given gene in each cluster 
    
    #find the cluster id that has the lowest rank for a given gene
    
    #first convert the human gene to mouse gene 
    cleaned_gene_list <- convert_genes(cleaned_gene_list)
    
    #rownames(linnarssonMatrixHumanReachable) are the gene names 
    unique_genes_all <- rownames(linnarssonMatrixMouse)
    
    #get the genes that have a mouse homolog 
    unique_genes_human_reachable <- mouse2human(unique_genes_all)$mouseGene
    unique_genes <- unique_genes_human_reachable
    linnarssonMatrix <- linnarssonMatrixHumanReachable
    
    
    forIndices <- linnarssonMatrix[,1, drop=F]
    forIndices$Gene <- rownames(linnarssonMatrix)
    forIndices %<>% mutate(isTargetGene = Gene %in% cleaned_gene_list)
    targetIndices <- forIndices$isTargetGene
    
    if(length(which(targetIndices == TRUE)) > 0) {
      
      wilcoxTests <- foreach(oneCol=iter(linnarssonMatrix, by='col'), .combine=rbind) %dopar% {
        data.frame(auc = auroc_analytic(oneCol, as.numeric(targetIndices)), 
                   pValue=wilcox.test(oneCol[targetIndices], oneCol[!targetIndices], conf.int = F)$p.value)
      }
      wilcoxTests$cluster_id <- colnames(linnarssonMatrix)
      wilcoxTests <- wilcoxTests[,c(3,1,2)]
      wilcoxTests %<>% mutate(pValue = signif(pValue, digits=3), auc = signif(auc, digits=3), adjusted_P = signif(p.adjust(pValue), digits=3))
      wilcoxTests <- inner_join(descriptions, wilcoxTests)
      wilcoxTests %<>% arrange(-auc)
      wilcoxTests %<>% left_join(cell_types, by = c('cluster_id' = 'Cluster name'))
      wilcoxTests %<>% mutate(Taxonomy_group = coalesce(Taxonomy_group, description))
    } else {
      wilcoxTests <- data.frame(cluster_id=character(),
                                description = character(),
                                auc = double(),
                                pValue = double(),
                                adjusted_P = double(),
                                Taxonomy_group = character(),
                                stringsAsFactors = FALSE)
    }
    
    cat(paste("\nGenes found in data:",sum(cleaned_gene_list %in% unique_genes), "of", total))
    # if(sum(cleaned_gene_list %in% unique_genes) != length(cleaned_gene_list)){
    #   cat(paste("\nGenes not found in data: ", mouse2human(cleaned_gene_list[which(!cleaned_gene_list %in% unique_genes)]))$humanGene)
    # }
    
    return(wilcoxTests)
}






