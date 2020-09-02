library(here)
library(readr)
library(magrittr)
library(dplyr)
library(tidyr)
library(pheatmap)

howard <- read_csv(here("Processed_Data/ZeiselEtAl/HowardCellTypes_zscore.csv"))
cell_zscore <- read_csv(here("Processed_Data/ZeiselEtAl/full_cell_taxon_zscore.csv"))


howard_gene_symbols <- howard %>% dplyr::select(gene_symbol, mouseGene)
top_genes <- c("MANEA","UBE2M","CKB", "ITPR3","SPRY2","SAMD5","TMEM106B","ZC3H7B","ASXL3","HSPA1A","ZNF184") 
howard_gene_symbols %<>% filter(gene_symbol %in% top_genes) %>% mutate()
top_genes_zscore <- howard_gene_symbols %>% left_join(cell_zscore %>% dplyr::select(mouse_gene, cell_type_taxon, cell_type_taxon_zscore), by = c('mouseGene' = 'mouse_gene'))
top_genes_zscore %<>% dplyr::select(-mouseGene)

long_top_genes_zscore <- spread(top_genes_zscore, gene_symbol, cell_type_taxon_zscore)
long_top_genes_zscore %<>% filter(cell_type_taxon != "Gene detected; No expression measured")

mat <- long_top_genes_zscore %>% dplyr::select(-cell_type_taxon) %>% as.matrix()
rownames(mat) <- long_top_genes_zscore$cell_type_taxon
mat <- t(mat)

#add in LST1 (missing data) as a category 
lst1_df = data.frame("LST1" = c(rep("NA",39)))
rownames(lst1_df) <- long_top_genes_zscore$cell_type_taxon

#change colour to gray
ann_color = list(LST1 = c("NA" = 'gray'))

pheatmap(mat, annotation_col = lst1_df, 
         annotation_colors = ann_color, 
         annotation_legend = FALSE,
         legend_breaks = c(-1, 0,1,2,3,4,5, max(mat)),
         legend_labels = c("-1", "0", "1","2","3","4", "5", "Z-Score\n"),
         angle_col = "45", 
         cellwidth = 9,
         cellheight = 12, 
         fontsize = 4.8,
         filename = here("Results/Heatmaps/cell_type_zscore_pheatmap.png"))

