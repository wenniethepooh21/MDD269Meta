library(cowplot)
library(ggplot2)
library(googledrive)
library(googlesheets4)
library(here)
library(readr)
library(dplyr)
library(tidyr)
#This script creates a heatmap showing the direction of expression of the top 11 genes in the Ramaker raw results
drawRamaker <- function() {
  top_genes <- c("MANEA","UBE2M","CKB","ITPR3","SPRY2","SAMD5","TMEM106B","ZC3H7B","LST1","ASXL3","HSPA1A") 
  
  #read in Ramaker data
  Ramaker <- read_csv(here("Processed_Data/RamakerEtAl/CombinedCompleteFemaleMaleRamakerTableMagma.csv"))
  Ramaker %<>% filter(gene_symbol %in% top_genes)
  Ramaker %<>% rowwise() %>% mutate(target_region = if_else(target_region == "AnCg", "ACC", target_region),
                                    expression_direction = t*log(P.Value)*-1)
  Ramaker_male <- Ramaker %>% filter(sex == "male")
  Ramaker_female <- Ramaker %>% filter(sex == "female")
  
  #arrange these genes with significant ones together 
  Ramaker_male %<>% mutate(gene_symbol = factor(gene_symbol, levels=rev(c("HSPA1A","ZC3H7B", "SAMD5", "SPRY2", "ITPR3", "MANEA", "UBE2M", "CKB", "TMEM106B","ASXL3", "LST1"))))
  Ramaker_female %<>% mutate(gene_symbol = factor(gene_symbol, levels = rev(c("HSPA1A","ZC3H7B", "SAMD5", "SPRY2", "ITPR3", "MANEA", "UBE2M", "CKB", "TMEM106B","ASXL3", "LST1"))))
  
  ramaker_male_regions <- Ramaker_male$target_region
  ramaker_male_symbol <- Ramaker_male$gene_symbol
  ramaker_male_direction <- Ramaker_male$expression_direction
  
  ramaker_female_regions <- Ramaker_female$target_region
  ramaker_female_symbol <- Ramaker_female$gene_symbol
  ramaker_female_direction <- Ramaker_female$expression_direction
  
  #source file that has the drawExpressionHeat function 
  source(here("R/visualization/heatmaps.R"))
  #function returns the heatmaps of both sexes combined
  Ramaker_plots <- drawExpressionHeat(Ramaker_male, Ramaker_female, ramaker_male_regions, ramaker_male_symbol, ramaker_female_regions, ramaker_female_symbol, ramaker_male_direction, ramaker_female_direction, "ramaker")
  # add title to the combined plots 
  title <- ggdraw() + draw_label("Ramaker Gene Expression", fontface = "bold",size = 15, hjust = 0.35)
  # plot heatmap
  ramaker_heat <- plot_grid(title, Ramaker_plots,ncol = 1, rel_heights = c(0.1, 1))
  ggsave(filename = here('Processed_Data/Meta_Analysis_Results/Heatmaps/top_genes_Ramaker_expression_heatmap.png'), dpi=300, width=12, height=8)
  return(ramaker_heat)
  
}