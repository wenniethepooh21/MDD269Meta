library(cowplot)
library(ggplot2)
library(googledrive)
library(googlesheets4)
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(magrittr)
#This script creates a heatmap showing the direction of expression of the top 11 genes in the Labonte raw results
drawLabonte <- function() {
  top_genes <- c("MANEA","UBE2M","CKB","ITPR3","SPRY2","SAMD5","TMEM106B","ZC3H7B","LST1","ASXL3","HSPA1A") 
  
  #read in Labonte data
  Labonte <- read_csv(here("Processed_Data/LabonteEtAl/CompleteLabonteTableMagma.csv"))
  Labonte %<>% filter(gene_symbol %in% top_genes)
  Labonte %<>% rowwise() %>% mutate(brain_region = if_else(brain_region == "Nac", "nAcc", 
                                                  if_else(brain_region == "Subic", "Sub", 
                                                  if_else(brain_region == "Anterior_Insula", "Ins", brain_region))),
                                    expression_direction = logFC*log(pvalue)*-1)
                      
  Labonte_male <- Labonte %>% filter(sex == "male")
  Labonte_female <- Labonte %>% filter(sex == "female")
  
  
  #arrange these genes with significant ones together 
  Labonte_male %<>% mutate(gene_symbol = factor(gene_symbol, levels=c("HSPA1A","MANEA","TMEM106B","SPRY2","SAMD5","ASXL3","ZC3H7B","UBE2M", "CKB" ,"LST1" ,"ITPR3")))
  Labonte_female %<>% mutate(gene_symbol = factor(gene_symbol, levels=c("HSPA1A","MANEA","TMEM106B","SPRY2","SAMD5","ASXL3","ZC3H7B","UBE2M", "CKB" ,"LST1" ,"ITPR3")))
  
  labonte_male_regions <- Labonte_male$brain_region
  labonte_male_symbol <- Labonte_male$gene_symbol
  labonte_male_direction <- Labonte_male$expression_direction
  
  labonte_female_regions <- Labonte_female$brain_region
  labonte_female_symbol <- Labonte_female$gene_symbol
  labonte_female_direction <- Labonte_female$expression_direction
  
  #source file that has the drawExpressionHeat function 
  source(here("R/visualization/heatmaps.R"))
  #function returns the heatmaps of both sexes combined
  labonte_plots <- drawExpressionHeat(Labonte_male, Labonte_female, labonte_male_regions, labonte_male_symbol, labonte_female_regions, labonte_female_symbol, labonte_male_direction, labonte_female_direction)
  # add title to the combined plots 
  title <- ggdraw() + draw_label("Labonte Gene Expression", fontface = "bold",size = 15)
  # plot heatmap
  labonte_heat <- plot_grid(title, labonte_plots,ncol = 1,rel_heights = c(0.1, 1)) 
  ggsave(filename = here('Processed_Data/Meta_Analysis_Results/Heatmaps/top_genes_Labonte_expression_heatmap.png'), dpi=300, width=12, height=8)
  return(labonte_heat)
  
  
}