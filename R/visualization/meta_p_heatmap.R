library(ggplot2)
library(googledrive)
library(googlesheets4)
library(dplyr)
library(tidyr)
library(viridis)
#This script creates a heatmap to display the corrected meta p-values of our top 11 genes and how they ranked in each meta-analysis

drawMetaHeat <- function() {
  top_genes <- c("MANEA","UBE2M","CKB","ITPR3","SPRY2","SAMD5","TMEM106B","ZC3H7B","LST1","ASXL3","HSPA1A") 
  
  # drive_auth() #authenticate gmail
  sheets_auth(token = drive_token())
  
  full <- read_sheet(drive_get("Meta_Analysis"), sheet = 'Full_Meta_Analysis') 
  female<-read_sheet(drive_get("Meta_Analysis"), sheet = 'Female_Meta_Analysis')
  male <-read_sheet(drive_get("Meta_Analysis"), sheet = 'Male_Meta_Analysis')
  cortical <- read_sheet(drive_get("Meta_Analysis"), sheet = 'Cortical_Meta_Analysis')
  SIfull <- read_sheet(drive_get("Sex_Interaction_Meta_Analysis"), sheet = 'Full_Meta_Analysis')
  SIcortical <- read_sheet(drive_get("Sex-Interaction-Meta-Analysis"), sheet = 'Cortical_Meta_Analysis')
  
  fulltable <- full %>% select(gene_symbol, Bonferroni_Correction) %>% rename(Full = Bonferroni_Correction)
  femaletable <-female %>% select(gene_symbol, Bonferroni_Correction)  %>% rename(Female = Bonferroni_Correction)
  maletable <- male %>% select(gene_symbol, Bonferroni_Correction) %>% rename(Male = Bonferroni_Correction)
  corticaltable <- cortical %>% select(gene_symbol, Bonferroni_Correction) %>% rename(Cortical = Bonferroni_Correction)
  SIfulltable <- SIfull%>% select(gene_symbol, Bonferroni_Correction) %>% rename(SI_Full = Bonferroni_Correction)
  SIcorticaltable <- SIcortical%>% select(gene_symbol, Bonferroni_Correction) %>% rename(SI_Cortical = Bonferroni_Correction)
  
  #join the results into one table
  full_meta_analysis_results <- left_join(fulltable, femaletable) %>% left_join(corticaltable) %>% left_join(maletable) %>% left_join(SIfulltable) %>% left_join(SIcorticaltable)
  full_meta_analysis_results %<>% gather("Analysis", "Bonferroni_Correction", 2:ncol(full_meta_analysis_results))
  full_meta_analysis_results %<>% filter(gene_symbol %in% top_genes)
  full_meta_analysis_results %<>% mutate(Bonferroni_Correction = as.numeric(Bonferroni_Correction))
  
  # test <- full_meta_analysis_results %>% mutate(Bonferroni_Correction = if_else(Bonferroni_Correction > 0.05, 0.05, Bonferroni_Correction))
  
  full_meta_analysis_results %<>% mutate(gene_symbol = factor(gene_symbol, levels=rev(c("HSPA1A","ZC3H7B", "SAMD5", "SPRY2", "ITPR3", "MANEA", "UBE2M", "CKB", "TMEM106B","ASXL3", "LST1"))),
                    Analysis = factor(Analysis, levels = c("Full", "Cortical", "Male", "Female", "SI_Cortical", "SI_Full")))
  
  #set the scale to range from 0 to 0.1
  full_meta_analysis_results %<>% mutate(Bonferroni_Correction = Bonferroni_Correction/10)
  ## used for breaking the corrected p-values into ranges
  # full_meta_analysis_results$bin <- cut(full_meta_analysis_results$Bonferroni_Correction, breaks = c(0, 0.0005, 0.001,0.005,0.01, 0.05, 0.1, 0.5, 5, 20,220),
  #                labels = c("0.-0.0005", "0.0005-0.001","0.001-0.005","0.005-0.01", "0.01-0.05", "0.05 - 0.1", "0.1-0.5", "0.5-5", "5-20","20-220"),
  #                include.lowest = TRUE)
  
  #arrange these genes with significant ones together 
  analysis_type <- full_meta_analysis_results$Analysis
  symbol <- full_meta_analysis_results$gene_symbol
  p_val <- full_meta_analysis_results$Bonferroni_Correction
  
  #heatmap
  meta_plot <- ggplot(full_meta_analysis_results, aes(x=analysis_type,y=symbol)) + 
    geom_tile(aes(fill = p_val), colour='black') +
    labs(x="Meta-Analysis", fill = "Corrected\np-value") +
    ylab('') +
    scale_fill_viridis(option = "D", direction = -1) + 
    theme(axis.title.x = element_text(size = 11),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 10,angle = 20, hjust=0.5,vjust=1), 
          axis.text.y = element_text(face = "italic",size = 10), 
          plot.title = element_text(size = 20,face = "bold",hjust = 0.5),
          panel.background = element_blank())  # center the title 

  print(meta_plot)
  
  
  
  
  title <- ggdraw() + draw_label("Meta p-values of Top 11 Genes", fontface = "bold",size = 20, hjust = 0.5)
  plot_grid(title, meta_plot,ncol = 1, rel_heights = c(0.1, 1))
  
  ggsave(filename = here('Processed_Data/Meta_Analysis_Results/Heatmaps/top_genes_meta_p_heatmap.png'), dpi=300, width=8, height=8)
  return(meta_plot)
}