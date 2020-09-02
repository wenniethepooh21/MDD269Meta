library(ggplot2)
library(googledrive)
library(googlesheets4)
library(dplyr)
library(tidyr)
library(viridis)
library(cowplot)
#This script creates a heatmap to display the corrected meta p-values of our top 12 genes and how they ranked in each meta-analysis

drawMetaHeat <- function() {
  top_genes <- c("MANEA","UBE2M","CKB","ITPR3","SPRY2","SAMD5","TMEM106B","ZC3H7B","LST1","ASXL3","HSPA1A","ZNF184") 

  full <- read_csv(here('Results', 'Tables', 'Meta_Analysis', "Full_Meta_Analysis.csv"))
  female <- read_csv(here('Results', 'Tables', 'Meta_Analysis', 'Female_Meta_Analysis.csv'))
  male <- read_csv(here('Results', 'Tables', 'Meta_Analysis', 'Male_Meta_Analysis.csv'))
  cortical <- read_csv(here('Results', 'Tables', 'Meta_Analysis', 'Cortical_Meta_Analysis.csv'))
  subcortical <- read_csv(here('Results', 'Tables', 'Meta_Analysis', 'Subcortical_Meta_Analysis.csv'))
  SIfull <- read_csv(here('Results', 'Tables', 'Meta_Analysis','Sex_interaction_Full_Meta_Analysis.csv'))
  SIcortical <- read_csv(here('Results', 'Tables', 'Meta_Analysis', 'Sex_interaction_Cortical_Meta_Analysis.csv'))
  SIsubcortical <- read_csv(here('Results', 'Tables', 'Meta_Analysis', 'Sex_interaction_Subcortical_Meta_Analysis.csv'))
  
  fulltable <- full %>% dplyr::select(gene_symbol, Bonferroni_meta_p) %>% dplyr::rename(Full = Bonferroni_meta_p)
  femaletable <-female %>% dplyr::select(gene_symbol, Bonferroni_meta_p)  %>% dplyr::rename(Female = Bonferroni_meta_p)
  maletable <- male %>% dplyr::select(gene_symbol, Bonferroni_meta_p) %>% dplyr::rename(Male = Bonferroni_meta_p)
  corticaltable <- cortical %>% dplyr::select(gene_symbol, Bonferroni_meta_p) %>% dplyr::rename(Cortical = Bonferroni_meta_p)
  subcorticaltable <- subcortical %>% dplyr::select(gene_symbol, Bonferroni_meta_p) %>% dplyr::rename(Subcortical = Bonferroni_meta_p)
  SIfulltable <- SIfull%>% dplyr::select(gene_symbol, Bonferroni_meta_p) %>% dplyr::rename(SI_Full = Bonferroni_meta_p)
  SIcorticaltable <- SIcortical%>% dplyr::select(gene_symbol, Bonferroni_meta_p) %>% dplyr::rename(SI_Cortical = Bonferroni_meta_p)
  SIsubcorticaltable <- SIsubcortical%>% dplyr::select(gene_symbol, Bonferroni_meta_p) %>% dplyr::rename(SI_Subcortical = Bonferroni_meta_p)
  
  #join the results into one table
  full_meta_analysis_results <- left_join(fulltable, femaletable) %>% left_join(corticaltable) %>% left_join(maletable) %>% left_join(SIfulltable) %>% left_join(SIcorticaltable) %>% left_join(subcorticaltable) %>% left_join(SIsubcorticaltable)
  full_meta_analysis_results %<>% gather("Analysis", "Bonferroni_meta_p", 2:ncol(full_meta_analysis_results))
  full_meta_analysis_results %<>% filter(gene_symbol %in% top_genes)
  full_meta_analysis_results %<>% mutate(Bonferroni_meta_p = as.numeric(Bonferroni_meta_p))
  

  full_meta_analysis_results %<>% mutate(gene_symbol = factor(gene_symbol, levels=rev(c("HSPA1A","ZC3H7B", "ITPR3", "UBE2M", "CKB", "SAMD5", "SPRY2", "TMEM106B", "LST1","ASXL3", "MANEA","ZNF184"))),
                                         Analysis = factor(Analysis, levels = c("Full", "Cortical","Subcortical", "Male", "Female", "SI_Cortical","SI_Subcortical", "SI_Full")))
  #set the scale to range from 0 to 0.1
  full_meta_analysis_results %<>% mutate(Bonferroni_meta_p = Bonferroni_meta_p/10)
 
  
  #arrange these genes with significant ones together 
  analysis_type <- full_meta_analysis_results$Analysis
  symbol <- full_meta_analysis_results$gene_symbol
  p_val <- full_meta_analysis_results$Bonferroni_meta_p
  
  #heatmap
  meta_plot <- ggplot(full_meta_analysis_results, aes(x=analysis_type,y=symbol)) + 
    geom_tile(aes(fill = p_val), colour='black') +
    labs(x="Meta-Analysis", fill = "Corrected\np-value") +
    ylab('') +
    scale_fill_viridis(option = "D", direction = -1) + 
    theme(axis.title.x = element_text(size = 11),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 8,angle = 20, hjust=0.8,vjust=0.9), 
          axis.text.y = element_text(face = "italic",size = 10), 
          plot.title = element_text(size = 20,face = "bold",hjust = 0.5),
          panel.background = element_blank(),
          plot.margin = unit(c(t=0.5,r=1,l=1,b=2),"cm"),
          legend.margin = margin(unit(c(t=0,r=-10,l=-0.5,b = 0), "cm"))) # center the title 

  
  title <- ggdraw() + draw_label("Meta p-values of Top 12 Genes", fontface = "bold",size = 20, hjust = 0.4,vjust = 0.5)
  full_meta_plot <- plot_grid(title, meta_plot,ncol = 1, rel_heights = c(0.1, 1), rel_widths = c(1,0.5))
  
  ggsave(filename = here('Results/Heatmaps/top_genes_meta_p_heatmap.png'), dpi=300, width=8, height=8)
  return(full_meta_plot)
}
