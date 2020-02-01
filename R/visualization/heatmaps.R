library(cowplot)
library(here)
library(readr)
library(readxl)
library(homologene)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(ggplot2)
library(googlesheets4)
library(googledrive)

top_genes <- c("MANEA","UBE2M","CKB","ITPR3","SPRY2","SAMD5","TMEM106B","ZC3H7B","LST1","ASXL3","HSPA1A") 

full <- read_sheet(drive_get("Meta-Analysis"), sheet = 'Full_Meta_Analysis') 
female<-read_sheet(drive_get("Meta-Analysis"), sheet = 'Female_Meta_Analysis')
male <-read_sheet(drive_get("Meta-Analysis"), sheet = 'Male_Meta_Analysis')
cortical <- read_sheet(drive_get("Meta-Analysis"), sheet = 'Cortical_Meta_Analysis')
SIfull <- read_sheet(drive_get("Sex-Interaction-Meta-Analysis"), sheet = 'Full_Meta_Analysis')
SIcortical <- read_sheet(drive_get("Sex-Interaction-Meta-Analysis"), sheet = 'Cortical_Meta_Analysis')

fulltable <- full %>% select(gene_symbol, Corrected_p) %>% rename(Full = Corrected_p)
femaletable <-female %>% select(gene_symbol, Corrected_p)  %>% rename(Female = Corrected_p)
maletable <- male %>% select(gene_symbol, Corrected_p) %>% rename(Male = Corrected_p)
corticaltable <- cortical %>% select(gene_symbol, Corrected_p) %>% rename(Cortical = Corrected_p)
SIfulltable <- SIfull%>% select(gene_symbol, Corrected_p) %>% rename(Sex_Interaction_Full = Corrected_p)
SIcorticaltable <- SIcortical%>% select(gene_symbol, Corrected_p) %>% rename(Sex_Interaction_Cortical = Corrected_p)

full_meta_analysis_results <- left_join(fulltable, femaletable) %>% left_join(corticaltable) %>% left_join(maletable) %>% left_join(SIfulltable) %>% left_join(SIcorticaltable)
full_meta_analysis_results %<>% gather("Analysis", "meta_p", 2:ncol(full_meta_analysis_results))
full_meta_analysis_results %<>% filter(gene_symbol %in% top_genes)
full_meta_analysis_results %<>% mutate(meta_p = as.numeric(meta_p))

# test <- full_meta_analysis_results %>% mutate(meta_p = if_else(meta_p > 0.05, 0.05, meta_p))

full_meta_analysis_results %<>% mutate(gene_symbol = factor(gene_symbol, levels=rev(c("HSPA1A","ZC3H7B", "SAMD5", "SPRY2", "ITPR3", "MANEA", "UBE2M", "CKB", "TMEM106B","ASXL3", "LST1"))),
                 Analysis = factor(Analysis, levels = c("Full", "Cortical", "Male", "Female", "Sex_Interaction_Cortical", "Sex_Interaction_Full")))
full_meta_analysis_results <- cut(full_meta_analysis_results$meta_p, breaks = c(0, 0.0005, 0.001,0.005,0.01, 0.05, 0.1, 0.5, 5, 20,220),
               labels = c("0.-0.0005", "0.0005-0.001","0.001-0.005","0.005-0.01", "0.01-0.05", "0.05 - 0.1", "0.1-0.5", "0.5-5", "5-20","20-220"),
               include.lowest = TRUE)
#arrange these genes with significant ones together 
analysis_type <- full_meta_analysis_results$Analysis
symbol <- full_meta_analysis_results$gene_symbol
p_val <- full_meta_analysis_results$meta_p

#heatmap
ggplot(full_meta_analysis_results, aes(x=analysis_type,y=symbol, fill = p_val)) + 
  geom_tile(aes(x=analysis_name, fill = Y1), colour='black') +
  labs(x="Meta Analysis", fill = "Corrected\np-value") +
  scale_fill_brewer(palette = "RdBu") +
  # scale_fill_gradient2(low="green",mid = "white", high="blue", guide="colorbar")+
  # scale_fill_distiller(palette = "RdBu")+
  ylab('') +
  theme(axis.text.x = element_text(angle = 45, hjust=1,vjust=1))


# plot the direction of expressions
fullLabonte <- read_csv( here("ProcessedData", "LabonteEtAl", "fullLabonteTable_magma.csv")) 
femaleLabonte <- read_csv(here("ProcessedData", "LabonteEtAl", "FemaleLabonteTable_magma.csv")) %>% select(gene_symbol, meta_direction)
maleLabonte <- read_csv(here("ProcessedData", "LabonteEtAl", "MaleLabonteTable_magma.csv"))%>% select(gene_symbol, meta_direction)
corticalLabonte <- read_csv( here("ProcessedData", "LabonteEtAl", "CorticalLabonteTableMagma.csv"))%>% select(gene_symbol, meta_direction)
fullDing <- read_csv(here("ProcessedData", "DingEtAl", "DingTableFisher.slim_magma.csv"))%>% select(gene_symbol, meta_direction)
femaleDing <-read_csv(here("ProcessedData", "DingEtAl", "DingTableFisher.slim.Female_magma.csv"))%>% select(gene_symbol, meta_direction)
maleDing <-read_csv(here("ProcessedData", "DingEtAl", "DingTableFisher.slim.Male_magma.csv"))%>% select(gene_symbol, meta_direction)
corticalDing <- read_csv( here("ProcessedData", "DingEtAl", "CorticalDingTableFisherMagma.slim.csv"))%>% select(gene_symbol, meta_direction)
fullRamaker <- read_csv(here("ProcessedData", "RamakerEtAl", "fullRamakerTable_magma.csv"))%>% select(gene_symbol, meta_direction)
femaleRamaker <- read_csv(here("ProcessedData", "RamakerEtAl", "FemaleRamakerTable_magma.csv"))%>% select(gene_symbol, meta_direction)
maleRamaker <- read_csv(here("ProcessedData", "RamakerEtAl", "MaleRamakerTable_magma.csv"))%>% select(gene_symbol, meta_direction)
corticalRamaker <- read_csv(here("ProcessedData", "RamakerEtAl", "fullCorticalRamakerTable.erTable.csv"))%>% select(gene_symbol, meta_direction)
Labonte <- read_csv(here("ProcessedData", "LabonteEtAl","CompleteLabonteTable.csv"))
Labonte %<>% select(symbol, brain_region, Male.logFC, Female.logFC)
Labonte %<>% filter(symbol %in% top_genes)

Labonte_male <- Labonte %>% select(-Female.logFC)
Labonte_female <- Labonte %>% select(-Male.logFC)
#arrange these genes with significant ones together 
Labonte %<>% mutate(symbol = factor(symbol, levels=c("HSPA1A","MANEA","TMEM106B","SPRY2","SAMD5","ASXL3","ZC3H7B","UBE2M", "CKB" ,"LST1" ,"ITPR3")))

Labonte_female %<>% mutate(symbol = factor(symbol, levels=rev(c("ZC3H7B","SAMD5","ASXL3","ITPR3","CKB","LST1","UBE2M","SPRY2","MANEA", "TMEM106B" ,"HSPA1A"))))

regions <- Labonte$brain_region
symbol <- Labonte$symbol
male_symbol <- Labonte_male$symbol
male_direction <- Labonte$Male.logFC
female_symbol <- Labonte_female$symbol
female_direction <- Labonte$Female.logFC

#heatmap
male_labonte <- ggplot(Labonte, aes(x=regions,y=symbol, fill = male_direction)) + 
  geom_tile(colour='black') +
  labs(x="Male", fill = "Expression\nlogFC") +
  scale_fill_distiller(palette = "RdBu") +
  # scale_fill_gradient2(low="green",mid = "white", high="blue", guide="colorbar")+
  # scale_fill_distiller(palette = "RdBu")+
  ylab('') +
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust=1,vjust=1))

#heatmap
female_labonte <- ggplot(Labonte, aes(x=regions,y=symbol, fill = female_direction)) + 
  geom_tile(colour='black') +
  labs(x=NULL, fill = "Expression\nlogFC") +
  scale_fill_distiller(palette = "RdBu") +
  ylab('') +
  #Themes are a powerful way to customize the non-data components of your plots: i.e. titles, labels, fonts, background, gridlines, and legends.
  # change the angle of the x-axis text to 45 for readability 
  # position the x-axis text so it doesn't go into the heatmap
  theme(axis.text.x = element_text(angle = 45, hjust=1,vjust=1))

labonte_plots <- plot_grid(male_labonte, female_labonte)

title <- ggdraw() + draw_label("Labonte", fontface = "bold",size = 28,x = 0, hjust = -3.7)

plot_grid(title, labonte_plots,ncol = 1, rel_heights = c(0.1, 1)) 
  
clustHeatHuman <- 
  ggplot(forClustHeatmap, aes(x=sample_origin, y=reorder(human_gene_symbol, rank_order))) +
  facet_wrap(~species, scales='free_x') +#, space="free_x") + 
  geom_tile(aes(x=sample_origin, fill = z_expression), colour='black') +
  labs(x=NULL, fill = "Expression\n(z-score)") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(position = 'right') +
  theme(axis.text.y  = element_text(size=8)) +
  ylab('') +
  scale_fill_distiller(palette = 'RdBu') +
  theme(axis.text.x = element_text(angle = 90, hjust=1,vjust=0.5), 
        plot.margin = grid::unit(c(t=0,r=0,b=2, l=0), "mm"),
        strip.text = element_blank())

clustHeatHuman
ggsave(filename = here('results', 'mouse', "DAM_dataset", 'clustered_heatmap_human.png'), dpi=300, width=4, height=12)