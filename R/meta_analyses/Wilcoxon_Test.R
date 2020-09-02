library(dplyr)
library(here)
library(readr)
library(magrittr)



############################################################################################################################################################
################ WILKOXON RANK SUM TEST ON THE SEXINTERACTION ANALYSES VS FULL & CORTICAL META-ANALYSES plot top genes ############################

full <- read_csv(here('Results', 'Tables', 'Meta_Analysis', "Full_Meta_Analysis.csv"))
cortical <- read_csv(here('Results', 'Tables', 'Meta_Analysis', 'Cortical_Meta_Analysis.csv'))
SI_Full <- read_csv(here('Results', 'Tables', 'Meta_Analysis','Sex_interaction_Full_Meta_Analysis.csv'))
SI_Cortical <- read_csv(here('Results', 'Tables', 'Meta_Analysis', 'Sex_interaction_Cortical_Meta_Analysis.csv'))

#compare full with full
full_wilkoxon <- full %>% filter(Corrected_p < 0.05) %>% mutate(group = "meta") %>% dplyr::select(group, Corrected_p)
SI_Full_wilkoxon <- SI_Full %>% filter(Corrected_p < 0.05) %>% mutate(group = "SI") %>% dplyr::select(group, Corrected_p)

table_full <- rbind(full_wilkoxon, SI_Full_wilkoxon)

boxplot(Corrected_p ~ group, data = table_full)


#compare cortical with cortical 
cortical_wilkoxon <- cortical %>% filter(Corrected_p < 0.05) %>% mutate(group = "meta") %>% dplyr::select(group, Corrected_p)
SI_cortical_wilkoxon <- SI_Cortical %>% filter(Corrected_p < 0.05) %>% mutate(group = "SI") %>% dplyr::select(group, Corrected_p)

table_cortical <- rbind(cortical_wilkoxon, SI_cortical_wilkoxon)

boxplot(Corrected_p ~ group, data = table_cortical)


#compare SI will meta
combined_table <- rbind(table_full,table_cortical)
boxplot(Corrected_p ~ group, data = combined_table)

########################################################################################################################################################
###################################################################### PAIRED 269 genes ################################################################
full <- read_csv(here('Results', 'Tables', 'Meta_Analysis', "Full_Meta_Analysis.csv"))
cortical <- read_csv(here('Results', 'Tables', 'Meta_Analysis', 'Cortical_Meta_Analysis.csv'))
SI_Full <- read_csv(here('Results', 'Tables', 'Meta_Analysis','Sex_interaction_Full_Meta_Analysis.csv'))
SI_Cortical <- read_csv(here('Results', 'Tables', 'Meta_Analysis', 'Sex_interaction_Cortical_Meta_Analysis.csv'))

full %<>% dplyr::select(gene_symbol, meta_p)
cortical %<>% dplyr::select(gene_symbol, meta_p)
SI_Full %<>% dplyr::select(gene_symbol, meta_p)
SI_Cortical %<>% dplyr::select(gene_symbol, meta_p)


full_analysis <- full %>% inner_join(SI_Full, by = c('gene_symbol' = 'gene_symbol'))
wilcox.test(full_analysis$meta_p.x, full_analysis$meta_p.y, paired=TRUE)

plot <- full %>% mutate(group = "meta") %>% dplyr::select(-gene_symbol)
full_table <- rbind(plot, SI_Full %>% mutate(group = "SI")%>% dplyr::select(-gene_symbol))
boxplot(meta_p ~ group, data = full_table)

cortical_analysis <- cortical %>% inner_join(SI_Cortical, by = c('gene_symbol' = 'gene_symbol'))
wilcox.test(cortical_analysis$meta_p.x, cortical_analysis$meta_p.y, paired=TRUE)
plot_c <- cortical %>% mutate(group = "meta") %>% dplyr::select(-gene_symbol)
cortical_table <- rbind(plot_c, SI_Cortical %>% mutate(group = "SI")%>% dplyr::select(-gene_symbol))
boxplot(meta_p ~ group, data = cortical_table)


############################################################################################################################################################
################ WILKOXON RANK SUM TEST ON THE 269 VS MAGMA UNIVERSE ############################

runWilcox <- function(table,howard_genes){
  
  #Label the genes as Howard or Genome 
  magma_table<- table %>% mutate(group = ifelse(gene_symbol %in% howard_genes, "Howard", "Genome")) %>% dplyr::select(meta_p, group)
  
  #visualize the meta p values in a box plot (no visible difference)
  boxplot(meta_p ~ group, data = magma_table)
  
  #compare the howard_table meta p values of genome vs howard using wilcox rank sum test
  return(wilcox.test(meta_p ~ group, alternative = "two.sided", paired = FALSE, mu = 0, data = magma_table))
  #null hypothesis is that there is no difference in meta p values between the genome genes and howard genes
  #p-value = 0.9985 a Wilcoxon rank sum test (equivalent to the Mann-Whitney test) 
  #meta p values are not significantly different Howard vs rest of genome, accept the null hypothesis 
}

#get the meta_p values for howard 269 
howard_genes <- read_csv(here("Processed_Data/HowardEtAl/fullHowardTable.csv")) %>% dplyr::select(gene_symbol) %>% pull()

#get the meta_p values for the rest of the magma genome -- get from running mergeMetaOnSlims.R
fullmagmatable <- read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/FullTableMetaMagma.csv"))
fullmagmatablerank <-read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/FullTableMetaMagmaRank.csv"))
fullmagmatablerank %<>% dplyr::rename(meta_p = meta_empirical_p)
SIfullmagmatable <-read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/SIFullTableMetaMagma.csv"))
SIfullmagmatablerank <-read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/SIFullTableMetaMagmaRank.csv"))
SIfullmagmatablerank %<>% dplyr::rename(meta_p = meta_empirical_p)
feamlemagmatable <-read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/FemaleTableMetaMagma.csv"))
femalemagmatablerank <-read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/FemaleTableMetaMagmaRank.csv"))
femalemagmatablerank %<>% dplyr::rename(meta_p = meta_empirical_p)
malemagmatable <-read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/MaleTableMetaMagma.csv"))
malemagmatablerank <-read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/MaleTableMetaMagmaRank.csv"))
malemagmatablerank %<>% dplyr::rename(meta_p = meta_empirical_p)
corticalmagmatable <-read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/CorticalTableMetaMagma.csv"))
corticalmagmatablerank <-read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/CorticalTableMetaMagmaRank.csv"))
corticalmagmatablerank %<>% dplyr::rename(meta_p = meta_empirical_p)
subcorticalmagmatable <-read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/SubcorticalTableMetaMagma.csv"))
subcorticalmagmatablerank <-read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/SubcorticalTableMetaMagmaRank.csv"))
subcorticalmagmatablerank %<>% dplyr::rename(meta_p = meta_empirical_p)
SIcorticalmagmatable <-read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/SICorticalTableMetaMagma.csv"))
SIcorticalmagmatablerank <-read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/SICorticalTableMetaMagmaRank.csv"))
SIcorticalmagmatablerank %<>% dplyr::rename(meta_p = meta_empirical_p)
SIsubcorticalmagmatable <-read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/SISubcorticalTableMetaMagma.csv"))
SIsubcorticalmagmatablerank <-read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/SISubcorticalTableMetaMagmaRank.csv"))
SIsubcorticalmagmatablerank %<>% dplyr::rename(meta_p = meta_empirical_p)


fullmagma_wilcox<- runWilcox(fullmagmatable, howard_genes)
fullmagmarank_wilcox<-runWilcox(fullmagmatablerank, howard_genes)
SIfullmagma_wilcox<-runWilcox(SIfullmagmatable, howard_genes)
SIfullmagmarank_wilcox<-runWilcox(SIfullmagmatablerank, howard_genes)
femalemagma_wilcox<-runWilcox(feamlemagmatable, howard_genes)
femalemagmarank_wilcox<-runWilcox(femalemagmatablerank, howard_genes)
malemagma_wilcox<-runWilcox(malemagmatable, howard_genes)
malemagmarank_wilcox<-runWilcox(malemagmatablerank, howard_genes)
corticalmagma_wilcox<-runWilcox(corticalmagmatable, howard_genes) 
corticalmagmarank_wilcox<-runWilcox(corticalmagmatablerank, howard_genes)
SIcorticalmagma_wilcox<-runWilcox(SIcorticalmagmatable, howard_genes)
SIcorticalmagmarank_wilcox<-runWilcox(SIcorticalmagmatablerank, howard_genes)

