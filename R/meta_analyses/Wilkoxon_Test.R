library(dplyr)
library(here)
library(readr)
library(magrittr)

runWilcox <- function(table,howard_genes){
  
  #Label the genes as Howard or Genome 
  magma_table<- table %>% mutate(group = ifelse(gene_symbol %in% howard_genes, "Howard", "Genome")) %>% select(meta_p, group)
  
  #visualize the meta p values in a box plot (no visible difference)
  boxplot(meta_p ~ group, data = magma_table)
  
  #compare the howard_table meta p values of genome vs howard using wilcox rank sum test
  return(wilcox.test(meta_p ~ group, alternative = "two.sided", paired = FALSE, mu = 0, data = magma_table))
  #null hypothesis is that there is no difference in meta p values between the genome genes and howard genes
  #p-value = 0.9985 a Wilcoxon rank sum test (equivalent to the Mann-Whitney test) 
  #meta p values are not significantly different Howard vs rest of genome, accept the null hypothesis 
}

#get the meta_p values for howard 269 
howard_genes <- read_csv(here("Processed_Data/HowardEtAl/fullHowardTable.csv")) %>% select(gene_symbol) %>% pull()

#get the meta_p values for the rest of the magma genome -- get from running mergeMetaOnSlims.R
fullmagmatable <- read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/FullTableMetaMagma.csv"))
fullmagmatablerank <-read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/FullTableMetaMagmaRank.csv"))
fullmagmatablerank %<>% rename(meta_p = meta_empirical_p)
SIfullmagmatable <-read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/SIFullTableMetaMagma.csv"))
SIfullmagmatablerank <-read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/SIFullTableMetaMagmaRank.csv"))
SIfullmagmatablerank %<>% rename(meta_p = meta_empirical_p)
feamlemagmatable <-read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/FemaleTableMetaMagma.csv"))
femalemagmatablerank <-read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/FemaleTableMetaMagmaRank.csv"))
femalemagmatablerank %<>% rename(meta_p = meta_empirical_p)
malemagmatable <-read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/MaleTableMetaMagma.csv"))
malemagmatablerank <-read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/MaleTableMetaMagmaRank.csv"))
malemagmatablerank %<>% rename(meta_p = meta_empirical_p)
corticalmagmatable <-read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/CorticalTableMetaMagma.csv"))
corticalmagmatablerank <-read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/CorticalTableMetaMagmaRank.csv"))
corticalmagmatablerank %<>% rename(meta_p = meta_empirical_p)
SIcorticalmagmatable <-read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/SICorticalTableMetaMagma.csv"))
SIcorticalmagmatablerank <-read_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/SICorticalTableMetaMagmaRank.csv"))
SIcorticalmagmatablerank %<>% rename(meta_p = meta_empirical_p)

fullmagma_wilcox<- runWilcox(fullmagmatable, howard_genes)
fullmagmarank_wilcox<-runWilcox(fullmagmatablerank, howard_genes)
SIfullmagma_wilcox<-runWilcox(SIfullmagmatable, howard_genes)
SIfullmagmarank_wilcox<-runWilcox(SIfullmagmatablerank, howard_genes)
femalemagma_wilcox<-runWilcox(feamlemagmatable, howard_genes)
femalemagmarank_wilcox<-runWilcox(femalemagmatablerank, howard_genes)
malemagma_wilcox<-runWilcox(malemagmatable, howard_genes)
malemagmarank_wilcox<-runWilcox(malemagmatablerank, howard_genes)
corticalmagma_wilcox<-runWilcox(corticalmagmatable, howard_genes) #population has a median distinct from the hypothetical value you entered
corticalmagmarank_wilcox<-runWilcox(corticalmagmatablerank, howard_genes)
SIcorticalmagma_wilcox<-runWilcox(SIcorticalmagmatable, howard_genes)
SIcorticalmagmarank_wilcox<-runWilcox(SIcorticalmagmatablerank, howard_genes)
# 
# fullMetaTable <- read_sheet(drive_get("Meta-Analysis"), sheet = 'Full_Meta_Analysis')
# femaleMetaTable<-read_sheet(drive_get("Meta-Analysis"), sheet = 'Female_Meta_Analysis')
# maleMetaTable <-read_sheet(drive_get("Meta-Analysis"), sheet = 'Male_Meta_Analysis')
# corticalMetaTable <- read_sheet(drive_get("Meta-Analysis"), sheet = 'Cortical_Meta_Analysis')
# 
# fullMetaTableRank <- read_sheet(drive_get("Meta-Analysis-Rank"), sheet = 'Full_Meta_Analysis_Rank')
# femaleMetaTableRank<-read_sheet(drive_get("Meta-Analysis-Rank"), sheet = 'Female_Meta_Analysis_Rank')
# maleMetaTableRank <-read_sheet(drive_get("Meta-Analysis-Rank"), sheet = 'Male_Meta_Analysis_Rank')
# corticalMetaTableRank <- read_sheet(drive_get("Meta-Analysis-Rank"), sheet = 'Cortical_Meta_Analysis_Rank')
# 
# 
# SIfullMetaTable <- read_sheet(drive_get("Sex-Interaction-Meta-Analysis"), sheet = 'Full_Meta_Analysis')
# SIcorticalMetaTable <- read_sheet(drive_get("Sex-Interaction-Meta-Analysis"), sheet = 'Cortical_Meta_Analysis')
# 
# SIfullMetaTableRank <- read_sheet(drive_get("Sex-Interaction-Meta-Analysis-Rank"), sheet = 'Full_Meta_Analysis_Rank')
# SIcorticalMetaTableRank <- read_sheet(drive_get("Sex-Interaction-Meta-Analysis-Rank"), sheet = 'Cortical_Meta_Analysis_Rank')
