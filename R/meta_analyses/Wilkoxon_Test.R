library(dplyr)
library(here)
library(readr)
library(magrittr)

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
corticalmagma_wilcox<-runWilcox(corticalmagmatable, howard_genes) #population has a median distinct from the hypothetical value you entered
corticalmagmarank_wilcox<-runWilcox(corticalmagmatablerank, howard_genes)
SIcorticalmagma_wilcox<-runWilcox(SIcorticalmagmatable, howard_genes)
SIcorticalmagmarank_wilcox<-runWilcox(SIcorticalmagmatablerank, howard_genes)

