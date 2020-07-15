library(gtools)
library(dplyr)
library(here)
library(readr)
library(magrittr)
library(tidyverse)

magma_genes <-  read_csv(here("Raw_Data/HowardEtAl/FullMagmaGenes.csv")) 
labonte_magma <- magma_genes %>% select(Labonte_genes) %>% distinct() %>% na.omit()
ding_magma <- magma_genes %>% select(Ding_genes) %>% distinct() %>% na.omit()
ramaker_magma <- magma_genes %>% select(Ramaker_genes) %>% distinct() %>% na.omit()


Howard <- read_csv(here("Processed_Data/HowardEtAl/fullHowardTable.csv"))

Labonte_table <- read_csv(here("Processed_Data/LabonteEtAl/CompleteLabonteTable.csv"))
Labonte_cortical <- Labonte_table %>% filter(brain_region != "Nac") %>% filter( brain_region != "Subic") %>% dplyr::select(brain_region)%>% distinct() %>%pull()
Ramaker_table <- read_csv(here("Processed_Data/RamakerEtAl/CompleteRamakerTable.csv")) 
Ramaker_cortical <- Ramaker_table %>% filter(target_region != "nAcc") %>% dplyr::select(target_region) %>% distinct() %>% pull()
Ding_table <- read_csv(here("Processed_Data/DingEtAl/CompleteDingTable.csv"))  
Ding_cortical <- Ding_table %>% filter(brain_region != "AMY") %>% dplyr::select(brain_region) %>% distinct() %>% pull()

l_combinations <- combinations(length(Labonte_cortical),2,Labonte_cortical) %>% as_tibble() #all possible combinations of 2
r_combinations <- combinations(length(Ramaker_cortical),1,Ramaker_cortical) %>% as_tibble()#all possible combinations of 1
d_combinations <- combinations(length(Ding_cortical),1,Ding_cortical) %>% as_tibble() #all possible combinations of 1

combined_df <- cbind(l_combinations, r_combinations, d_combinations) 
new_cols <- c('l_1', 'l_2', 'r', 'd')
colnames(combined_df) <- new_cols
combined_df %<>% as_tibble()

#make all possible combinations to match subcortical analysis
combined_df %<>% complete(nesting(l_1,l_2), nesting(r,d))

dir.create(here('Results/Meta_Analysis_Results/cortical_combinations'), recursive = TRUE)

source(here("R/transcriptomic_meta/Labonte_Meta_Analysis.R"))
source(here("R/transcriptomic_meta/Ramaker_Meta_Analysis.R"))
source(here("R/transcriptomic_meta/Ding_Meta_Analysis.R"))
source(here("R/transcriptomic_meta/Percentile_Rank_Analysis.R"))
source(here("R/meta_analyses/mergeAnalysis.R"))

#re-do cortical analyses with the defined brain regions 
#number of p-values we need to filter = 4 (2 regions * 2 sexes)
for (i in 1:nrow(combined_df)){
  file_name <- paste0('cortical_combinations_run_',i,'.csv')
  
  labonte_region_1 <- nth(combined_df,1)[i]
  labonte_region_2 <- nth(combined_df,2)[i]
  order <- sort(c(labonte_region_1, labonte_region_2))
  labonte_dir <-paste0('female_',order[1],'_female_',order[2],'_male_',order[1],'_male_',order[2])
  labonte_summary_table <- Labonte_table %>% filter(brain_region == labonte_region_1 | brain_region == labonte_region_2) 
  labonte_summary_table %<>% LabonteMetaAnalysis(4,2)
  
  labonte_summary_table%<>% right_join(labonte_magma %>% select(Labonte_genes) %>% distinct() %>% na.omit(), by = c('gene_symbol' = 'Labonte_genes'))
  labonte_summary_table %<>% getRank()
  # labonte_summary_table %>% write_csv(here("Results/LabonteEtAl/CorticalLabonteTableMagma.csv"))
  
  
  ramaker_region <- nth(combined_df,3)[i]
  ramaker_colname_male <- paste0(ramaker_region,'_Male_directions')
  ramaker_colname_female <- paste0(ramaker_region,'_Female_directions')
  ramaker_colname <- paste0(ramaker_region,'_directions')
  ramaker_united_colname <- paste0(ramaker_region,'.F_',ramaker_region,'.M')
  cortical_ramaker <- Ramaker_table %>% filter(target_region == ramaker_region) %>% RamakerMetaAnalysis(ramaker_region)
  #Extract the cortical region data & run analysis on female data
  female_ramaker_cortical <- read_csv(here("Processed_Data/RamakerEtAl/CompleteFemaleRamakerTable.csv"))
  female_ramaker_cortical %<>% filter(target_region == ramaker_region)
  #run the meta-analysis on only the cortical regions sampled in female data
  female_ramaker_cortical %<>% RamakerMetaAnalysis(ramaker_region)
  female_ramaker_cortical %<>% rename(!!ramaker_colname_female := ramaker_colname)
  # female_ramaker_cortical %>% write_csv(here("Processed_Data/RamakerEtAl/CorticalFemaleRamakerTable.csv"))
  
  #Extract the cortical region data & run analysis on male data
  male_ramaker_cortical<- read_csv(here("Processed_Data/RamakerEtAl/CompleteMaleRamakerTable.csv"))
  male_ramaker_cortical %<>% filter(target_region == ramaker_region)
  #run the meta-analysis on only the cortical regions sampled in male data
  male_ramaker_cortical %<>% RamakerMetaAnalysis(ramaker_region)
  male_ramaker_cortical %<>% rename(!!ramaker_colname_male := !!ramaker_colname)
  # male_ramaker_cortical %>% write_csv(here("Processed_Data/RamakerEtAl/CorticalMaleRamakerTable.csv"))
  
  #merge all directions from male and female data to visualize cortical directions across sexes
  ramaker_summary_table <- left_join(female_ramaker_cortical %>% select(gene_symbol,!!ramaker_colname_female), male_ramaker_cortical %>% select(gene_symbol,!!ramaker_colname_male ))
  ramaker_summary_table %<>% unite(col = !!ramaker_united_colname, !!ramaker_colname_female, !!ramaker_colname_male, sep = "")
  ramaker_summary_table %<>% left_join(cortical_ramaker%>% select(-sex, -!!ramaker_colname) %>% distinct())
  
  ramaker_summary_table %<>% right_join(ramaker_magma, by = c('gene_symbol' = 'Ramaker_genes'))
  ramaker_summary_table %<>% getRank()
  # ramaker_summary_table %>% write_csv(here("Processed_Data/RamakerEtAl/CorticalRamakerTableMagma.csv"))
  
  
  ding_region <- nth(combined_df,4)[i]
  ding_dir <-paste0('female_',ding_region,'_male_',ding_region)
  #Perform meta-analysis on cortical brain regions in both sexes
  #number of p-values we need to filter = 2 (1 region * 2 sex) 
  #number of brain regions = 1
  ding_summary_table <- Ding_table %>% filter(brain_region == ding_region) %>% DingMetaAnalysis(2,1)
  
  ding_summary_table%<>% right_join(ding_magma, by = c('gene_symbol' = 'Ding_genes'))
  ding_summary_table %<>% getRank()
  
  
  corticalTable <- mergeMetaStudies(Howard, labonte_summary_table, labonte_dir, ding_summary_table, ding_dir, ramaker_summary_table, ramaker_united_colname)
  corticalTable %<>% MetaAnalysis()
  
  corticalTable %>% write_csv(here('Results/Meta_Analysis_Results/cortical_combinations',file_name))
}
 
#import the meta-analysis files
 
 
 
 
 
 
 