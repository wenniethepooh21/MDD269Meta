library(readr)
library(magrittr)
library(dplyr)
library(tidyr)
library(here)

fullLabonte <- read_csv(here("Processed_Data/LabonteEtAl/FullLabonteTableMagma.csv"))
femaleLabonte <- read_csv(here("Processed_Data/LabonteEtAl/FemaleLabonteTableMagma.csv")) 
maleLabonte <- read_csv(here("Processed_Data/LabonteEtAl/MaleLabonteTableMagma.csv"))
corticalLabonte <- read_csv( here("Processed_Data/LabonteEtAl/CorticalLabonteTableMagma.csv"))
subcorticalLabonte <- read_csv( here("Processed_Data/LabonteEtAl/SubcorticalLabonteTableMagma.csv"))
fullLabonte_flipped <- read_csv(here("Processed_Data/LabonteEtAl/FullLabonteTableMagma_flipped.csv"))
corticalLabonte_flipped <- read_csv( here("Processed_Data/LabonteEtAl/CorticalLabonteTableMagma_flipped.csv"))
subcorticalLabonte_flipped <- read_csv( here("Processed_Data/LabonteEtAl/SubcorticalLabonteTableMagma_flipped.csv"))

fullDing <- read_csv(here("Processed_Data/DingEtAl/FullDingTableMagma.csv"))
femaleDing <-read_csv(here("Processed_Data/DingEtAl/FemaleDingTableMagma.csv"))
maleDing <-read_csv(here("Processed_Data/DingEtAl/MaleDingTableMagma.csv"))
corticalDing <- read_csv( here("Processed_Data/DingEtAl/CorticalDingTableMagma.csv"))
subcorticalDing <- read_csv( here("Processed_Data/DingEtAl/SubcorticalDingTableMagma.csv"))
fullDing_flipped <- read_csv(here("Processed_Data/DingEtAl/FullDingTableMagma_flipped.csv"))
corticalDing_flipped <- read_csv(here("Processed_Data/DingEtAl/CorticalDingTableMagma_flipped.csv"))
subcorticalDing_flipped <- read_csv(here("Processed_Data/DingEtAl/SubcorticalDingTableMagma_flipped.csv"))


fullRamaker <- read_csv(here("Processed_Data/RamakerEtAl/FullRamakerTableMagma.csv"))
femaleRamaker <- read_csv(here("Processed_Data/RamakerEtAl/FemaleRamakerTableMagma.csv"))
maleRamaker <- read_csv(here("Processed_Data/RamakerEtAl/MaleRamakerTableMagma.csv"))
corticalRamaker <- read_csv(here("Processed_Data/RamakerEtAl/CorticalRamakerTableMagma.csv"))
subcorticalRamaker <- read_csv(here("Processed_Data/RamakerEtAl/SubcorticalRamakerTableMagma.csv"))
fullRamaker_flipped <- read_csv(here("Processed_Data/RamakerEtAl/FullRamakerTableMagma_flipped.csv"))
corticalRamaker_flipped <- read_csv(here("Processed_Data/RamakerEtAl/CorticalRamakerTableMagma_flipped.csv"))
subcorticalRamaker_flipped <- read_csv(here("Processed_Data/RamakerEtAl/SubcorticalRamakerTableMagma_flipped.csv"))


#read in the equivalent Howard genes used in each transcriptomic study
Howard <- read_csv(here("Processed_Data/HowardEtAl/fullHowardTable.csv"))
Howard %<>% rename("Howard.pvalue" = "Howard_pvalue")
#read in the differentially expressed "prior ranking" for the genes identified in Crow, et al. 
DE_Prior <- read_tsv(here("Raw_Data/CrowEtAl/pnas.1802973116.sd02.txt"))
DE_Prior$DE_Prior_Rank <- signif(as.numeric(DE_Prior$DE_Prior_Rank),digits=3)
#Read in the associated transcriptomic cell types and brain regions that maximally expressed each Howard gene
Howard_Polygenics <- read_csv(here("Processed_Data/HowardEtAl/HowardRegionsPolygenicCellTypesZScore_four.csv"))

# Column names - for the direction of expression in each brain region
fullLabonteDir <- "female_Anterior_Insula_female_BA11_female_BA25_female_BA8/9_female_Nac_female_Subic_male_Anterior_Insula_male_BA11_male_BA25_male_BA8/9_male_Nac_male_Subic"
femaleLabonteDir <- "female_Anterior_Insula_female_BA11_female_BA25_female_BA8/9_female_Nac_female_Subic"
maleLabonteDir <- "male_Anterior_Insula_male_BA11_male_BA25_male_BA8/9_male_Nac_male_Subic"
corticalLabonteDir <- "female_Anterior_Insula_female_BA11_female_BA25_female_BA8/9_male_Anterior_Insula_male_BA11_male_BA25_male_BA8/9"
subcorticalLabonteDir <- "female_Nac_female_Subic_male_Nac_male_Subic"

fullDingDir <- "female_ACC.1_female_ACC.2_female_AMY_female_DLPFC_male_ACC.1_male_ACC.2_male_AMY_male_DLPFC"
femaleDingDir <-"female_ACC.1_female_ACC.2_female_AMY_female_DLPFC"
maleDingDir <-"male_ACC.1_male_ACC.2_male_AMY_male_DLPFC"
corticalDingDir <-"female_ACC.1_female_ACC.2_female_DLPFC_male_ACC.1_male_ACC.2_male_DLPFC"
subcorticalDingDir <-"female_AMY_male_AMY"

fullRamakerDir <- "AnCg.F_nAcc.F_DLPFC.F_AnCg.M_nAcc.M_DLPFC.M"
femaleRamakerDir <- "AnCg_nAcc_DLPFC_Female_directions"
maleRamakerDir <- "AnCg_nAcc_DLPFC_Male_directions"
corticalRamakerDir <- "AnCg.F_DLPFC.F_AnCg.M_DLPFC.M"
subcorticalRamakerDir <- "nAcc.F_nAcc.M"

#Perform the meta-analysis combining all study-specific meta analyses together
#Functions used for combining the data from each study and meta-analysis calculations can be found here
source(here("R/meta_analyses/mergeAnalysis.R"))

##### Perform meta-analysis on the entire 17,842 MAGMA genes to run the Wilcox Test (files will be used in Wilkoxon_Test.R) ####--------------------------------------------------
magma <- read_csv(here("Processed_Data/HowardEtAl/FullMagmaGenes.csv"))
dir.create(here('Processed_Data/Meta_Analysis_Results/MAGMA'), recursive = TRUE)
# Full meta-analysis
mergeMagmaMetaRank(magma, fullLabonte, fullLabonteDir, fullDing, fullDingDir, fullRamaker, fullRamakerDir,"meta") %>% write_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/FullTableMetaMagma.csv"))
# Full meta-analysis rank
mergeMagmaMetaRank(magma, fullLabonte, fullLabonteDir, fullDing, fullDingDir, fullRamaker, fullRamakerDir,"rank") %>% write_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/FullTableMetaMagmaRank.csv"))
# Female meta-analysis
mergeMagmaMetaRank(magma, femaleLabonte, femaleLabonteDir, femaleDing, femaleDingDir, femaleRamaker, femaleRamakerDir,"meta") %>% write_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/FemaleTableMetaMagma.csv"))
# Female meta-analysis rank
mergeMagmaMetaRank(magma, femaleLabonte, femaleLabonteDir, femaleDing, femaleDingDir, femaleRamaker, femaleRamakerDir,"rank") %>% write_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/FemaleTableMetaMagmaRank.csv"))
# Male meta-analysis
mergeMagmaMetaRank(magma, maleLabonte, maleLabonteDir, maleDing, maleDingDir, maleRamaker, maleRamakerDir,"meta") %>% write_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/MaleTableMetaMagma.csv"))
# Male meta-analysis rank
mergeMagmaMetaRank(magma, maleLabonte, maleLabonteDir, maleDing, maleDingDir, maleRamaker, maleRamakerDir,"rank") %>% write_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/MaleTableMetaMagmaRank.csv"))
# Cortical meta-analysis
mergeMagmaMetaRank(magma, corticalLabonte, corticalLabonteDir, corticalDing, corticalDingDir, corticalRamaker, corticalRamakerDir,"meta") %>% write_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/CorticalTableMetaMagma.csv"))
# Cortical meta-analysis rank
mergeMagmaMetaRank(magma, corticalLabonte, corticalLabonteDir, corticalDing, corticalDingDir, corticalRamaker, corticalRamakerDir,"rank") %>% write_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/CorticalTableMetaMagmaRank.csv"))
# Subcortical meta-analysis
mergeMagmaMetaRank(magma, subcorticalLabonte, subcorticalLabonteDir, subcorticalDing, subcorticalDingDir, subcorticalRamaker, subcorticalRamakerDir,"meta") %>% write_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/SubcorticalTableMetaMagma.csv"))
# Subcortical meta-analysis rank
mergeMagmaMetaRank(magma, subcorticalLabonte, subcorticalLabonteDir, subcorticalDing, subcorticalDingDir, subcorticalRamaker, subcorticalRamakerDir,"rank") %>% write_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/SubcorticalTableMetaMagmaRank.csv"))
# Sex-interaction Full meta-analysis
mergeMagmaMetaRank(magma, fullLabonte_flipped, fullLabonteDir, fullDing_flipped, fullDingDir, fullRamaker_flipped, fullRamakerDir,"meta") %>% write_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/SIFullTableMetaMagma.csv"))
# Sex-interaction Full meta-analysis rank
mergeMagmaMetaRank(magma, fullLabonte_flipped, fullLabonteDir, fullDing_flipped, fullDingDir, fullRamaker_flipped, fullRamakerDir,"rank") %>% write_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/SIFullTableMetaMagmaRank.csv"))
# Sex-interaction Cortical meta-analysis
mergeMagmaMetaRank(magma, corticalLabonte_flipped, corticalLabonteDir, corticalDing_flipped, corticalDingDir, corticalRamaker_flipped, corticalRamakerDir,"meta") %>% write_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/SICorticalTableMetaMagma.csv"))
# Sex-interaction Cortical meta-analysis rank
mergeMagmaMetaRank(magma, corticalLabonte_flipped, corticalLabonteDir, corticalDing_flipped, corticalDingDir, corticalRamaker_flipped, corticalRamakerDir,"rank") %>% write_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/SICorticalTableMetaMagmaRank.csv"))
# Sex-interaction Subcortical meta-analysis
mergeMagmaMetaRank(magma, subcorticalLabonte_flipped, subcorticalLabonteDir, subcorticalDing_flipped, subcorticalDingDir, subcorticalRamaker_flipped, subcorticalRamakerDir,"meta") %>% write_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/SISubcorticalTableMetaMagma.csv"))
# Sex-interaction Subcortical meta-analysis rank
mergeMagmaMetaRank(magma, subcorticalLabonte_flipped, subcorticalLabonteDir, subcorticalDing_flipped, subcorticalDingDir, subcorticalRamaker_flipped, subcorticalRamakerDir,"rank") %>% write_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/SISubcorticalTableMetaMagmaRank.csv"))

#-----------------------------------------------------------------------------------------
# Get min p values
# Combine study-specific meta-analysis
full_lab <- get_p(Howard, fullLabonte, 'Labonte')
full_ding <- get_p(Howard, fullDing, 'Ding')
full_ramaker <- get_p(Howard, fullRamaker, 'Ramaker')
full_p <- rbind(full_lab, full_ding, full_ramaker)

female_lab <- get_p(Howard, femaleLabonte, 'Labonte')
female_ding <- get_p(Howard, femaleDing, 'Ding')
female_ramaker <- get_p(Howard, femaleRamaker, 'Ramaker')
female_p <- rbind(female_lab, female_ding, female_ramaker)

male_lab <- get_p(Howard, maleLabonte, 'Labonte')
male_ding <- get_p(Howard, maleDing, 'Ding')
male_ramaker <- get_p(Howard, maleRamaker, 'Ramaker')
male_p <- rbind(male_lab, male_ding, male_ramaker)

cortical_lab <- get_p(Howard, corticalLabonte, 'Labonte')
cortical_ding <- get_p(Howard, corticalDing, 'Ding')
cortical_ramaker <- get_p(Howard, corticalRamaker, 'Ramaker')
cortical_p <- rbind(cortical_lab, cortical_ding, cortical_ramaker)

subcortical_lab <- get_p(Howard, subcorticalLabonte, 'Labonte')
subcortical_ding <- get_p(Howard, subcorticalDing, 'Ding')
subcortical_ramaker <- get_p(Howard, subcorticalRamaker, 'Ramaker')
subcortical_p <- rbind(subcortical_lab, subcortical_ding, subcortical_ramaker)

full_flip_lab <- get_p(Howard, fullLabonte_flipped, 'Labonte')
full_flip_ding <- get_p(Howard, fullDing_flipped, 'Ding')
full_flip_ramaker <- get_p(Howard, fullRamaker_flipped, 'Ramaker')
full_p_flip <- rbind(full_flip_lab, full_flip_ding, full_flip_ramaker)

cortical_flip_lab <- get_p(Howard, corticalLabonte_flipped, 'Labonte')
cortical_flip_ding <- get_p(Howard, corticalDing_flipped, 'Ding')
cortical_flip_ramaker <- get_p(Howard, corticalRamaker_flipped, 'Ramaker')
cortical_p_flip <- rbind(cortical_flip_lab, cortical_flip_ding, cortical_flip_ramaker)

subcortical_flip_lab <- get_p(Howard, subcorticalLabonte_flipped, 'Labonte')
subcortical_flip_ding <- get_p(Howard, subcorticalDing_flipped, 'Ding')
subcortical_flip_ramaker <- get_p(Howard, subcorticalRamaker_flipped, 'Ramaker')
subcortical_p_flip <- rbind(subcortical_flip_lab, subcortical_flip_ding, subcortical_flip_ramaker)

#-----------------------------------------------------------------------------------------
#Perform meta-analysis on the 269 GWAS identified depression genes
#full analysis
fullTable <- mergeMetaStudies(Howard, fullLabonte, fullLabonteDir, fullDing, fullDingDir, fullRamaker, fullRamakerDir)
fullTable %<>% MetaAnalysis(full_p)
fullTable %<>% left_join(DE_Prior %>% dplyr::select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
#fullTable %<>% left_join(Howard_Polygenics %>% dplyr::select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
fullTable %<>% left_join(Howard_Polygenics %>% dplyr::select(-mouseGene), by = c('gene_symbol' = 'gene_symbol', 'gene_name'='gene_name'))
fullTable %<>% mutate_all(replace_na, replace = "Not_Available")

#Female analysis
femaleTable <- mergeMetaStudies(Howard,femaleLabonte, femaleLabonteDir, femaleDing, femaleDingDir, femaleRamaker, femaleRamakerDir)
femaleTable %<>% MetaAnalysis(female_p)
femaleTable %<>% left_join(DE_Prior %>% dplyr::select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
#femaleTable %<>% left_join(Howard_Polygenics%>% dplyr::select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
femaleTable %<>% left_join(Howard_Polygenics%>% dplyr::select(-mouseGene), by = c('gene_symbol' = 'gene_symbol', 'gene_name'='gene_name'))
femaleTable %<>% mutate_all(replace_na, replace = "Not_Available")

#male analysis
maleTable <- mergeMetaStudies(Howard, maleLabonte, maleLabonteDir,maleDing, maleDingDir, maleRamaker, maleRamakerDir)
maleTable %<>% MetaAnalysis(male_p)
maleTable %<>% left_join(DE_Prior %>% dplyr::select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
#maleTable %<>% left_join(Howard_Polygenics%>% dplyr::select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
maleTable %<>% left_join(Howard_Polygenics%>% dplyr::select(-mouseGene), by = c('gene_symbol' = 'gene_symbol', 'gene_name'='gene_name'))
maleTable %<>% mutate_all(replace_na, replace = "Not_Available")

# #cortical analysis
corticalTable <- mergeMetaStudies(Howard, corticalLabonte, corticalLabonteDir, corticalDing, corticalDingDir, corticalRamaker, corticalRamakerDir)
corticalTable %<>% MetaAnalysis(cortical_p)
corticalTable %<>% left_join(DE_Prior %>% dplyr::select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name')) 
#corticalTable %<>% left_join(Howard_Polygenics%>% dplyr::select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
corticalTable %<>% left_join(Howard_Polygenics%>% dplyr::select(-mouseGene), by = c('gene_symbol' = 'gene_symbol', 'gene_name'='gene_name'))
corticalTable %<>% mutate_all(replace_na, replace = "Not_Available")

# subcortical analysis
subcorticalTable <- mergeMetaStudies(Howard, subcorticalLabonte, subcorticalLabonteDir, subcorticalDing, subcorticalDingDir, subcorticalRamaker, subcorticalRamakerDir)
subcorticalTable %<>% MetaAnalysis(subcortical_p)
subcorticalTable %<>% left_join(DE_Prior %>% dplyr::select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name')) 
#csuborticalTable %<>% left_join(Howard_Polygenics%>% dplyr::select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
subcorticalTable %<>% left_join(Howard_Polygenics%>% dplyr::select(-mouseGene), by = c('gene_symbol' = 'gene_symbol', 'gene_name'='gene_name'))
subcorticalTable %<>% mutate_all(replace_na, replace = "Not_Available")

#Sex-interaction analyses
fullTable_Flip <- mergeMetaStudies(Howard, fullLabonte_flipped, fullLabonteDir, fullDing_flipped, fullDingDir, fullRamaker_flipped, fullRamakerDir)
fullTable_Flip %<>% MetaAnalysis(full_p_flip)
fullTable_Flip %<>% left_join(DE_Prior %>% dplyr::select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
#fullTable_Flip %<>% left_join(Howard_Polygenics%>% dplyr::select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
fullTable_Flip %<>% left_join(Howard_Polygenics%>% dplyr::select( -mouseGene), by = c('gene_symbol' = 'gene_symbol', 'gene_name'='gene_name'))
fullTable_Flip %<>% mutate_all(replace_na, replace = "Not_Available")

corticalTable_Flip <- mergeMetaStudies(Howard, corticalLabonte_flipped, corticalLabonteDir, corticalDing_flipped, corticalDingDir, corticalRamaker_flipped, corticalRamakerDir)
corticalTable_Flip %<>% MetaAnalysis(cortical_p_flip)
corticalTable_Flip %<>% left_join(DE_Prior %>% dplyr::select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name')) 
#corticalTable_Flip %<>% left_join(Howard_Polygenics%>% dplyr::select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
corticalTable_Flip %<>% left_join(Howard_Polygenics%>% dplyr::select(-mouseGene), by = c('gene_symbol' = 'gene_symbol', 'gene_name'='gene_name'))
corticalTable_Flip %<>% mutate_all(replace_na, replace = "Not_Available")

subcorticalTable_Flip <- mergeMetaStudies(Howard, subcorticalLabonte_flipped, subcorticalLabonteDir, subcorticalDing_flipped, subcorticalDingDir, subcorticalRamaker_flipped, subcorticalRamakerDir)
subcorticalTable_Flip %<>% MetaAnalysis(subcortical_p_flip)
subcorticalTable_Flip %<>% left_join(DE_Prior %>% dplyr::select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name')) 
#subcorticalTable_Flip %<>% left_join(Howard_Polygenics%>% dplyr::select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
subcorticalTable_Flip %<>% left_join(Howard_Polygenics%>% dplyr::select(-mouseGene), by = c('gene_symbol' = 'gene_symbol', 'gene_name'='gene_name'))
subcorticalTable_Flip %<>% mutate_all(replace_na, replace = "Not_Available")


########################################################
### Perform the Ranking analysis on all above tables ##
#######################################################

fullTableRank <- mergeMetaStudies(Howard, fullLabonte, fullLabonteDir, fullDing, fullDingDir, fullRamaker, fullRamakerDir)
fullTableRank %<>% GenomeRank()
fullTableRank %<>% left_join(DE_Prior %>% dplyr::select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name')) 
#fullTableRank %<>% left_join(Howard_Polygenics %>% dplyr::select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
fullTableRank %<>% left_join(Howard_Polygenics %>% dplyr::select(-mouseGene), by = c('gene_symbol' = 'gene_symbol', 'gene_name'='gene_name'))
fullTableRank %<>% mutate_all(replace_na, replace = "Not_Available")

femaleTableRank <- mergeMetaStudies(Howard,femaleLabonte, femaleLabonteDir, femaleDing, femaleDingDir, femaleRamaker, femaleRamakerDir)
femaleTableRank %<>% GenomeRank()
femaleTableRank %<>% left_join(DE_Prior %>% dplyr::select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
femaleTableRank %<>% left_join(Howard_Polygenics%>% dplyr::select( -mouseGene), by = c('gene_symbol' = 'gene_symbol', 'gene_name'='gene_name'))
femaleTableRank %<>% mutate_all(replace_na, replace = "Not_Available")

maleTableRank <- mergeMetaStudies(Howard, maleLabonte, maleLabonteDir,maleDing, maleDingDir, maleRamaker, maleRamakerDir)
maleTableRank %<>% GenomeRank()
maleTableRank %<>% left_join(DE_Prior %>% dplyr::select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
maleTableRank %<>% left_join(Howard_Polygenics%>% dplyr::select( -mouseGene), by = c('gene_symbol' = 'gene_symbol', 'gene_name'='gene_name'))
maleTableRank %<>% mutate_all(replace_na, replace = "Not_Available")

corticalTableRank <- mergeMetaStudies(Howard, corticalLabonte, corticalLabonteDir, corticalDing, corticalDingDir, corticalRamaker, corticalRamakerDir)
corticalTableRank %<>% GenomeRank()
corticalTableRank %<>% left_join(DE_Prior %>% dplyr::select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
corticalTableRank %<>% left_join(Howard_Polygenics%>% dplyr::select( -mouseGene), by = c('gene_symbol' = 'gene_symbol', 'gene_name'='gene_name'))
corticalTableRank %<>% mutate_all(replace_na, replace = "Not_Available")

subcorticalTableRank <- mergeMetaStudies(Howard, subcorticalLabonte, subcorticalLabonteDir, subcorticalDing, subcorticalDingDir, subcorticalRamaker, subcorticalRamakerDir)
subcorticalTableRank %<>% GenomeRank()
subcorticalTableRank %<>% left_join(DE_Prior %>% dplyr::select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
subcorticalTableRank %<>% left_join(Howard_Polygenics%>% dplyr::select( -mouseGene), by = c('gene_symbol' = 'gene_symbol', 'gene_name'='gene_name'))
subcorticalTableRank %<>% mutate_all(replace_na, replace = "Not_Available")


fullTableRank_Flip <- mergeMetaStudies(Howard, fullLabonte_flipped, fullLabonteDir, fullDing_flipped, fullDingDir, fullRamaker_flipped, fullRamakerDir)
fullTableRank_Flip %<>% GenomeRank()
fullTableRank_Flip %<>% left_join(DE_Prior %>% dplyr::select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
fullTableRank_Flip %<>% left_join(Howard_Polygenics%>% dplyr::select( -mouseGene), by = c('gene_symbol' = 'gene_symbol', 'gene_name'='gene_name'))
fullTableRank_Flip %<>% mutate_all(replace_na, replace = "Not_Available")

corticalTableRank_Flip <- mergeMetaStudies(Howard, corticalLabonte_flipped, corticalLabonteDir, corticalDing_flipped, corticalDingDir, corticalRamaker_flipped, corticalRamakerDir)
corticalTableRank_Flip %<>% GenomeRank()
corticalTableRank_Flip %<>% left_join(DE_Prior %>% dplyr::select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
corticalTableRank_Flip %<>% left_join(Howard_Polygenics%>% dplyr::select( -mouseGene), by = c('gene_symbol' = 'gene_symbol', 'gene_name'='gene_name'))
corticalTableRank_Flip %<>% mutate_all(replace_na, replace = "Not_Available")

subcorticalTableRank_Flip <- mergeMetaStudies(Howard, subcorticalLabonte_flipped, subcorticalLabonteDir, subcorticalDing_flipped, subcorticalDingDir, subcorticalRamaker_flipped, subcorticalRamakerDir)
subcorticalTableRank_Flip %<>% GenomeRank()
subcorticalTableRank_Flip %<>% left_join(DE_Prior %>% dplyr::select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
subcorticalTableRank_Flip %<>% left_join(Howard_Polygenics%>% dplyr::select( -mouseGene), by = c('gene_symbol' = 'gene_symbol', 'gene_name'='gene_name'))
subcorticalTableRank_Flip %<>% mutate_all(replace_na, replace = "Not_Available")

#Save files
dir.create(here('Results', 'Tables', 'Meta_Analysis'), recursive = TRUE)

write_csv(fullTable, here('Results', 'Tables', 'Meta_Analysis', "Full_Meta_Analysis.csv"))
write_csv(femaleTable, here('Results', 'Tables', 'Meta_Analysis', 'Female_Meta_Analysis.csv'))
write_csv(maleTable, here('Results', 'Tables', 'Meta_Analysis', 'Male_Meta_Analysis.csv'))
write_csv(corticalTable, here('Results', 'Tables', 'Meta_Analysis', 'Cortical_Meta_Analysis.csv'))
write_csv(subcorticalTable, here('Results', 'Tables', 'Meta_Analysis', 'Subcortical_Meta_Analysis.csv'))
write_csv(fullTable_Flip, here('Results', 'Tables', 'Meta_Analysis','Sex_interaction_Full_Meta_Analysis.csv'))
write_csv(corticalTable_Flip, here('Results', 'Tables', 'Meta_Analysis', 'Sex_interaction_Cortical_Meta_Analysis.csv'))
write_csv(subcorticalTable_Flip, here('Results', 'Tables', 'Meta_Analysis', 'Sex_interaction_Subcortical_Meta_Analysis.csv'))

dir.create(here('Results', 'Tables', 'Genome_Percentile_Rank'), recursive = TRUE)
write_csv(fullTableRank, here('Results', 'Tables', 'Genome_Percentile_Rank', "Full_Meta_Analysis_Rank.csv"))
write_csv(femaleTableRank, here('Results', 'Tables', 'Genome_Percentile_Rank', 'Female_Meta_Analysis_Rank.csv'))
write_csv(maleTableRank, here('Results', 'Tables', 'Genome_Percentile_Rank', 'Male_Meta_Analysis_Rank.csv'))
write_csv(corticalTableRank, here('Results', 'Tables', 'Genome_Percentile_Rank', 'Cortical_Meta_Analysis_Rank.csv'))
write_csv(subcorticalTableRank, here('Results', 'Tables', 'Genome_Percentile_Rank', 'Subcortical_Meta_Analysis_Rank.csv'))
write_csv(fullTableRank_Flip, here('Results', 'Tables', 'Genome_Percentile_Rank', 'Sex_interaction_Full_Meta_Analysis_Rank.csv'))
write_csv(corticalTableRank_Flip, here('Results', 'Tables', 'Genome_Percentile_Rank', 'Sex_interaction_Cortical_Meta_Analysis_Rank.csv'))
write_csv(subcorticalTableRank_Flip, here('Results', 'Tables', 'Genome_Percentile_Rank', 'Sex_interaction_Subcortical_Meta_Analysis_Rank.csv'))



