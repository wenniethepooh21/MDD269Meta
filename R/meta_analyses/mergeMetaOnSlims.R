library(readr)
library(magrittr)
library(dplyr)
library(tidyr)
library(here)
library(googledrive)
library(googlesheets4)

use_gdrive <- FALSE

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
full_p <- rbind(fullLabonte %>% dplyr::select(gene_symbol, min_p_across_regions) %>% mutate(min_p_study = 'Labonte'), fullDing %>% dplyr::select(gene_symbol, min_p_across_regions) %>% mutate(min_p_study = 'Ding'),
                     fullRamaker %>% dplyr::select(gene_symbol, min_p_across_regions) %>% mutate(min_p_study = 'Ramaker')) 
female_p <- rbind(femaleLabonte %>% dplyr::select(gene_symbol, min_p_across_regions) %>% mutate(min_p_study = 'Labonte'), femaleDing %>% dplyr::select(gene_symbol, min_p_across_regions) %>% mutate(min_p_study = 'Ding'),
                femaleRamaker %>% dplyr::select(gene_symbol, min_p_across_regions) %>% mutate(min_p_study = 'Ramaker')) 
male_p <- rbind(maleLabonte %>% dplyr::select(gene_symbol, min_p_across_regions) %>% mutate(min_p_study = 'Labonte'), maleDing %>% dplyr::select(gene_symbol, min_p_across_regions) %>% mutate(min_p_study = 'Ding'),
                maleRamaker %>% dplyr::select(gene_symbol, min_p_across_regions) %>% mutate(min_p_study = 'Ramaker')) 
cortical_p <- rbind(corticalLabonte %>% dplyr::select(gene_symbol, min_p_across_regions) %>% mutate(min_p_study = 'Labonte'), corticalDing %>% dplyr::select(gene_symbol, min_p_across_regions) %>% mutate(min_p_study = 'Ding'),
                corticalRamaker %>% dplyr::select(gene_symbol, min_p_across_regions) %>% mutate(min_p_study = 'Ramaker')) 
subcortical_p <- rbind(subcorticalLabonte %>% dplyr::select(gene_symbol, min_p_across_regions) %>% mutate(min_p_study = 'Labonte'), subcorticalDing %>% dplyr::select(gene_symbol, min_p_across_regions) %>% mutate(min_p_study = 'Ding'),
                subcorticalRamaker %>% dplyr::select(gene_symbol, min_p_across_regions) %>% mutate(min_p_study = 'Ramaker')) 
full_p_flip <- rbind(fullLabonte_flipped %>% dplyr::select(gene_symbol, min_p_across_regions) %>% mutate(min_p_study = 'Labonte'), fullDing_flipped %>% dplyr::select(gene_symbol, min_p_across_regions) %>% mutate(min_p_study = 'Ding'),
                fullRamaker_flipped %>% dplyr::select(gene_symbol, min_p_across_regions) %>% mutate(min_p_study = 'Ramaker')) 
cortical_p_flip <- rbind(corticalLabonte_flipped %>% dplyr::select(gene_symbol, min_p_across_regions) %>% mutate(min_p_study = 'Labonte'), corticalDing_flipped %>% dplyr::select(gene_symbol, min_p_across_regions) %>% mutate(min_p_study = 'Ding'),
                    corticalRamaker_flipped %>% dplyr::select(gene_symbol, min_p_across_regions) %>% mutate(min_p_study = 'Ramaker')) 
subcortical_p_flip <- rbind(subcorticalLabonte_flipped %>% dplyr::select(gene_symbol, min_p_across_regions) %>% mutate(min_p_study = 'Labonte'), subcorticalDing_flipped %>% dplyr::select(gene_symbol, min_p_across_regions) %>% mutate(min_p_study = 'Ding'),
                       subcorticalRamaker_flipped %>% dplyr::select(gene_symbol, min_p_across_regions) %>% mutate(min_p_study = 'Ramaker')) 


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
fullTableRank %<>% left_join(Howard_Polygenics %>% dplyr::select(-mouseGene), by = c('gene_symbol' = 'gene_symbol'))
fullTableRank %<>% mutate_all(replace_na, replace = "Not_Available")

femaleTableRank <- mergeMetaStudies(Howard,femaleLabonte, femaleLabonteDir, femaleDing, femaleDingDir, femaleRamaker, femaleRamakerDir)
femaleTableRank %<>% GenomeRank()
femaleTableRank %<>% left_join(DE_Prior %>% dplyr::select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
femaleTableRank %<>% left_join(Howard_Polygenics%>% dplyr::select( -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
femaleTableRank %<>% mutate_all(replace_na, replace = "Not_Available")

maleTableRank <- mergeMetaStudies(Howard, maleLabonte, maleLabonteDir,maleDing, maleDingDir, maleRamaker, maleRamakerDir)
maleTableRank %<>% GenomeRank()
maleTableRank %<>% left_join(DE_Prior %>% dplyr::select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
maleTableRank %<>% left_join(Howard_Polygenics%>% dplyr::select( -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
maleTableRank %<>% mutate_all(replace_na, replace = "Not_Available")

corticalTableRank <- mergeMetaStudies(Howard, corticalLabonte, corticalLabonteDir, corticalDing, corticalDingDir, corticalRamaker, corticalRamakerDir)
corticalTableRank %<>% GenomeRank()
corticalTableRank %<>% left_join(DE_Prior %>% dplyr::select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
corticalTableRank %<>% left_join(Howard_Polygenics%>% dplyr::select( -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
corticalTableRank %<>% mutate_all(replace_na, replace = "Not_Available")

subcorticalTableRank <- mergeMetaStudies(Howard, subcorticalLabonte, subcorticalLabonteDir, subcorticalDing, subcorticalDingDir, subcorticalRamaker, subcorticalRamakerDir)
subcorticalTableRank %<>% GenomeRank()
subcorticalTableRank %<>% left_join(DE_Prior %>% dplyr::select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
subcorticalTableRank %<>% left_join(Howard_Polygenics%>% dplyr::select( -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
subcorticalTableRank %<>% mutate_all(replace_na, replace = "Not_Available")


fullTableRank_Flip <- mergeMetaStudies(Howard, fullLabonte_flipped, fullLabonteDir, fullDing_flipped, fullDingDir, fullRamaker_flipped, fullRamakerDir)
fullTableRank_Flip %<>% GenomeRank()
fullTableRank_Flip %<>% left_join(DE_Prior %>% dplyr::select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
fullTableRank_Flip %<>% left_join(Howard_Polygenics%>% dplyr::select( -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
fullTableRank_Flip %<>% mutate_all(replace_na, replace = "Not_Available")

corticalTableRank_Flip <- mergeMetaStudies(Howard, corticalLabonte_flipped, corticalLabonteDir, corticalDing_flipped, corticalDingDir, corticalRamaker_flipped, corticalRamakerDir)
corticalTableRank_Flip %<>% GenomeRank()
corticalTableRank_Flip %<>% left_join(DE_Prior %>% dplyr::select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
corticalTableRank_Flip %<>% left_join(Howard_Polygenics%>% dplyr::select( -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
corticalTableRank_Flip %<>% mutate_all(replace_na, replace = "Not_Available")

subcorticalTableRank_Flip <- mergeMetaStudies(Howard, subcorticalLabonte_flipped, subcorticalLabonteDir, subcorticalDing_flipped, subcorticalDingDir, subcorticalRamaker_flipped, subcorticalRamakerDir)
subcorticalTableRank_Flip %<>% GenomeRank()
subcorticalTableRank_Flip %<>% left_join(DE_Prior %>% dplyr::select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
subcorticalTableRank_Flip %<>% left_join(Howard_Polygenics%>% dplyr::select( -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
subcorticalTableRank_Flip %<>% mutate_all(replace_na, replace = "Not_Available")


#FULL DATA SHEETS
#Access googlesheets to upload the tables online for an interactive experience

if (use_gdrive == TRUE) {
    sheets_auth(token = drive_token())
    
    ma <- drive_get("~/Thesis/Manuscript/Tables/Full_Tables/Meta_Analysis")
    if(nrow(ma) != 0) {
      drive_rm(ma)
    }
    #create the google worksheet
    ma <- sheets_create("Meta_Analysis", sheets = c('Full_Meta_Analysis', 'Male_Meta_Analysis','Female_Meta_Analysis','Cortical_Meta_Analysis','Subcortical_Meta_Analysis', 'Sex_interaction_Full_Meta_Analysis', 'Sex_interaction_Cortical_Meta_Analysis', 'Sex_interaction_Subcortical_Meta_Analysis'))
    sheets_write(fullTable, ma, sheet = "Full_Meta_Analysis")
    sheets_write(femaleTable, ma, sheet = 'Female_Meta_Analysis')
    sheets_write(maleTable, ma,'Male_Meta_Analysis')
    sheets_write(corticalTable, ma, 'Cortical_Meta_Analysis')
    sheets_write(subcorticalTable, ma, 'Subcortical_Meta_Analysis')
    sheets_write(fullTable_Flip, ma, 'Sex_interaction_Full_Meta_Analysis')
    sheets_write(corticalTable_Flip, ma, 'Sex_interaction_Cortical_Meta_Analysis')
    sheets_write(subcorticalTable_Flip, ma, 'Sex_interaction_Subcortical_Meta_Analysis')
    
    drive_mv(file = "Meta_Analysis", path = "~/Thesis/Manuscript/Tables/Full_Tables/")  # move Sheets file
    
    mar <- drive_get("Thesis/Manuscript/Tables/Full_Tables/Genome_Percentile_Rank")
    if(nrow(mar) != 0) {
      drive_rm(mar)
    }
    #create the google worksheet
    mar <- sheets_create("Genome_Percentile_Rank", sheets = c('Full_Meta_Analysis_Rank', 'Male_Meta_Analysis_Rank','Female_Meta_Analysis_Rank', 'Cortical_Meta_Analysis_Rank','Subcortical_Meta_Analysis_Rank', 'Sex_interaction_Full_Meta_Analysis_Rank',  'Sex_interaction_Cortical_Meta_Analysis_Rank',  'Sex_interaction_Subcortical_Meta_Analysis_Rank'))
    sheets_write(fullTableRank, mar, sheet = "Full_Meta_Analysis_Rank")
    sheets_write(femaleTableRank, mar, sheet = 'Female_Meta_Analysis_Rank')
    sheets_write(maleTableRank, mar,'Male_Meta_Analysis_Rank')
    sheets_write(corticalTableRank, mar, 'Cortical_Meta_Analysis_Rank')
    sheets_write(subcorticalTableRank, mar, 'Subcortical_Meta_Analysis_Rank')
    sheets_write(fullTableRank_Flip, mar, 'Sex_interaction_Full_Meta_Analysis_Rank')
    sheets_write(corticalTableRank_Flip, mar, 'Sex_interaction_Cortical_Meta_Analysis_Rank')
    sheets_write(subcorticalTableRank_Flip, mar, 'Sex_interaction_Subcortical_Meta_Analysis_Rank')
    
    drive_mv(file = "Genome_Percentile_Rank", path = "~/Thesis/Manuscript/Tables/Full_Tables/")  # move Sheets file
    
    ma <- drive_get("~/Thesis/Manuscript/Tables/Slim_Tables/Official_Meta_Analysis")
    if(nrow(ma) != 0) {
      drive_rm(ma)
    }
    #create the google worksheet
    ma <- sheets_create("Official_Meta_Analysis", sheets = c('First_model_Full_Meta_Analysis', 'First_model_Male_Meta_Analysis','First_model_Female_Meta_Analysis','First_model_Cortical_Meta_Analysis','First_model_Subcortical_Meta_Analysis','Sex_interaction_Full_Meta_Analysis', 'Sex_interaction_Cortical_Meta_Analysis', 'Sex_interaction_Subcortical_Meta_Analysis'))
    sheets_write(fullTable %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_p, Bonferroni_meta_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), ma, sheet = "First_model_Full_Meta_Analysis")
    sheets_write(femaleTable %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_p, Bonferroni_meta_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), ma, sheet = 'First_model_Female_Meta_Analysis')
    sheets_write(maleTable %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_p, Bonferroni_meta_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), ma,'First_model_Male_Meta_Analysis')
    sheets_write(corticalTable %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_p, Bonferroni_meta_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), ma, 'First_model_Cortical_Meta_Analysis')
    sheets_write(subcorticalTable %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_p, Bonferroni_meta_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), ma, 'First_model_Subcortical_Meta_Analysis')
    sheets_write(fullTable_Flip  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_p, Bonferroni_meta_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), ma, sheet = "Sex_interaction_Full_Meta_Analysis")
    sheets_write(corticalTable_Flip  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_p, Bonferroni_meta_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), ma, 'Sex_interaction_Cortical_Meta_Analysis')
    sheets_write(subcorticalTable_Flip  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_p, Bonferroni_meta_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), ma, 'Sex_interaction_Subcortical_Meta_Analysis')
    
    drive_mv(file = "Official_Meta_Analysis", path = "~/Thesis/Manuscript/Tables/Slim_Tables/")  # move Sheets file
    
    mar <- drive_get("Thesis/Manuscript/Tables/Slim_Tables/Official_Genome_Percentile_Rank")
    if(nrow(mar) != 0) {
      drive_rm(mar)
    }
    #create the google worksheet
      mar <- sheets_create("Official_Genome_Percentile_Rank", sheets = c('First_model_Full_Meta_Analysis_Rank', 'First_model_Male_Meta_Analysis_Rank','First_model_Female_Meta_Analysis_Rank', 'First_model_Cortical_Meta_Analysis_Rank', 'First_model_Subcortical_Meta_Analysis_Rank','Sex_interaction_Full_Meta_Analysis_Rank', 'Sex_interaction_Cortical_Meta_Analysis_Rank', 'Sex_interaction_Subcortical_Meta_Analysis_Rank'))
      sheets_write(fullTableRank  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), mar, sheet = "First_model_Full_Meta_Analysis_Rank")
      sheets_write(femaleTableRank  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), mar, sheet = 'First_model_Female_Meta_Analysis_Rank')
      sheets_write(maleTableRank  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), mar,'First_model_Male_Meta_Analysis_Rank')
      sheets_write(corticalTableRank  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), mar, 'First_model_Cortical_Meta_Analysis_Rank')
      sheets_write(subcorticalTableRank  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), mar, 'First_model_Subcortical_Meta_Analysis_Rank')
      sheets_write(fullTableRank_Flip  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), mar, sheet = "Sex_interaction_Full_Meta_Analysis_Rank")
      sheets_write(corticalTableRank_Flip  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), mar, 'Sex_interaction_Cortical_Meta_Analysis_Rank')
      sheets_write(subcorticalTableRank_Flip  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), mar, 'Sex_interaction_Subcortical_Meta_Analysis_Rank')
      drive_mv(file = "Official_Genome_Percentile_Rank", path = "~/Thesis/Manuscript/Tables/Slim_Tables/")  # move Sheets file
      
    } else {
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
      
  
}

