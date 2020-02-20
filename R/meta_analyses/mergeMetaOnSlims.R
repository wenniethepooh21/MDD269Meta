library(readr)
library(magrittr)
library(dplyr)
library(tidyr)
library(here)
library(googledrive)
library(googlesheets4)

fullLabonte <- read_csv(here("Processed_Data/LabonteEtAl/FullLabonteTableMagma.csv"))
femaleLabonte <- read_csv(here("Processed_Data/LabonteEtAl/FemaleLabonteTableMagma.csv")) 
maleLabonte <- read_csv(here("Processed_Data/LabonteEtAl/MaleLabonteTableMagma.csv"))
corticalLabonte <- read_csv( here("Processed_Data/LabonteEtAl/CorticalLabonteTableMagma.csv"))
fullLabonte_flipped <- read_csv(here("Processed_Data/LabonteEtAl/FullLabonteTableMagma_flipped.csv"))
corticalLabonte_flipped <- read_csv( here("Processed_Data/LabonteEtAl/CorticalLabonteTableMagma_flipped.csv"))

fullDing <- read_csv(here("Processed_Data/DingEtAl/FullDingTableMagma.csv"))
femaleDing <-read_csv(here("Processed_Data/DingEtAl/FemaleDingTableMagma.csv"))
maleDing <-read_csv(here("Processed_Data/DingEtAl/MaleDingTableMagma.csv"))
corticalDing <- read_csv( here("Processed_Data/DingEtAl/CorticalDingTableMagma.csv"))
fullDing_flipped <- read_csv(here("Processed_Data/DingEtAl/FullDingTableMagma_flipped.csv"))
corticalDing_flipped <- read_csv(here("Processed_Data/DingEtAl/CorticalDingTableMagma_flipped.csv"))

fullRamaker <- read_csv(here("Processed_Data/RamakerEtAl/FullRamakerTableMagma.csv"))
femaleRamaker <- read_csv(here("Processed_Data/RamakerEtAl/FemaleRamakerTableMagma.csv"))
maleRamaker <- read_csv(here("Processed_Data/RamakerEtAl/MaleRamakerTableMagma.csv"))
corticalRamaker <- read_csv(here("Processed_Data/RamakerEtAl/CorticalRamakerTableMagma.csv"))
fullRamaker_flipped <- read_csv(here("Processed_Data/RamakerEtAl/FullRamakerTableMagma_flipped.csv"))
corticalRamaker_flipped <- read_csv(here("Processed_Data/RamakerEtAl/CorticalRamakerTableMagma_flipped.csv"))

#read in the equivalent Howard genes used in each transcriptomic study
Howard <- read_csv(here("Processed_Data/HowardEtAl/fullHowardTable.csv"))
#read in the differentially expressed "prior ranking" for the genes identified in Crow, et al. 
DE_Prior <- read_tsv(here("Raw_Data/CrowEtAl/pnas.1802973116.sd02.txt"))
#Read in the associated transcriptomic cell types and brain regions that maximally expressed each Howard gene
Howard_Polygenics <- read_csv(here("Processed_Data/HowardEtAl/HowardRegionsPolygenicCellTypes_four.csv"))

# Column names - for the direction of expression in each brain region
fullLabonteDir <- "female_Anterior_Insula_female_BA11_female_BA25_female_BA8/9_female_Nac_female_Subic_male_Anterior_Insula_male_BA11_male_BA25_male_BA8/9_male_Nac_male_Subic"
femaleLabonteDir <- "female_Anterior_Insula_female_BA11_female_BA25_female_BA8/9_female_Nac_female_Subic"
maleLabonteDir <- "male_Anterior_Insula_male_BA11_male_BA25_male_BA8/9_male_Nac_male_Subic"
corticalLabonteDir <- "female_Anterior_Insula_female_BA11_female_BA25_female_BA8/9_male_Anterior_Insula_male_BA11_male_BA25_male_BA8/9"

fullDingDir <- "female_ACC.1_female_ACC.2_female_AMY_female_DLPFC_male_ACC.1_male_ACC.2_male_AMY_male_DLPFC"
femaleDingDir <-"female_ACC.1_female_ACC.2_female_AMY_female_DLPFC"
maleDingDir <-"male_ACC.1_male_ACC.2_male_AMY_male_DLPFC"
corticalDingDir <-"female_ACC.1_female_ACC.2_female_DLPFC_male_ACC.1_male_ACC.2_male_DLPFC"

fullRamakerDir <- "AnCg.F_nAcc.F_DLPFC.F_AnCg.M_nAcc.M_DLPFC.M"
femaleRamakerDir <- "AnCg_nAcc_DLPFC_Female_directions"
maleRamakerDir <- "AnCg_nAcc_DLPFC_Male_directions"
corticalRamakerDir <- "AnCg.F_DLPFC.F_AnCg.M_DLPFC.M"

#Perform the meta-analysis combining all study-specific meta analyses together
#Functions used for combining the data from each study and meta-analysis calculations can be found here
source(here("R/meta_analyses/mergeAnalysis.R"))

##### Perform meta-analysis on the entire 17,842 MAGMA genes to run the Wilcox Test ####--------------------------------------------------
magma <- read_csv(here("Processed_Data/HowardEtAl/FullMagmaGenes.csv"))
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
# Sex-interaction Full meta-analysis
mergeMagmaMetaRank(magma, fullLabonte_flipped, fullLabonteDir, fullDing_flipped, fullDingDir, fullRamaker_flipped, fullRamakerDir,"meta") %>% write_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/SIFullTableMetaMagma.csv"))
# Sex-interaction Full meta-analysis rank
mergeMagmaMetaRank(magma, fullLabonte_flipped, fullLabonteDir, fullDing_flipped, fullDingDir, fullRamaker_flipped, fullRamakerDir,"rank") %>% write_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/SIFullTableMetaMagmaRank.csv"))
# Sex-interaction Cortical meta-analysis
mergeMagmaMetaRank(magma, corticalLabonte_flipped, corticalLabonteDir, corticalDing_flipped, corticalDingDir, corticalRamaker_flipped, corticalRamakerDir,"meta") %>% write_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/SICorticalTableMetaMagma.csv"))
# Sex-interaction Cortical meta-analysis rank
mergeMagmaMetaRank(magma, corticalLabonte_flipped, corticalLabonteDir, corticalDing_flipped, corticalDingDir, corticalRamaker_flipped, corticalRamakerDir,"rank") %>% write_csv(here("Processed_Data/Meta_Analysis_Results/MAGMA/SICorticalTableMetaMagmaRank.csv"))
#-----------------------------------------------------------------------------------------
#Perform meta-analysis on the 269 GWAS identified depression genes
#full analysis
fullTable <- mergeMetaStudies(Howard, fullLabonte, fullLabonteDir, fullDing, fullDingDir, fullRamaker, fullRamakerDir)
fullTable %<>% MetaAnalysis()
fullTable %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
fullTable %<>% left_join(Howard_Polygenics %>% select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
fullTable %<>% mutate_all(replace_na, replace = "Not_Available")

#Female analysis
femaleTable <- mergeMetaStudies(Howard,femaleLabonte, femaleLabonteDir, femaleDing, femaleDingDir, femaleRamaker, femaleRamakerDir)
femaleTable %<>% MetaAnalysis()
femaleTable %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
femaleTable %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
femaleTable %<>% mutate_all(replace_na, replace = "Not_Available")

#male analysis
maleTable <- mergeMetaStudies(Howard, maleLabonte, maleLabonteDir,maleDing, maleDingDir, maleRamaker, maleRamakerDir)
maleTable %<>% MetaAnalysis()
maleTable %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
maleTable %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
maleTable %<>% mutate_all(replace_na, replace = "Not_Available")

# #cortical analysis
corticalTable <- mergeMetaStudies(Howard, corticalLabonte, corticalLabonteDir, corticalDing, corticalDingDir, corticalRamaker, corticalRamakerDir)
corticalTable %<>% MetaAnalysis()
corticalTable %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name')) 
corticalTable %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
corticalTable %<>% mutate_all(replace_na, replace = "Not_Available")

#Sex-interaction analyses
fullTable_Flip <- mergeMetaStudies(Howard, fullLabonte_flipped, fullLabonteDir, fullDing_flipped, fullDingDir, fullRamaker_flipped, fullRamakerDir)
fullTable_Flip %<>% MetaAnalysis()
fullTable_Flip %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
fullTable_Flip %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
fullTable_Flip %<>% mutate_all(replace_na, replace = "Not_Available")

corticalTable_Flip <- mergeMetaStudies(Howard, corticalLabonte_flipped, corticalLabonteDir, corticalDing_flipped, corticalDingDir, corticalRamaker_flipped, corticalRamakerDir)
corticalTable_Flip %<>% MetaAnalysis()
corticalTable_Flip %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name')) 
corticalTable_Flip %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
corticalTable_Flip %<>% mutate_all(replace_na, replace = "Not_Available")

########################################################
### Perform the Ranking analysis on all above tables ##
#######################################################

fullTableRank <- mergeMetaStudies(Howard, fullLabonte, fullLabonteDir, fullDing, fullDingDir, fullRamaker, fullRamakerDir)
fullTableRank %<>% GenomeRank()
fullTableRank %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name')) 
fullTableRank %<>% left_join(Howard_Polygenics %>% select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
fullTableRank %<>% mutate_all(replace_na, replace = "Not_Available")

femaleTableRank <- mergeMetaStudies(Howard,femaleLabonte, femaleLabonteDir, femaleDing, femaleDingDir, femaleRamaker, femaleRamakerDir)
femaleTableRank %<>% GenomeRank()
femaleTableRank %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
femaleTableRank %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
femaleTableRank %<>% mutate_all(replace_na, replace = "Not_Available")

maleTableRank <- mergeMetaStudies(Howard, maleLabonte, maleLabonteDir,maleDing, maleDingDir, maleRamaker, maleRamakerDir)
maleTableRank %<>% GenomeRank()
maleTableRank %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
maleTableRank %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
maleTableRank %<>% mutate_all(replace_na, replace = "Not_Available")

corticalTableRank <- mergeMetaStudies(Howard, corticalLabonte, corticalLabonteDir, corticalDing, corticalDingDir, corticalRamaker, corticalRamakerDir)
corticalTableRank %<>% GenomeRank()
corticalTableRank %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
corticalTableRank %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
corticalTableRank %<>% mutate_all(replace_na, replace = "Not_Available")

fullTableRank_Flip <- mergeMetaStudies(Howard, fullLabonte_flipped, fullLabonteDir, fullDing_flipped, fullDingDir, fullRamaker_flipped, fullRamakerDir)
fullTableRank_Flip %<>% GenomeRank()
fullTableRank_Flip %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
fullTableRank_Flip %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
fullTableRank_Flip %<>% mutate_all(replace_na, replace = "Not_Available")

corticalTableRank_Flip <- mergeMetaStudies(Howard, corticalLabonte_flipped, corticalLabonteDir, corticalDing_flipped, corticalDingDir, corticalRamaker_flipped, corticalRamakerDir)
corticalTableRank_Flip %<>% GenomeRank()
corticalTableRank_Flip %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
corticalTableRank_Flip %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'gene_symbol'))
corticalTableRank_Flip %<>% mutate_all(replace_na, replace = "Not_Available")

# #FULL DATA SHEETS
# #Access googlesheets to upload the tables online for an interactive experience
# 
# sheets_auth(token = drive_token())
# 
# ma <- drive_get("~/Thesis/Manuscript/Tables/Full_Tables/Meta_Analysis")
# if(nrow(ma) != 0) {
#   drive_rm(ma)
# }
# #create the google worksheet
# ma <- sheets_create("Meta_Analysis", sheets = c('Full_Meta_Analysis', 'Male_Meta_Analysis','Female_Meta_Analysis','Cortical_Meta_Analysis'))
# sheets_write(fullTable, ma, sheet = "Full_Meta_Analysis")
# sheets_write(femaleTable, ma, sheet = 'Female_Meta_Analysis')
# sheets_write(maleTable, ma,'Male_Meta_Analysis')
# sheets_write(corticalTable, ma, 'Cortical_Meta_Analysis')
# 
# drive_mv(file = "Meta_Analysis", path = "~/Thesis/Manuscript/Tables/Full_Tables/")  # move Sheets file
# 
# mar <- drive_get("Thesis/Manuscript/Tables/Full_Tables/Meta_Analysis_Rank")
# if(nrow(mar) != 0) {
#   drive_rm(mar)
# }
# #create the google worksheet
# mar <- sheets_create("Meta_Analysis_Rank", sheets = c('Full_Meta_Analysis_Rank', 'Male_Meta_Analysis_Rank','Female_Meta_Analysis_Rank', 'Cortical_Meta_Analysis_Rank'))
# sheets_write(fullTableRank, mar, sheet = "Full_Meta_Analysis_Rank")
# sheets_write(femaleTableRank, mar, sheet = 'Female_Meta_Analysis_Rank')
# sheets_write(maleTableRank, mar,'Male_Meta_Analysis_Rank')
# sheets_write(corticalTableRank, mar, 'Cortical_Meta_Analysis_Rank')
# 
# drive_mv(file = "Meta_Analysis_Rank", path = "~/Thesis/Manuscript/Tables/Full_Tables/")  # move Sheets file
# 
# 
# maf <- drive_get("Thesis/Manuscript/Tables/Full_Tables/Sex_Interaction_Meta_Analysis")
# if(nrow(maf) != 0) {
#   drive_rm(maf)
# }
# #create the google worksheet
# maf <- sheets_create("Sex_Interaction_Meta_Analysis", sheets = c('Full_Meta_Analysis', 'Cortical_Meta_Analysis'))
# sheets_write(fullTable_Flip, maf, sheet = "Full_Meta_Analysis")
# sheets_write(corticalTable_Flip, maf, 'Cortical_Meta_Analysis')
# 
# drive_mv(file = "Sex_Interaction_Meta_Analysis", path = "~/Thesis/Manuscript/Tables/Full_Tables/")  # move Sheets file
# 
# 
# marf <- drive_get("Thesis/Manuscript/Tables/Full_Tables/Sex_Interaction_Meta_Analysis_Rank")
# if(nrow(marf) != 0) {
#   drive_rm(marf)
# }
# #create the google worksheet
# marf <- sheets_create("Sex_Interaction_Meta_Analysis_Rank", sheets = c('Full_Meta_Analysis_Rank', 'Cortical_Meta_Analysis_Rank'))
# sheets_write(fullTableRank_Flip, marf, sheet = "Full_Meta_Analysis_Rank")
# sheets_write(corticalTableRank_Flip, marf, 'Cortical_Meta_Analysis_Rank')
# 
# drive_mv(file = "Sex_Interaction_Meta_Analysis_Rank", path = "~/Thesis/Manuscript/Tables/Full_Tables/")  # move Sheets file

##################################################################################
##################################SLIM SHEETS
##################################################################################
ma <- drive_get("~/Thesis/Manuscript/Tables/Slim_Tables/Official_Meta_Analysis")
if(nrow(ma) != 0) {
  drive_rm(ma)
}
#create the google worksheet
ma <- sheets_create("Official_Meta_Analysis", sheets = c('Full_Meta_Analysis', 'Male_Meta_Analysis','Female_Meta_Analysis','Cortical_Meta_Analysis'))
sheets_write(fullTable %>% select(gene_symbol, Howard.pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_p, Bonferroni_Correction, slim_region_location, cell_type_taxon, DE_Prior_Rank), ma, sheet = "Full_Meta_Analysis")
sheets_write(femaleTable %>% select(gene_symbol, Howard.pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_p, Bonferroni_Correction, slim_region_location, cell_type_taxon, DE_Prior_Rank), ma, sheet = 'Female_Meta_Analysis')
sheets_write(maleTable %>% select(gene_symbol, Howard.pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_p, Bonferroni_Correction, slim_region_location, cell_type_taxon, DE_Prior_Rank), ma,'Male_Meta_Analysis')
sheets_write(corticalTable %>% select(gene_symbol, Howard.pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_p, Bonferroni_Correction, slim_region_location, cell_type_taxon, DE_Prior_Rank), ma, 'Cortical_Meta_Analysis')

drive_mv(file = "Official_Meta_Analysis", path = "~/Thesis/Manuscript/Tables/Slim_Tables/")  # move Sheets file

mar <- drive_get("Thesis/Manuscript/Tables/Slim_Tables/Official_Meta_Analysis_Rank")
if(nrow(mar) != 0) {
  drive_rm(mar)
}
#create the google worksheet
mar <- sheets_create("Official_Meta_Analysis_Rank", sheets = c('Full_Meta_Analysis_Rank', 'Male_Meta_Analysis_Rank','Female_Meta_Analysis_Rank', 'Cortical_Meta_Analysis_Rank'))
sheets_write(fullTableRank  %>% select(gene_symbol, Howard.pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Corrected_meta_empirical_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), mar, sheet = "Full_Meta_Analysis_Rank")
sheets_write(femaleTableRank  %>% select(gene_symbol, Howard.pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Corrected_meta_empirical_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), mar, sheet = 'Female_Meta_Analysis_Rank')
sheets_write(maleTableRank  %>% select(gene_symbol, Howard.pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Corrected_meta_empirical_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), mar,'Male_Meta_Analysis_Rank')
sheets_write(corticalTableRank  %>% select(gene_symbol, Howard.pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Corrected_meta_empirical_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), mar, 'Cortical_Meta_Analysis_Rank')

drive_mv(file = "Official_Meta_Analysis_Rank", path = "~/Thesis/Manuscript/Tables/Slim_Tables/")  # move Sheets file


maf <- drive_get("Thesis/Manuscript/Tables/Slim_Tables/Official_Sex_Interaction_Meta_Analysis")
if(nrow(maf) != 0) {
  drive_rm(maf)
}
#create the google worksheet
maf <- sheets_create("Official_Sex_Interaction_Meta_Analysis", sheets = c('Full_Meta_Analysis', 'Cortical_Meta_Analysis'))
sheets_write(fullTable_Flip  %>% select(gene_symbol, Howard.pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_p, Bonferroni_Correction, slim_region_location, cell_type_taxon, DE_Prior_Rank), maf, sheet = "Full_Meta_Analysis")
sheets_write(corticalTable_Flip  %>% select(gene_symbol, Howard.pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_p, Bonferroni_Correction, slim_region_location, cell_type_taxon, DE_Prior_Rank), maf, 'Cortical_Meta_Analysis')

drive_mv(file = "Official_Sex_Interaction_Meta_Analysis", path = "~/Thesis/Manuscript/Tables/Slim_Tables/")  # move Sheets file


marf <- drive_get("Thesis/Manuscript/Tables/Slim_Tables/Official_Sex_Interaction_Meta_Analysis_Rank")
if(nrow(marf) != 0) {
  drive_rm(marf)
}
#create the google worksheet
marf <- sheets_create("Official_Sex_Interaction_Meta_Analysis_Rank", sheets = c('Full_Meta_Analysis_Rank', 'Cortical_Meta_Analysis_Rank'))
sheets_write(fullTableRank_Flip  %>% select(gene_symbol, Howard.pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Corrected_meta_empirical_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), marf, sheet = "Full_Meta_Analysis_Rank")
sheets_write(corticalTableRank_Flip  %>% select(gene_symbol, Howard.pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Corrected_meta_empirical_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), marf, 'Cortical_Meta_Analysis_Rank')

drive_mv(file = "Official_Sex_Interaction_Meta_Analysis_Rank", path = "~/Thesis/Manuscript/Tables/Slim_Tables/")  # move Sheets file



