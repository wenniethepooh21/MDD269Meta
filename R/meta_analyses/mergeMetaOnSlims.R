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
corticalLabonte_flipped <- read_csv( here("Processed_Data/LabonteEtAl/CorticalLabonteTableMagma.csv"))

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
source(here("R", "loading datasets","mergeAnalysis.R"))

source(here("R", "loading datasets","ZeiselPolygenic.R"))

##### MAGMA ANALYSIS ON FULL GENOME Wilcox Test ####--------------------------------------------------
magma <- read_csv(here("data", "HowardEtAl", "FullMagmaGenes.csv"))

#full analysis of meta analysis on magma genes
mergeMagmaMetaRank(magma, fullLabonte, fullLabonteDir, fullDing, fullDingDir, fullRamaker, fullRamakerDir,"meta") %>% write_csv(here("Processed_Data/Meta_Results/MAGMA/fullTableMetaMagma.csv"))
mergeMagmaMetaRank(magma, fullLabonte, fullLabonteDir, fullDing, fullDingDir, fullRamaker, fullRamakerDir,"rank") %>% write_csv(path = here("Processed_Data/Meta_Results/MAGMA/fullTableMetaMagmaRank.csv"))
mergeMagmaMetaRank(magma, fullLabonte_flipped, fullLabonteDir, fullDing_flipped, fullDingDir, fullRamaker_flipped, fullRamakerDir,"meta") %>% write_csv(path = here("Processed_Data/Meta_Results/MAGMA/SIFullTableMetaMagma.csv"))
mergeMagmaMetaRank(magma, fullLabonte_flipped, fullLabonteDir, fullDing_flipped, fullDingDir, fullRamaker_flipped, fullRamakerDir,"rank") %>% write_csv(path = here("Processed_Data/Meta_Results/MAGMA/SIFullTableMetaMagmaRank.csv"))
mergeMagmaMetaRank(magma, femaleLabonte, femaleLabonteDir, femaleDing, femaleDingDir, femaleRamaker, femaleRamakerDir,"meta") %>% write_csv(path = here("Processed_Data/Meta_Results/MAGMA/femaleTableMetaMagma.csv"))
mergeMagmaMetaRank(magma, femaleLabonte, femaleLabonteDir, femaleDing, femaleDingDir, femaleRamaker, femaleRamakerDir,"rank") %>% write_csv(path = here("Processed_Data/Meta_Results/MAGMA/femaleTableMetaMagmaRank.csv"))
mergeMagmaMetaRank(magma, maleLabonte, maleLabonteDir, maleDing, maleDingDir, maleRamaker, maleRamakerDir,"meta") %>% write_csv(path = here("Processed_Data/Meta_Results/MAGMA/maleTableMetaMagma.csv"))
mergeMagmaMetaRank(magma, maleLabonte, maleLabonteDir, maleDing, maleDingDir, maleRamaker, maleRamakerDir,"rank") %>% write_csv(path = here("Processed_Data/Meta_Results/MAGMA/maleTableMetaMagmaRank.csv"))
mergeMagmaMetaRank(magma, corticalLabonte, corticalLabonteDir, corticalDing, corticalDingDir, corticalRamaker, corticalRamakerDir,"meta") %>% write_csv(path = here("Processed_Data/Meta_Results/MAGMA/corticalTableMetaMagma.csv"))
mergeMagmaMetaRank(magma, corticalLabonte, corticalLabonteDir, corticalDing, corticalDingDir, corticalRamaker, corticalRamakerDir,"rank") %>% write_csv(path = here("Processed_Data/Meta_Results/MAGMA/corticalTableMetaMagmaRank.csv"))
mergeMagmaMetaRank(magma, corticalLabonte_flipped, corticalLabonteDir, corticalDing_flipped, corticalDingDir, corticalRamaker_flipped, corticalRamakerDir,"meta") %>% write_csv(path = here("Processed_Data/Meta_Results/MAGMA/SICorticalTableMetaMagma.csv"))
mergeMagmaMetaRank(magma, corticalLabonte_flipped, corticalLabonteDir, corticalDing_flipped, corticalDingDir, corticalRamaker_flipped, corticalRamakerDir,"rank") %>% write_csv(path = here("Processed_Data/Meta_Results/MAGMA/SICorticalTableMetaMagmaRank.csv"))
#-----------------------------------------------------------------------------------------

#full analysis
fullTable <- mergeMeta(Howard, fullLabonte, fullLabonteDir, fullDing, fullDingDir, fullRamaker, fullRamakerDir)
fullTable %<>% MetaAnalysis()
fullTable %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
#for each gene in Howard, find the top cell type cluster from zeisel data
fullTable %<>% left_join(Howard_Polygenics %>% select(-Updated_Gene_Names, -mouseGene, -Howard.Genes_upper), by = c('gene_symbol' = 'Howard.Genes')) %>% write_csv(path = here("ProcessedData", "fullTableMeta.csv")) #write out for wilkoxon test
fullTable %<>% mutate_all(replace_na, replace = "Not_Available")

fullTableRank <- mergeMeta(Howard, fullLabonte, fullLabonteDir, fullDing, fullDingDir, fullRamaker, fullRamakerDir)
fullTableRank %<>% GenomeRank()
fullTableRank %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name')) 
#for each gene in Howard, find the top cell type cluster from zeisel data
fullTableRank %<>% left_join(Howard_Polygenics %>% select(-Updated_Gene_Names, -mouseGene,-Howard.Genes_upper), by = c('gene_symbol' = 'Howard.Genes'))
fullTableRank %<>% mutate_all(replace_na, replace = "Not_Available")


fullTable_Flip <- mergeMeta(Howard, fullLabonte_flipped, fullLabonteDir, fullDing_flipped, fullDingDir, fullRamaker_flipped, fullRamakerDir)
fullTable_Flip %<>% MetaAnalysis()
fullTable_Flip %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
#for each gene in Howard, find the top cell type cluster from zeisel data
fullTable_Flip %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene,-Howard.Genes_upper), by = c('gene_symbol' = 'Howard.Genes'))
fullTable_Flip %<>% mutate_all(replace_na, replace = "Not_Available")


fullTableRank_Flip <- mergeMeta(Howard, fullLabonte_flipped, fullLabonteDir, fullDing_flipped, fullDingDir, fullRamaker_flipped, fullRamakerDir)
fullTableRank_Flip %<>% GenomeRank()
fullTableRank_Flip %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
#for each gene in Howard, find the top cell type cluster from zeisel data
fullTableRank_Flip %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene,-Howard.Genes_upper), by = c('gene_symbol' = 'Howard.Genes'))
fullTableRank_Flip %<>% mutate_all(replace_na, replace = "Not_Available")


#Female analysis
femaleTable <- mergeMeta(Howard,femaleLabonte, femaleLabonteDir, femaleDing, femaleDingDir, femaleRamaker, femaleRamakerDir)
femaleTable %<>% MetaAnalysis()
femaleTable %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
#for each gene in Howard, find the top cell type cluster from zeisel data
femaleTable %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene,-Howard.Genes_upper), by = c('gene_symbol' = 'Howard.Genes'))
femaleTable %<>% mutate_all(replace_na, replace = "Not_Available")

femaleTableRank <- mergeMeta(Howard,femaleLabonte, femaleLabonteDir, femaleDing, femaleDingDir, femaleRamaker, femaleRamakerDir)
femaleTableRank %<>% GenomeRank()
femaleTableRank %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
#for each gene in Howard, find the top cell type cluster from zeisel data
femaleTableRank %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene,-Howard.Genes_upper), by = c('gene_symbol' = 'Howard.Genes'))
femaleTableRank %<>% mutate_all(replace_na, replace = "Not_Available")

#male analysis
maleTable <- mergeMeta(Howard, maleLabonte, maleLabonteDir,maleDing, maleDingDir, maleRamaker, maleRamakerDir)
maleTable %<>% MetaAnalysis()
maleTable %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
#for each gene in Howard, find the top cell type cluster from zeisel data
maleTable %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene,-Howard.Genes_upper), by = c('gene_symbol' = 'Howard.Genes'))
maleTable %<>% mutate_all(replace_na, replace = "Not_Available")


maleTableRank <- mergeMeta(Howard, maleLabonte, maleLabonteDir,maleDing, maleDingDir, maleRamaker, maleRamakerDir)
maleTableRank %<>% GenomeRank()
maleTableRank %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
#for each gene in Howard, find the top cell type cluster from zeisel data
maleTableRank %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene,-Howard.Genes_upper), by = c('gene_symbol' = 'Howard.Genes'))
maleTableRank %<>% mutate_all(replace_na, replace = "Not_Available")


# maleTable_Flip <- mergeMeta(Howard, MaleLabonte_flipped, MaleLabonteDir, MaleDing_flipped, MaleDingDir, MaleRamaker_flipped, maleRamakerDir)
# maleTable_Flip %<>% MetaAnalysis()
# maleTable_Flip %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
# #for each gene in Howard, find the top cell type cluster from zeisel data
# maleTable_Flip %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene,-Howard.Genes_upper), by = c('gene_symbol' = 'Howard.Genes'))
# 
# maleTableRank_Flip <- mergeMeta(Howard, MaleLabonte_flipped, MaleLabonteDir,MaleDing_flipped, MaleDingDir, MaleRamaker_flipped, maleRamakerDir)
# maleTableRank_Flip %<>% GenomeRank()
# maleTableRank_Flip %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
# #for each gene in Howard, find the top cell type cluster from zeisel data
# maleTableRank_Flip %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'Howard.Genes'))

# #cortical analysis
corticalTable <- mergeMeta(Howard, corticalLabonte, corticalLabonteDir, corticalDing, corticalDingDir, corticalRamaker, corticalRamakerDir)
corticalTable %<>% MetaAnalysis()
corticalTable %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name')) 
#for each gene in Howard, find the top cell type cluster from zeisel data
corticalTable %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene,-Howard.Genes_upper), by = c('gene_symbol' = 'Howard.Genes'))
corticalTable %<>% mutate_all(replace_na, replace = "Gene_Not_Detected")

corticalTableRank <- mergeMeta(Howard, corticalLabonte, corticalLabonteDir, corticalDing, corticalDingDir, corticalRamaker, corticalRamakerDir)
corticalTableRank %<>% GenomeRank()
corticalTableRank %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
#for each gene in Howard, find the top cell type cluster from zeisel data
corticalTableRank %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene,-Howard.Genes_upper), by = c('gene_symbol' = 'Howard.Genes'))
corticalTableRank %<>% mutate_all(replace_na, replace = "Gene_Not_Detected")

corticalTable_Flip <- mergeMeta(Howard, corticalLabonte_flipped, corticalLabonteDir, corticalDing_flipped, corticalDingDir, corticalRamaker_flipped, corticalRamakerDir)
corticalTable_Flip %<>% MetaAnalysis()
corticalTable_Flip %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name')) 
#for each gene in Howard, find the top cell type cluster from zeisel data
corticalTable_Flip %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene,-Howard.Genes_upper), by = c('gene_symbol' = 'Howard.Genes'))
corticalTable_Flip %<>% mutate_all(replace_na, replace = "Gene_Not_Detected")

corticalTableRank_Flip <- mergeMeta(Howard, corticalLabonte_flipped, corticalLabonteDir, corticalDing_flipped, corticalDingDir, corticalRamaker_flipped, corticalRamakerDir)
corticalTableRank_Flip %<>% GenomeRank()
corticalTableRank_Flip %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
#for each gene in Howard, find the top cell type cluster from zeisel data
corticalTableRank_Flip %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene,-Howard.Genes_upper), by = c('gene_symbol' = 'Howard.Genes'))
corticalTableRank_Flip %<>% mutate_all(replace_na, replace = "Gene_Not_Detected")

# #femalecortical
# femaleCortical <- mergeMeta(Howard, femaleCorticalLabonte, coritcalFemaleLabonteDir, femaleCorticalDing, corticalFemaleDingDir,femaleCorticalRamaker, corticalFemaleRamakerDir)
# femaleCortical %<>% MetaAnalysis()
# femaleCortical %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
# #for each gene in Howard, find the top cell type cluster from zeisel data
# femaleCortical %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'Howard.Genes'))
# # 
# femaleCorticalRank <- mergeMeta(Howard, femaleCorticalLabonte, coritcalFemaleLabonteDir, femaleCorticalDing, corticalFemaleDingDir,femaleCorticalRamaker, corticalFemaleRamakerDir)
# femaleCorticalRank %<>% GenomeRank()
# femaleCorticalRank %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
# #for each gene in Howard, find the top cell type cluster from zeisel data
# femaleCorticalRank %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'Howard.Genes'))
# # 
# #malecortical
# maleCortical <- mergeMeta(Howard, maleCorticalLabonte, corticalMaleLabonteDir,maleCorticalDing, corticalMaleDingDir, maleCorticalRamaker, corticalMaleRamakerDir)
# maleCortical %<>% MetaAnalysis()
# maleCortical %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
# #for each gene in Howard, find the top cell type cluster from zeisel data
# maleCortical %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'Howard.Genes'))
# 
# maleCorticalRank <- mergeMeta(Howard, maleCorticalLabonte, corticalMaleLabonteDir,maleCorticalDing, corticalMaleDingDir, maleCorticalRamaker, corticalMaleRamakerDir)
# maleCorticalRank %<>% GenomeRank()
# maleCorticalRank %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
# #for each gene in Howard, find the top cell type cluster from zeisel data
# maleCorticalRank %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'Howard.Genes'))
# 
# maleCortical_Flip <- mergeMeta(Howard, maleCorticalLabonte_flipped, corticalMaleLabonteDir,maleCorticalDing_flipped, corticalMaleDingDir, maleCorticalRamaker_flipped, corticalMaleRamakerDir)
# maleCortical_Flip %<>% MetaAnalysis()
# maleCortical_Flip %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
# #for each gene in Howard, find the top cell type cluster from zeisel data
# maleCortical_Flip %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'Howard.Genes'))
# 
# 
# maleCorticalRank_Flip <- mergeMeta(Howard, maleCorticalLabonte_flipped, corticalMaleLabonteDir,maleCorticalDing_flipped, corticalMaleDingDir, maleCorticalRamaker_flipped, corticalMaleRamakerDir)
# maleCorticalRank_Flip %<>% GenomeRank()
# maleCorticalRank_Flip %<>% left_join(DE_Prior %>% select(Gene_Name, DE_Prior_Rank), by = c('gene_symbol' = 'Gene_Name'))
# #for each gene in Howard, find the top cell type cluster from zeisel data
# maleCorticalRank_Flip %<>% left_join(Howard_Polygenics%>% select(-Updated_Gene_Names, -mouseGene), by = c('gene_symbol' = 'Howard.Genes'))
# 
# 
#Access googlesheets
drive_auth()
sheets_auth(token = drive_token())

ma <- drive_get("Meta-Analysis")
if(nrow(ma) != 0) {
  drive_rm(ma)
}
#create the google worksheet
ma <- sheets_create("Meta-Analysis", sheets = c('Full_Meta_Analysis', 'Male_Meta_Analysis','Female_Meta_Analysis','Cortical_Meta_Analysis'))
sheets_write(fullTable, ma, sheet = "Full_Meta_Analysis")
sheets_write(femaleTable, ma, sheet = 'Female_Meta_Analysis')
sheets_write(maleTable, ma,'Male_Meta_Analysis')
sheets_write(corticalTable, ma, 'Cortical_Meta_Analysis')

mar <- drive_get("Meta-Analysis-Rank")
if(nrow(mar) != 0) {
  drive_rm(mar)
}
#create the google worksheet
mar <- sheets_create("Meta-Analysis-Rank", sheets = c('Full_Meta_Analysis_Rank', 'Male_Meta_Analysis_Rank','Female_Meta_Analysis_Rank', 'Cortical_Meta_Analysis_Rank'))
sheets_write(fullTableRank, mar, sheet = "Full_Meta_Analysis_Rank")
sheets_write(femaleTableRank, mar, sheet = 'Female_Meta_Analysis_Rank')
sheets_write(maleTableRank, mar,'Male_Meta_Analysis_Rank')
sheets_write(corticalTableRank, mar, 'Cortical_Meta_Analysis_Rank')

maf <- drive_get("Gender-Interaction-Meta-Analysis")
if(nrow(maf) != 0) {
  drive_rm(maf)
}
#create the google worksheet
maf <- sheets_create("Sex-Interaction-Meta-Analysis", sheets = c('Full_Meta_Analysis', 'Cortical_Meta_Analysis'))
sheets_write(fullTable_Flip, maf, sheet = "Full_Meta_Analysis")
sheets_write(corticalTable_Flip, maf, 'Cortical_Meta_Analysis')

marf <- drive_get("Sex-Interaction-Meta-Analysis-Rank")
if(nrow(marf) != 0) {
  drive_rm(marf)
}
#create the google worksheet
marf <- sheets_create("Sex-Interaction-Meta-Analysis-Rank", sheets = c('Full_Meta_Analysis_Rank', 'Cortical_Meta_Analysis_Rank'))
sheets_write(fullTableRank_Flip, marf, sheet = "Full_Meta_Analysis_Rank")
sheets_write(corticalTableRank_Flip, marf, 'Cortical_Meta_Analysis_Rank')
 
