library(readr)
library(magrittr)
library(dplyr)
library(tidyr)
library(here)
library(googledrive)
library(googlesheets4)

source(here("R", "loading datasets","mergeAnalysis.R"))
source(here("R", "loading datasets","ZeiselPolygenic.R"))

fullLabonte <- read_csv( here("ProcessedData", "LabonteEtAl", "fullLabonteTable_magma.csv"))
fullLabonte_flipped <- read_csv(here("ProcessedData", "LabonteEtAl", "fullLabonteTableMagma_flipped.csv"))
femaleLabonte <- read_csv(here("ProcessedData", "LabonteEtAl", "FemaleLabonteTable_magma.csv")) 
maleLabonte <- read_csv(here("ProcessedData", "LabonteEtAl", "MaleLabonteTable_magma.csv"))
maleLabonte_flipped <- read_csv(here("ProcessedData", "LabonteEtAl", "MaleLabonteTableMagma_flipped.csv"))

corticalLabonte <- read_csv( here("ProcessedData", "LabonteEtAl", "CorticalLabonteTableMagma.csv"))
corticalLabonte_flipped <- read_csv( here("ProcessedData", "LabonteEtAl", "CorticalLabonteTableMagma_flipped.csv"))
femaleCorticalLabonte <- read_csv(here("ProcessedData", "LabonteEtAl", "CorticalFemaleLabonteTableMagma.csv"))
maleCorticalLabonte <- read_csv(here("ProcessedData", "LabonteEtAl", "CorticalMaleLabonteTableMagma.csv"))
maleCorticalLabonte_flipped <- read_csv(here("ProcessedData", "LabonteEtAl", "CorticalMaleLabonteTableMagma_flipped.csv"))

fullDing <- read_csv(here("ProcessedData", "DingEtAl", "DingTableFisher.slim_magma.csv"))
fullDing_flipped <- read_csv(here("ProcessedData", "DingEtAl", "DingTableFisherMagma.slim_flipped.csv"))
femaleDing <-read_csv(here("ProcessedData", "DingEtAl", "DingTableFisher.slim.Female_magma.csv"))
maleDing <-read_csv(here("ProcessedData", "DingEtAl", "DingTableFisher.slim.Male_magma.csv"))
maleDing_flipped <-read_csv(here("ProcessedData", "DingEtAl", "DingTableFisherMagma.slim.Male_flipped.csv"))

corticalDing <- read_csv( here("ProcessedData", "DingEtAl", "CorticalDingTableFisherMagma.slim.csv"))
corticalDing_flipped <- read_csv( here("ProcessedData", "DingEtAl", "CorticalDingTableFisherMagma.slim_flipped.csv"))
femaleCorticalDing <- read_csv(here("ProcessedData", "DingEtAl", "CorticalDingTableFisherMagma.slim.Female.csv"))
maleCorticalDing <- read_csv(here("ProcessedData", "DingEtAl", "CorticalDingTableFisherMagma.slim.Male.csv"))
maleCorticalDing_flipped <- read_csv(here("ProcessedData", "DingEtAl", "CorticalDingTableFisherMagma.slim.Male_flipped.csv"))

fullRamaker <- read_csv(here("ProcessedData", "RamakerEtAl", "fullRamakerTable_magma.csv"))
fullRamaker_flipped <- read_csv(here("ProcessedData", "RamakerEtAl", "fullRamakerTableMagma_flipped.csv"))
femaleRamaker <- read_csv(here("ProcessedData", "RamakerEtAl", "FemaleRamakerTable_magma.csv"))
maleRamaker <- read_csv(here("ProcessedData", "RamakerEtAl", "MaleRamakerTable_magma.csv"))
maleRamaker_flipped <- read_csv(here("ProcessedData", "RamakerEtAl", "MaleRamakerTableMagma_flipped.csv"))

corticalRamaker <- read_csv(here("ProcessedData", "RamakerEtAl", "fullCorticalRamakerTable.csv"))
corticalRamaker_flipped <- read_csv(here("ProcessedData", "RamakerEtAl", "fullCorticalRamakerTable_flipped.csv"))
femaleCorticalRamaker <- read_csv(here("ProcessedData", "RamakerEtAl", "fullCorticalFemaleRamakerTable.csv"))
maleCorticalRamaker <- read_csv(here("ProcessedData", "RamakerEtAl", "fullCorticalMaleRamakerTable.csv"))
maleCorticalRamaker_flipped <- read_csv(here("ProcessedData", "RamakerEtAl", "fullCorticalMaleRamakerTable_flipped.csv"))

Howard <- read_csv(here("ProcessedData", "HowardEtAl", "fullHowardTable.csv"))
# Howard %<>% dplyr::rename(gene_symbol = Howard.Genes)
DE_Prior <- read_tsv(here("data", "CrowEtAl","pnas.1802973116.sd02.txt"))
Howard_Polygenics <- read_csv(here("ProcessedData","HowardEtAl", "HowardRegionsPolygenicCellTypes_four.csv"))


fullLabonteDir <- "female_Anterior_Insula_female_BA11_female_BA25_female_BA8/9_female_Nac_female_Subic_male_Anterior_Insula_male_BA11_male_BA25_male_BA8/9_male_Nac_male_Subic"
femaleLabonteDir <- "female_Anterior_Insula_female_BA11_female_BA25_female_BA8/9_female_Nac_female_Subic"
maleLabonteDir <- "male_Anterior_Insula_male_BA11_male_BA25_male_BA8/9_male_Nac_male_Subic"

corticalLabonteDir <- "female_Anterior_Insula_female_BA11_female_BA25_female_BA8/9_male_Anterior_Insula_male_BA11_male_BA25_male_BA8/9"
coritcalFemaleLabonteDir <-"female_Anterior_Insula_female_BA11_female_BA25_female_BA8/9"
corticalMaleLabonteDir <-"male_Anterior_Insula_male_BA11_male_BA25_male_BA8/9"

fullDingDir <- "female_ACC.1_female_ACC.2_female_AMY_female_DLPFC_male_ACC.1_male_ACC.2_male_AMY_male_DLPFC"
femaleDingDir <-"female_ACC.1_female_ACC.2_female_AMY_female_DLPFC"
maleDingDir <-"male_ACC.1_male_ACC.2_male_AMY_male_DLPFC"

corticalDingDir <-"female_ACC.1_female_ACC.2_female_DLPFC_male_ACC.1_male_ACC.2_male_DLPFC"
corticalFemaleDingDir <-"female_ACC.1_female_ACC.2_female_DLPFC"
corticalMaleDingDir<-"male_ACC.1_male_ACC.2_male_DLPFC"

fullRamakerDir <- "AnCg.F_nAcc.F_DLPFC.F_AnCg.M_nAcc.M_DLPFC.M"
femaleRamakerDir <- "AnCg_nAcc_DLPFC_Female_directions"
maleRamakerDir <- "AnCg_nAcc_DLPFC_Male_directions"

corticalRamakerDir <- "AnCg.F_DLPFC.F_AnCg.M_DLPFC.M"
corticalFemaleRamakerDir <- "AnCg_DLPFC_Female_directions"
corticalMaleRamakerDir <- "AnCg_DLPFC_Male_directions"

##### MAGMA ANALYSIS ON FULL GENOME Wilcox Test ####--------------------------------------------------
magma <- read_csv(here("data", "HowardEtAl", "FullMagmaGenes.csv"))

#full analysis of meta analysis on magma genes
mergeMetaMagma(magma, fullLabonte, fullLabonteDir, fullDing, fullDingDir, fullRamaker, fullRamakerDir) %>% write_csv(here("ProcessedData","fullTableMetaMagma.csv"))
mergeMetaMagmaRank(magma, fullLabonte, fullLabonteDir, fullDing, fullDingDir, fullRamaker, fullRamakerDir) %>% write_csv(path = here("ProcessedData", "fullTableMetaMagmaRank.csv"))
mergeMetaMagma(magma, fullLabonte_flipped, fullLabonteDir, fullDing_flipped, fullDingDir, fullRamaker_flipped, fullRamakerDir) %>% write_csv(path = here("ProcessedData", "SIFullTableMetaMagma.csv"))
mergeMetaMagmaRank(magma, fullLabonte_flipped, fullLabonteDir, fullDing_flipped, fullDingDir, fullRamaker_flipped, fullRamakerDir) %>% write_csv(path = here("ProcessedData", "SIFullTableMetaMagmaRank.csv"))
mergeMetaMagma(magma, femaleLabonte, femaleLabonteDir, femaleDing, femaleDingDir, femaleRamaker, femaleRamakerDir) %>% write_csv(path = here("ProcessedData", "femaleTableMetaMagma.csv"))
mergeMetaMagmaRank(magma, femaleLabonte, femaleLabonteDir, femaleDing, femaleDingDir, femaleRamaker, femaleRamakerDir) %>% write_csv(path = here("ProcessedData", "femaleTableMetaMagmaRank.csv"))
mergeMetaMagma(magma, maleLabonte, maleLabonteDir, maleDing, maleDingDir, maleRamaker, maleRamakerDir) %>% write_csv(path = here("ProcessedData", "maleTableMetaMagma.csv"))
mergeMetaMagmaRank(magma, maleLabonte, maleLabonteDir, maleDing, maleDingDir, maleRamaker, maleRamakerDir) %>% write_csv(path = here("ProcessedData", "maleTableMetaMagmaRank.csv"))
mergeMetaMagma(magma, corticalLabonte, corticalLabonteDir, corticalDing, corticalDingDir, corticalRamaker, corticalRamakerDir) %>% write_csv(path = here("ProcessedData", "corticalTableMetaMagma.csv"))
mergeMetaMagmaRank(magma, corticalLabonte, corticalLabonteDir, corticalDing, corticalDingDir, corticalRamaker, corticalRamakerDir) %>% write_csv(path = here("ProcessedData", "corticalTableMetaMagmaRank.csv"))
mergeMetaMagma(magma, corticalLabonte_flipped, corticalLabonteDir, corticalDing_flipped, corticalDingDir, corticalRamaker_flipped, corticalRamakerDir) %>% write_csv(path = here("ProcessedData", "SICorticalTableMetaMagma.csv"))
mergeMetaMagmaRank(magma, corticalLabonte_flipped, corticalLabonteDir, corticalDing_flipped, corticalDingDir, corticalRamaker_flipped, corticalRamakerDir) %>% write_csv(path = here("ProcessedData", "SICorticalTableMetaMagmaRank.csv"))
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
 
