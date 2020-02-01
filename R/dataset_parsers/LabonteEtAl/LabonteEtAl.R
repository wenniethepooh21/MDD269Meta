library(magrittr)
library(tidyverse)
library(here)
#library(biomaRt)
library(dplyr)
library(readxl)

source(here("R", "loading datasets", "LabonteMetaPAnalysis.R"))

#Import male data
LabonteMale <- read_xlsx(here("data", "LabonteEtAl", "Male DEG all.xlsx"))
LabonteMale %<>% rename(symbol  = "Gene name", brain_region = "Brain Region", Male.pvalue = "p-value", Male.logFC = logFC, Male.AveExpr = AveExpr)

#Import female data
#LabonteFemale <- read_csv(here("data","LabonteEtAl","nm.4386-S4.csv"))
LabonteFemale <- read_xlsx(here("data", "LabonteEtAl", "Female DEG all.xlsx"))
LabonteFemale %<>% rename(brain_region = `Brain region`, symbol = `Gene Name`, Female.pvalue = `p-value`, Female.logFC = logFC, Female.AveExpr = AveExpr)

Labonte <- left_join(LabonteMale %>% select(-Biotype, -Description), LabonteFemale %>% select(-Biotype, -Description), by = c("brain_region", "ENSG", "symbol"))

########################################### update gene symbols to remove excel dates from data#########
#convert ENSG to gene symbol to remove dates from symbols column
#ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

#This line was previously run to generate associated gene name, output was saved
#ENSG_Genes <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol'), filters = 'ensembl_gene_id', values = unique(Labonte$ENSG), mart = ensembl) %>% as_tibble() 
#ENSG_Genes %>% write_csv(here("ProcessedData", "LabonteEtAl", "getBMGenes.csv"))

ENSG_Genes <- read_csv(here("ProcessedData", "LabonteEtAl", "getBMGenes.csv"))
Labonte_genes <- dplyr::select(Labonte, ENSG, symbol)
Labonte_genes %<>% left_join(ENSG_Genes, by = c('ENSG' = 'ensembl_gene_id'))

#merge gene symbols into one column 
na_index <- which(is.na(Labonte_genes$hgnc_symbol))
Labonte_genes$hgnc_symbol[na_index] <- Labonte_genes$symbol[na_index]
Labonte_genes %<>% dplyr::select(-symbol)	
Labonte_genes <- Labonte_genes[!duplicated(Labonte_genes$ENSG), ]
Labonte %<>% left_join(Labonte_genes, by = c('ENSG' = 'ENSG'))
Labonte %<>% dplyr::select(-symbol)
Labonte %<>% rename(symbol = "hgnc_symbol")
############################################################################################
#re-order 
Labonte <- Labonte[,c(1:2,9,3:8)]
Labonte %<>% mutate(brain_region = ifelse(brain_region == "Anterior Insula", "Anterior_Insula",brain_region))

#Write Data to .csv files
Labonte %>% write_csv(path = here("ProcessedData", "LabonteEtAl","CompleteLabonteTable.csv"))

Labonte_long <- bind_rows(Labonte %>% dplyr::select(brain_region, symbol, logFC = Male.logFC, pvalue = Male.pvalue) %>% mutate(sex = "male"),
                          Labonte %>% dplyr::select(brain_region, symbol, logFC = Female.logFC, pvalue = Female.pvalue) %>% mutate(sex = "female"))

#combine 12 (6 regions * 2 sexes) + 12 pvalues for each gene
Labonte_summary_results <- LabonteMeta(Labonte_long,12,6)

#get the sex data separately
Labonte_Male <- Labonte %>% select(brain_region, symbol, logFC = Male.logFC, pvalue = Male.pvalue) %>% mutate(sex = "male")
Labonte_Female <- Labonte %>% select(brain_region, symbol, logFC = Female.logFC, pvalue = Female.pvalue) %>% mutate(sex = "female")

#combine 6 (6 regions * 1 sexes) 
Labonte_Female_results <- LabonteMeta(Labonte_Female,6,6)

#combine 6 (6 regions * 1 sexes) 
Labonte_Male_results <- LabonteMeta(Labonte_Male,6,6)

#needs writing out
Labonte_summary_results %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "fullLabonteTable.csv"))
Labonte_Female_results %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "FemaleLabonteTable.csv"))
Labonte_Male_results %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "MaleLabonteTable.csv"))


# #------- Genome ranking (MAGMA list)
magma <- read_csv(here("data", "HowardEtAl", "FullMagmaGenes.csv"))
full_Labonte <- read_csv( here("ProcessedData", "LabonteEtAl", "CompleteLabonteTable.csv")) %>% na.omit()
full_Labonte %<>% mutate(symbol = gsub("C([X0-9]+)ORF([0-9]+)", "C\\1orf\\2", symbol))

Labonte_magma <- left_join(magma %>% select(Labonte_genes) %>% distinct(), full_Labonte, by = c('Labonte_genes' = 'symbol')) %>% na.omit() #remove the one NA value
Labonte_magma %<>% rename(symbol = Labonte_genes) %>% write_csv(here("ProcessedData", "LabonteEtAl", "CompleteLabonteTableMagma.csv"))

Labonte_long_magma <- bind_rows(Labonte_magma %>% dplyr::select(brain_region, symbol, logFC = Male.logFC, pvalue = Male.pvalue) %>% mutate(sex = "male"),
                          Labonte_magma %>% dplyr::select(brain_region, ENSG, symbol, logFC = Female.logFC, pvalue = Female.pvalue) %>% mutate(sex = "female"))

#combine 12 (6 regions * 2 sexes) + 12 pvalues for each gene
Labonte_summary_magma <- LabonteMeta(Labonte_long_magma,12,6)

#get the sex data separately
Labonte_Male_magma <- Labonte_magma %>% select(brain_region, symbol, logFC = Male.logFC, pvalue = Male.pvalue) %>% mutate(sex = "male")
Labonte_Female_magma <- Labonte_magma %>% select(brain_region, symbol, logFC = Female.logFC, pvalue = Female.pvalue) %>% mutate(sex = "female")

#combine 6 (6 regions * 1 sexes) 
Labonte_Female_magma <- LabonteMeta(Labonte_Female_magma,6,6)

#combine 6 (6 regions * 1 sexes) 
Labonte_Male_magma <- LabonteMeta(Labonte_Male_magma,6,6)

Labonte_summary_magma %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "fullLabonteTable_magma.csv")) 
Labonte_Female_magma %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "FemaleLabonteTable_magma.csv"))
Labonte_Male_magma %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "MaleLabonteTable_magma.csv"))


fullLabonte_meta <- read_csv(here("ProcessedData", "LabonteEtAl", "fullLabonteTable_magma.csv"))
femaleLabonte_meta <- read_csv(here("ProcessedData", "LabonteEtAl", "FemaleLabonteTable_magma.csv"))
maleLabonte_meta <- read_csv(here("ProcessedData", "LabonteEtAl", "MaleLabonteTable_magma.csv"))

num_genes_meta <- length(fullLabonte_meta$gene_symbol)
num_genes_female <- length(femaleLabonte_meta$gene_symbol)
num_genes_males <- length(maleLabonte_meta$gene_symbol)

fullLabonte <- getRank(fullLabonte_meta, num_genes_meta)
merged_full <- left_join(fullLabonte_meta, fullLabonte %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_full %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "fullLabonteTable_magma.csv"))

femaleLabonte <- getRank(femaleLabonte_meta, num_genes_female)
merged_female <- left_join(femaleLabonte_meta, femaleLabonte %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_female %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "FemaleLabonteTable_magma.csv"))

maleLabonte <- getRank(maleLabonte_meta, num_genes_males)
merged_male <- left_join(maleLabonte_meta, maleLabonte %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_male %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "MaleLabonteTable_magma.csv"))


#--------------------- Non-magma list ---------------------------------
# #get the meta-p value from the analysis
# fullLabonte_meta <- read_csv(here("ProcessedData", "LabonteEtAl", "fullLabonteTable.csv"))
# femaleLabonte_meta <- read_csv(here("ProcessedData", "LabonteEtAl", "FemaleLabonteTable.csv"))
# maleLabonte_meta <- read_csv(here("ProcessedData", "LabonteEtAl", "MaleLabonteTable.csv"))
# 
# num_genes_meta <- length(fullLabonte_meta$gene_symbol)
# num_genes_female <- length(femaleLabonte_meta$gene_symbol)
# num_genes_males <- length(maleLabonte_meta$gene_symbol)
# 
# fullLabonte <- getRank(fullLabonte_meta, num_genes_meta)
# merged_full <- left_join(fullLabonte_meta, fullLabonte %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
# merged_full %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "fullLabonteTable.csv"))
# 
# femaleLabonte <- getRank(femaleLabonte_meta, num_genes_female)
# merged_female <- left_join(femaleLabonte_meta, femaleLabonte %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
# merged_female %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "FemaleLabonteTable.csv"))
# 
# maleLabonte <- getRank(maleLabonte_meta, num_genes_males)
# merged_male <- left_join(maleLabonte_meta, maleLabonte %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
# merged_male %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "MaleLabonteTable.csv"))
# 

# 
# 
# 
# 
# 
# 
# 
# 
# 







