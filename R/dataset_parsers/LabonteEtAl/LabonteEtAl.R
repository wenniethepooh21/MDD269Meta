library(magrittr)
library(tidy)
library(here)
#library(biomaRt)
library(dplyr)
library(readxl)


#Read in raw male expression data
Labonte_Male <- read_xlsx(here("Raw_Data/LabonteEtAl/Male DEG all.xlsx"))
Labonte_Male %<>% rename(gene_symbol  = "Gene name", brain_region = "Brain Region", Male.pvalue = "p-value", Male.logFC = logFC, Male.AveExpr = AveExpr)

#Read in raw female expression data
Labonte_Female <- read_xlsx(here("Raw_Data/LabonteEtAl/Female DEG all.xlsx"))
Labonte_Female %<>% rename(brain_region = `Brain region`, gene_symbol = `Gene Name`, Female.pvalue = `p-value`, Female.logFC = logFC, Female.AveExpr = AveExpr)

Labonte <- left_join(Labonte_Male %>% select(-Biotype, -Description), Labonte_Female %>% select(-Biotype, -Description), by = c("brain_region", "ENSG", "gene_symbol"))

########################################### update excel converted gene symbols from dates to gene symbols using ENSG number################
#convert ENSG to gene gene_symbol to remove dates from gene_symbols column
#ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

#This line was previously run to generate associated gene name, output was saved to save time
#ENSG_Genes <- getBM(attributes = c('ensembl_gene_id','hgnc_gene_symbol'), filters = 'ensembl_gene_id', values = unique(Labonte$ENSG), mart = ensembl) %>% as_tibble() 
#ENSG_Genes %>% write_csv(here("Processed_Data/LabonteEtAl/getBMGenes.csv"))

ENSG_Genes <- read_csv(here("Processed_Data/LabonteEtAl/getBMGenes.csv")) # read in gene mapping file 
Labonte_genes <- Labonte %>% select(ENSG, gene_symbol)
Labonte_genes %<>% left_join(ENSG_Genes, by = c('ENSG' = 'ensembl_gene_id'))

#merge gene gene_symbols into one column 
Labonte_genes %<>% rowwise() %>% mutate(hgnc_symbol = if_else(is.na(hgnc_symbol), gene_symbol, hgnc_symbol))
Labonte_genes %<>% select(-gene_symbol)	
Labonte_genes %<>% distinct()
#Find if all ENSG's are distinct
Labonte_genes %>% group_by(ENSG) %>% filter(n() > 1)
#keep the hgnc symbol == "LINC00856", remove for "LINC00595" for ENSG00000230417 == LINC00856
Labonte_genes %<>% mutate(hgnc_symbol = if_else(ENSG == "ENSG00000230417" & hgnc_symbol == "LINC00595", "LINC00856", hgnc_symbol)) %>% distinct() 
Labonte %<>% left_join(Labonte_genes)
Labonte %<>% select(-gene_symbol)
Labonte %<>% rename(gene_symbol = "hgnc_symbol")
############################################################################################
#re-order columns
Labonte <- Labonte[,c(1:2,9,3:8)]
Labonte %<>% mutate(brain_region = ifelse(brain_region == "Anterior Insula", "Anterior_Insula",brain_region))

#filter transcriptomic dataset for MAGMA genes tested by Howard, et al. 
magma_table <- read_csv(here("Raw_Data/HowardEtAl/FullMagmaGenes.csv")) %>% select(Labonte_genes) %>% distinct() %>% na.omit()
magma_labonte <- Labonte %>% right_join(magma_table , by = c('gene_symbol' = 'Labonte_genes'))
#Write out formatted expression data 
magma_labonte %>% write_csv(path = here("Processed_Data/LabonteEtAl/CompleteLabonteTableMagma.csv"))

#Expand table row wise to hold male and female data 
Labonte_long <- bind_rows(magma_labonte %>% select(brain_region, gene_symbol, logFC = Male.logFC, pvalue = Male.pvalue) %>% mutate(sex = "male"),
                          magma_labonte %>% select(brain_region, gene_symbol, logFC = Female.logFC, pvalue = Female.pvalue) %>% mutate(sex = "female"))


#Perform meta-analyses on Labonte dataset 
source(here("R/transcriptomic_meta/Labonte_Meta_Analysis.R"))
#Perform meta-analysis on all brain regions across both sexes
#number of p-values we need to filter = 12 (6 regions * 2 sexes) 
#number of brain regions = 6
Labonte_summary_results <- Labonte_long %>% LabonteMetaAnalysis(12,6)
Labonte_summary_results %>% write_csv(here("Processed_Data/LabonteEtAl/fullLabonteTableMagma.csv"))

#Perform the sex-specific meta-analyses
#Extract Female data 
Labonte_Female_Magma <- magma_labonte %>% select(brain_region, gene_symbol, logFC = Female.logFC, pvalue = Female.pvalue) %>% mutate(sex = "female")
#Perform meta-analysis on all brain regions in female data 
#number of p-values we need to filter = 6 (6 regions * 1 sex) 
#number of brain regions = 6
Labonte_Female_results <- Labonte_Female_Magma %>% LabonteMetaAnalysis(6,6)
Labonte_Female_results %>% write_csv(here("Processed_Data/LabonteEtAl/FemaleLabonteTableMagma.csv"))

#Extract Male data
Labonte_Male_Magma <- magma_labonte %>% select(brain_region, gene_symbol, logFC = Male.logFC, pvalue = Male.pvalue) %>% mutate(sex = "male")
#Perform meta-analysis on all brain regions in female data 
#number of p-values we need to filter = 6 (6 regions * 1 sex) 
#number of brain regions = 6
Labonte_Male_results <- Labonte_Male_Magma %>% LabonteMetaAnalysis(6,6)
Labonte_Male_results %>% write_csv(path = here("Processed_Data/LabonteEtAl/MaleLabonteTableMagma.csv"))


############
# CORTICAL ANALYSIS
############
#Pick out only the cortical regions that Labonte studied
Labonte_Cortical <- Labonte_long %>% filter(brain_region != "Nac") %>% filter( brain_region != "Subic")
#Perform meta-analysis on all brain regions across both sexes
#number of p-values we need to filter = 8 (4 regions * 2 sexes) 
#number of brain regions = 4
cortical_Labonte_summary_results <- Labonte_Cortical %>% LabonteMeta(8,4)
cortical_Labonte_summary_results %>% write_csv(path = here("Processed_Data/LabonteEtAl/CorticalLabonteTableMagma.csv"))


##################
# SEX-INTERACTION
##################

#Flip male logFC value
Labonte_long_flip <- Labonte_long %>% mutate(logFC = if_else(sex == "male",logFC *-1, logFC))


Labonte_long <- bind_rows(Labonte %>% dplyr::select(brain_region, symbol, logFC = Male.logFC, pvalue = Male.pvalue) %>% mutate(sex = "male"),
                          Labonte %>% dplyr::select(brain_region, symbol, logFC = Female.logFC, pvalue = Female.pvalue) %>% mutate(sex = "female"))

#combine 12 (6 regions * 2 sexes) + 12 pvalues for each gene
Labonte_summary_results <- LabonteMeta(Labonte_long,12,6)
#change directions back to original

#unflipped directions 
labonte_unfliped <- read_csv(here("ProcessedData", "LabonteEtAl", "fullLabonteTable.csv"))
labonte_unfliped %<>% select(1:2)
labonte_unfliped %<>% right_join(Labonte_summary_results %>% select(-2))

#needs writing out
labonte_unfliped %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "fullLabonteTableMagma_flipped.csv"))

# #------- Genome ranking (MAGMA list)
magma <- read_csv(here("data", "HowardEtAl", "FullMagmaGenes.csv"))
full_Labonte <- read_csv( here("ProcessedData", "LabonteEtAl", "CompleteLabonteTable.csv")) %>% na.omit()
full_Labonte %<>% mutate(gene_symbol = gsub("C([X0-9]+)ORF([0-9]+)", "C\\1orf\\2", gene_symbol))

Labonte_magma <- left_join(magma %>% select(Labonte_genes) %>% distinct(), full_Labonte, by = c('Labonte_genes' = 'gene_symbol')) %>% na.omit() #remove the one NA value
Labonte_magma %<>% rename(gene_symbol = Labonte_genes) %>% write_csv(here("ProcessedData", "LabonteEtAl", "CompleteLabonteTableMagma.csv"))

Labonte_long_magma <- bind_rows(Labonte_magma %>% dplyr::select(brain_region, gene_symbol, logFC = Male.logFC, pvalue = Male.pvalue) %>% mutate(sex = "male"),
                          Labonte_magma %>% dplyr::select(brain_region, ENSG, gene_symbol, logFC = Female.logFC, pvalue = Female.pvalue) %>% mutate(sex = "female"))

#combine 12 (6 regions * 2 sexes) + 12 pvalues for each gene
Labonte_summary_magma <- LabonteMeta(Labonte_long_magma,12,6)

#get the sex data separately
Labonte_Male_magma <- Labonte_magma %>% select(brain_region, gene_symbol, logFC = Male.logFC, pvalue = Male.pvalue) %>% mutate(sex = "male")
Labonte_Female_magma <- Labonte_magma %>% select(brain_region, gene_symbol, logFC = Female.logFC, pvalue = Female.pvalue) %>% mutate(sex = "female")

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

num_genes_meta <- length(fullLabonte_meta$gene_gene_symbol)
num_genes_female <- length(femaleLabonte_meta$gene_gene_symbol)
num_genes_males <- length(maleLabonte_meta$gene_gene_symbol)

fullLabonte <- getRank(fullLabonte_meta, num_genes_meta)
merged_full <- left_join(fullLabonte_meta, fullLabonte %>% dplyr::select(gene_gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_gene_symbol" = "gene_gene_symbol"))
merged_full %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "fullLabonteTable_magma.csv"))

femaleLabonte <- getRank(femaleLabonte_meta, num_genes_female)
merged_female <- left_join(femaleLabonte_meta, femaleLabonte %>% dplyr::select(gene_gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_gene_symbol" = "gene_gene_symbol"))
merged_female %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "FemaleLabonteTable_magma.csv"))

maleLabonte <- getRank(maleLabonte_meta, num_genes_males)
merged_male %>% write_csv(path = here("ProcessedData", "LabonteEtAl", "MaleLabonteTable_magma.csv"))


###############################################################
###### GENOME PERCENTILE RANKING FOR ALL ABOVE ANALYSES ######
##############################################################
#load script that holds the genome percentile ranking function used by all transcriptomic studies
source(here("R/transcriptomic_meta/Percentile_Rank_Analysis.R"))







