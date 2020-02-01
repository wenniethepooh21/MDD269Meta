library(readr)

library(tidyr)
library(magrittr)
library(dplyr)
library(GEOquery)

setwd('../../../school/thesis/')
library(here)
#based on signature.R by the biojupies team
#https://github.com/MaayanLab/biojupies-plugins/blob/1024a6ed702ad8b0958d4ccdd2afe89cbe493a51/library/core_scripts/signature/signature.R
#from https://amp.pharm.mssm.edu/biojupies/notebook/tlNmgbMxY (just load biojupies with GSE80655)

#Method for getting sequence counts
#"Raw RNA-seq data for GEO dataset GSE80655 was downloaded from the SRA database (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80655) and quantified to gene-level counts using the ARCHS4 pipeline (Lachmann et al., 2017). Gene counts were downloaded from the ARCHS4 gene expression matrix v6. For more information about ARCHS4, as well as free access to the quantified gene expression matrix, visit the project home page at the following URL: http://amp.pharm.mssm.edu/archs4/download.html."
#(Biojupies)

source(here("R/transcriptomic_meta/RamakerFullMetaAnalysis.R"))

metadata <- read_tsv(here("Raw_Data/RamakerEtAl/GSE80655/GSE80655-metadata.txt.txt"))
regions <- unique(metadata$`brain region`)
full_results <- tibble()

#counts from Biojupies
rawcount_dataframe <- read.csv(here("data", "RamakerEtAl", "GSE80655 from biojupies", "GSE80655-expression.txt.txt"),sep="\t", row.names = 1)

gds <- getGEO("GSE80655")
metadata_for_pH <- phenoData(gds$GSE80655_series_matrix.txt.gz)
metadata_for_pH <- as_tibble(as(metadata_for_pH, "data.frame"))
metadata_for_pH %<>% select(Sample_geo_accession = geo_accession, brain_ph = 'brain ph:ch1')
metadata_for_pH %<>% mutate(brain_ph = as.numeric(brain_ph))
print(paste("NA pH values:" , nrow(metadata_for_pH %>% filter(is.na(brain_ph)))))
#replacing 3 missing values with mean brain ph
metadata_for_pH %<>% mutate(brain_ph = if_else(is.na(brain_ph), mean(brain_ph, na.rm=T), brain_ph))

#read count data
read_counts <- as.data.frame(colSums(rawcount_dataframe))
read_counts$Sample_geo_accession = rownames(read_counts)
read_counts <- as_tibble(read_counts) %>% select(Sample_geo_accession , read_counts = `colSums(rawcount_dataframe)`)

metadata %<>% rename(clinical_diagnosis = `clinical diagnosis`)
metadata %<>% mutate(clinical_diagnosis = gsub(" ", "_", clinical_diagnosis))
metadata <- inner_join(metadata, metadata_for_pH)
metadata <- inner_join(metadata, read_counts)

female_metadata <- metadata[which(metadata$gender == "F"),]
male_metadata <- metadata[which(metadata$gender == "M"),]

#full Ramaker data results
summary_results <- RamakerMeta(metadata, read_counts, rawcount_dataframe, regions, full_results)
write_csv(summary_results, path = here("ProcessedData", "RamakerEtAl", "CompleteRamakerTable.csv"))
summary_results %<>% RamakerAnalysis(regions)

#full female Ramaker data results
female_summary_results <- RamakerMeta(female_metadata, read_counts, rawcount_dataframe, regions, full_results)
write_csv(female_summary_results, path = here("ProcessedData", "RamakerEtAl", "CompleteFemaleRamakerTable.csv"))
female_summary_results %<>% RamakerAnalysis(regions)
female_summary_results %<>% rename(AnCg_nAcc_DLPFC_Female_directions = AnCg_nAcc_DLPFC_directions)

#full male Ramaker data results
male_summary_results <- RamakerMeta(male_metadata, read_counts, rawcount_dataframe, regions, full_results)
write_csv(male_summary_results, path = here("ProcessedData", "RamakerEtAl", "CompleteMaleRamakerTable.csv"))
male_summary_results %<>% RamakerAnalysis(regions)
male_summary_results %<>% rename(AnCg_nAcc_DLPFC_Male_directions = AnCg_nAcc_DLPFC_directions)


#merge all gender directions into summary_results table for better visualization
#summary <- summary_results %>% select(gene_symbol, AnCg_nAcc_DLPFC_directions)

summary <- left_join(female_summary_results %>% select(gene_symbol,AnCg_nAcc_DLPFC_Female_directions), male_summary_results %>% select(gene_symbol,AnCg_nAcc_DLPFC_Male_directions ))
summary %<>% unite(AnCg.F_nAcc.F_DLPFC.F_AnCg.M_nAcc.M_DLPFC.M, AnCg_nAcc_DLPFC_Female_directions, AnCg_nAcc_DLPFC_Male_directions, sep = "")
summary %<>% left_join(summary_results) %>% select(-AnCg_nAcc_DLPFC_directions) %>% distinct()

#Save data in .csv files
summary %>% write_csv(path = here("ProcessedData", "RamakerEtAl", "fullRamakerTable.csv"))
female_summary_results %>% write_csv(path = here("ProcessedData", "RamakerEtAl", "FemaleRamakerTable.csv"))
male_summary_results %>% write_csv(path = here("ProcessedData", "RamakerEtAl", "MaleRamakerTable.csv"))


#---------------------- Genome Ranking MAGMA --------------------------#
magma <- read_csv(here("data", "HowardEtAl", "FullMagmaGenes.csv"))
full_Ramaker <- read_csv( here("ProcessedData", "RamakerEtAl","CompleteRamakerTable.csv"))
full_Ramaker %<>% mutate(gene_symbol = gsub("C([X0-9]+)ORF([0-9]+)", "C\\1orf\\2", gene_symbol)) #change from upper case to lower case 
full_Ramaker %<>% right_join(magma %>% select(Ramaker_genes) %>% distinct(), by = c('gene_symbol' = 'Ramaker_genes')) %>% na.omit() %>% write_csv(here("ProcessedData","RamakerEtAl", "CompleteRamakerTableMagma.csv"))
full_Ramaker %<>% RamakerAnalysis(regions)


female_Ramaker_magma <- read_csv(here("ProcessedData", "RamakerEtAl", "CompleteFemaleRamakerTable.csv"))
female_Ramaker_magma %<>% mutate(gene_symbol = gsub("C([X0-9]+)ORF([0-9]+)", "C\\1orf\\2", gene_symbol))
female_Ramaker_magma  %<>% right_join(magma %>% select(Ramaker_genes) %>% distinct(), by = c('gene_symbol' = 'Ramaker_genes')) %>% na.omit() %>% write_csv(here("ProcessedData", "RamakerEtAl","CompleteFemaleRamakerTableMagma.csv"))
female_Ramaker_magma %<>% RamakerAnalysis(regions)
female_Ramaker_magma %<>% rename(AnCg_nAcc_DLPFC_Female_directions = AnCg_nAcc_DLPFC_directions)

male_Ramaker_magma <- read_csv(here("ProcessedData", "RamakerEtAl", "CompleteMaleRamakerTable.csv"))
male_Ramaker_magma %<>% mutate(gene_symbol = gsub("C([X0-9]+)ORF([0-9]+)", "C\\1orf\\2", gene_symbol))
male_Ramaker_magma  %<>% right_join(magma %>% select(Ramaker_genes) %>% distinct(), by = c( 'gene_symbol' = 'Ramaker_genes')) %>% na.omit() %>% write_csv(here("ProcessedData","RamakerEtAl", "CompleteMaleRamakerTableMagma.csv"))
male_Ramaker_magma %<>% RamakerAnalysis(regions)
male_Ramaker_magma %<>% rename(AnCg_nAcc_DLPFC_Male_directions = AnCg_nAcc_DLPFC_directions)

summary_magma <- left_join(female_Ramaker_magma %>% select(gene_symbol,AnCg_nAcc_DLPFC_Female_directions), male_Ramaker_magma %>% select(gene_symbol,AnCg_nAcc_DLPFC_Male_directions ))
summary_magma %<>% unite(AnCg.F_nAcc.F_DLPFC.F_AnCg.M_nAcc.M_DLPFC.M, AnCg_nAcc_DLPFC_Female_directions, AnCg_nAcc_DLPFC_Male_directions, sep = "")
summary_magma %<>% left_join(full_Ramaker) %>% select(-AnCg_nAcc_DLPFC_directions) %>% distinct()

#Save data in .csv files
summary_magma %>% write_csv(path = here("ProcessedData", "RamakerEtAl", "fullRamakerTable_magma.csv"))
female_Ramaker_magma %>% write_csv(path = here("ProcessedData", "RamakerEtAl", "FemaleRamakerTable_magma.csv"))
male_Ramaker_magma %>% write_csv(path = here("ProcessedData", "RamakerEtAl", "MaleRamakerTable_magma.csv"))

#get the meta-p value from the analysis
fullRamaker_meta <- read_csv(here("ProcessedData", "RamakerEtAl", "fullRamakerTable_magma.csv"))
femaleRamaker_meta <- read_csv(here("ProcessedData", "RamakerEtAl", "FemaleRamakerTable_magma.csv"))
maleRamaker_meta <- read_csv(here("ProcessedData", "RamakerEtAl", "MaleRamakerTable_magma.csv"))


num_genes_meta <- length(fullRamaker_meta$gene_symbol)
num_genes_female <- length(femaleRamaker_meta$gene_symbol)
num_genes_males <- length(maleRamaker_meta$gene_symbol)


fullRamaker <- getRank(fullRamaker_meta, num_genes_meta)
merged_full <- left_join(fullRamaker_meta, fullRamaker %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_full %>% write_csv(path = here("ProcessedData", "RamakerEtAl", "fullRamakerTable_magma.csv"))

femaleRamaker <- getRank(femaleRamaker_meta, num_genes_female)
merged_female <- left_join(femaleRamaker_meta, femaleRamaker %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_female %>% write_csv(path = here("ProcessedData", "RamakerEtAl", "FemaleRamakerTable_magma.csv"))

maleRamaker <- getRank(maleRamaker_meta, num_genes_males)
merged_male <- left_join(maleRamaker_meta, maleRamaker %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
merged_male %>% write_csv(path = here("ProcessedData", "RamakerEtAl", "MaleRamakerTable_magma.csv"))



#---------------------- Genome Ranking NON MAGMA --------------------------#

# #get the meta-p value from the analysis
# fullRamaker_meta <- read_csv(here("ProcessedData", "RamakerEtAl", "fullRamakerTable.csv"))
# femaleRamaker_meta <- read_csv(here("ProcessedData", "RamakerEtAl", "FemaleRamakerTable.csv"))
# maleRamaker_meta <- read_csv(here("ProcessedData", "RamakerEtAl", "MaleRamakerTable.csv"))
# 
# 
# num_genes_meta <- length(fullRamaker_meta$gene_symbol)
# num_genes_female <- length(femaleRamaker_meta$gene_symbol)
# num_genes_males <- length(maleRamaker_meta$gene_symbol)
# 
# 
# fullRamaker <- getRank(fullRamaker_meta, num_genes_meta)
# merged_full <- left_join(fullRamaker_meta, fullRamaker %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
# merged_full %>% write_csv(path = here("ProcessedData", "RamakerEtAl", "fullRamakerTable.csv"))
# 
# femaleRamaker <- getRank(femaleRamaker_meta, num_genes_female)
# merged_female <- left_join(femaleRamaker_meta, femaleRamaker %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
# merged_female %>% write_csv(path = here("ProcessedData", "RamakerEtAl", "FemaleRamakerTable.csv"))
# 
# maleRamaker <- getRank(maleRamaker_meta, num_genes_males)
# merged_male <- left_join(maleRamaker_meta, maleRamaker %>% dplyr::select(gene_symbol, meta_Up, meta_Down, genome_percentile_rank), by = c("gene_symbol" = "gene_symbol"))
# merged_male %>% write_csv(path = here("ProcessedData", "RamakerEtAl", "MaleRamakerTable.csv"))

