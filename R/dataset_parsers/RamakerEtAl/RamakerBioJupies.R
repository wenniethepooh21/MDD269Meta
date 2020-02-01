library(readr)
library(tidyr)
library(GEOquery)
library(magrittr)
library(dplyr)


setwd('../../../school/thesis/')
library(here)
#based on signature.R by the biojupies team
#https://github.com/MaayanLab/biojupies-plugins/blob/1024a6ed702ad8b0958d4ccdd2afe89cbe493a51/library/core_scripts/signature/signature.R
#from https://amp.pharm.mssm.edu/biojupies/notebook/tlNmgbMxY (just load biojupies with GSE80655)

#Method for getting sequence counts
#"Raw RNA-seq data for GEO dataset GSE80655 was downloaded from the SRA database (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80655) and quantified to gene-level counts using the ARCHS4 pipeline (Lachmann et al., 2017). Gene counts were downloaded from the ARCHS4 gene expression matrix v6. For more information about ARCHS4, as well as free access to the quantified gene expression matrix, visit the project home page at the following URL: http://amp.pharm.mssm.edu/archs4/download.html."
#(Biojupies)

gds <- getGEO("GSE80655")
metadata_for_pH <- phenoData(gds$GSE80655_series_matrix.txt.gz) 
metadata_for_pH <- as_tibble(as(metadata_for_pH, "data.frame"))
metadata_for_pH %<>% select(Sample_geo_accession = geo_accession, brain_ph = 'brain ph:ch1')
#count the number of missing brain pH values from samples
print(paste("NA pH values:" , nrow(metadata_for_pH %>% filter(is.na(brain_ph)))))
# convert column into numeric column, will get warning from NA's
metadata_for_pH %<>% mutate(brain_ph = as.numeric(brain_ph))
#replacing 3 missing values with mean brain ph
metadata_for_pH %<>% mutate(brain_ph = if_else(is.na(brain_ph), mean(brain_ph, na.rm=T), brain_ph))


#read rawcounts data downloaded from their GSE
rawcount_dataframe <- read.csv(here("Raw_Data/RamakerEtAl/GSE80655/GSE80655-expression.txt.txt"),sep="\t", row.names = 1)
read_counts <- as.data.frame(colSums(rawcount_dataframe))
#create a new column that copied the GSM id's to convert to tibble 
read_counts$Sample_geo_accession = rownames(read_counts)
#convert to tibble (removes the first column)
read_counts <- as_tibble(read_counts) %>% select(Sample_geo_accession , read_counts = `colSums(rawcount_dataframe)`)

#read in the meta data downloaded from their GSE
metadata <- read_tsv(here("Raw_Data/RamakerEtAl/GSE80655/GSE80655-metadata.txt.txt"))
metadata %<>% rename(clinical_diagnosis = `clinical diagnosis`)
metadata %<>% mutate(clinical_diagnosis = gsub(" ", "_", clinical_diagnosis))
#merge data together
fullmetadata <- inner_join(metadata, metadata_for_pH) %>% inner_join(read_counts)
#separate the data for sex-specific analysis 
female_metadata <- fullmetadata %>% filter(gender == "F")
male_metadata <- fullmetadata %>% filter(gender == "M") 

source(here("R/transcriptomic_meta/Ramaker_Meta_Analysis.R"))
#Read in the associated list of MAGMA genes that howard tested (17,842)
magma_table <- read_csv(here("Raw_Data/HowardEtAl/FullMagmaGenes.csv"))
#get the unique brain regions
regions <- unique(metadata$`brain region`)
#create empty tibble for data to be populated 
full_results <- tibble()

#Perform Ramaker meta-analysis functions in RamakerMetaAnalysis.R
R_summary_results <- RamakerDEModel(fullmetadata, read_counts, rawcount_dataframe, regions, full_results)
R_summary_results %<>% mutate(gene_symbol = gsub("C([X0-9]+)ORF([0-9]+)", "C\\1orf\\2", gene_symbol)) #change from upper case to lower case for open reading frame genes
#Filter Ramaker genome for the MAGMA genes tested by Howard
R_summary_results %<>% right_join(magma_table %>% select(Ramaker_genes) %>% distinct(), by = c('gene_symbol' = 'Ramaker_genes')) %>% na.omit()
R_summary_results %>% write_csv(here("Processed_Data/RamakerEtAl/CompleteRamakerTableMagma.csv")) #write out for easier access in other analyses (cortical)
#perform meta-analysis on full (female and male all brain regions) data
R_summary_results %<>% RamakerMetaAnalysis(regions)

#create empty tibble for data to be populated 
female_results <- tibble()
#full female Ramaker data results
R_female_summary_results <- RamakerDEModel(female_metadata, read_counts, rawcount_dataframe, regions, female_results)
R_female_summary_results %<>% mutate(gene_symbol = gsub("C([X0-9]+)ORF([0-9]+)", "C\\1orf\\2", gene_symbol)) #change from upper case to lower case for open reading frame genes
#Filter Ramaker genome for the MAGMA genes tested by Howard
R_female_summary_results %<>% right_join(magma_table %>% select(Ramaker_genes) %>% distinct(), by = c('gene_symbol' = 'Ramaker_genes')) %>% na.omit()
R_female_summary_results %>% write_csv(here("Processed_Data/RamakerEtAl/CompleteFemaleRamakerTableMagma.csv"))#write out for easier access in other analyses
#perform meta-analysis on female data
R_female_summary_results %<>% RamakerMetaAnalysis(regions)
R_female_summary_results %<>% rename(AnCg_nAcc_DLPFC_Female_directions = AnCg_nAcc_DLPFC_directions)
#Save full female meta-analysis results
R_female_summary_results %>% write_csv(here("Processed_Data/RamakerEtAl/FemaleRamakerTableMagma.csv"))

#create empty tibble for data to be populated 
male_results <- tibble()
#full male Ramaker data results
R_male_summary_results <- RamakerDEModel(male_metadata, read_counts, rawcount_dataframe, regions, male_results)
R_male_summary_results %<>% mutate(gene_symbol = gsub("C([X0-9]+)ORF([0-9]+)", "C\\1orf\\2", gene_symbol)) #change from upper case to lower case for open reading frame genes
R_male_summary_results %<>% right_join(magma_table %>% select(Ramaker_genes) %>% distinct(), by = c('gene_symbol' = 'Ramaker_genes')) %>% na.omit()
R_male_summary_results %>% write_csv(here("Processed_Data/RamakerEtAl/CompleteMaleRamakerTableMagma.csv"))#write out for easier access in other analyses
#perform meta-analysis on male data
R_male_summary_results %<>% RamakerMetaAnalysis(regions)
R_male_summary_results %<>% rename(AnCg_nAcc_DLPFC_Male_directions = AnCg_nAcc_DLPFC_directions)
#Save full male meta-analysis results
R_male_summary_results %>% write_csv(here("Processed_Data/RamakerEtAl/MaleRamakerTableMagma.csv"))

#merge all gender directions into summary_results table full meta analysis visualization
Ramaker_summary <- left_join(R_female_summary_results %>% select(gene_symbol,AnCg_nAcc_DLPFC_Female_directions), R_male_summary_results %>% select(gene_symbol,AnCg_nAcc_DLPFC_Male_directions ))
Ramaker_summary %<>% unite(AnCg.F_nAcc.F_DLPFC.F_AnCg.M_nAcc.M_DLPFC.M, AnCg_nAcc_DLPFC_Female_directions, AnCg_nAcc_DLPFC_Male_directions, sep = "")
Ramaker_summary %<>% left_join(R_summary_results) %>% select(-AnCg_nAcc_DLPFC_directions) %>% distinct()
#Save full meta-analysis results
Ramaker_summary %>% write_csv(here("Processed_Data/RamakerEtAl/fullRamakerTableMagma.csv")) #used for combining across transcriptomic studies 

########################################################
###### CORTICAL ANALYSIS (FULL, FEMALE AND MALE) ######
########################################################
#extract cortical data and re-run meta-analysis -- don't have to re-run model creation
Ramaker_cortical<- read_csv(here("Processed_Data/RamakerEtAl/CompleteRamakerTableMagma.csv"))
#Extract the cortical region data & run analysis
Ramaker_cortical %<>% filter(target_region != "nAcc")
#run the meta-analysis on only the cortical regions sampled
Ramaker_cortical %<>% RamakerMetaAnalysis(regions)

#Extract the cortical region data & run analysis on female data
female_ramaker_cortical <- read_csv(here("ProcessedData", "RamakerEtAl", "CompleteFemaleRamakerTableMagma.csv"))
female_ramaker_cortical %<>% filter(target_region != "nAcc")
#run the meta-analysis on only the cortical regions sampled in female data
female_ramaker_cortical %<>% RamakerMetaAnalysis(regions)
female_ramaker_cortical %<>% rename(AnCg_DLPFC_Female_directions = AnCg_DLPFC_directions)
female_ramaker_cortical %>% write_csv(here("ProcessedData", "RamakerEtAl", "fullCorticalFemaleRamakerTableMagma.csv"))

#Extract the cortical region data & run analysis on male data
male_ramaker_cortical<- read_csv(here("ProcessedData", "RamakerEtAl", "CompleteMaleRamakerTableMagma.csv"))
male_ramaker_cortical %<>% filter(target_region != "nAcc")
#run the meta-analysis on only the cortical regions sampled in male data
male_ramaker_cortical %<>% RamakerMetaAnalysis(regions)
male_ramaker_cortical %<>% rename(AnCg_DLPFC_Male_directions = AnCg_DLPFC_directions)
male_ramaker_cortical %>% write_csv(here("ProcessedData", "RamakerEtAl", "fullCorticalMaleRamakerTableMagma.csv"))

#merge all directions from male and female data to visualize cortical directions across sexes
cortical_summary <- left_join(female_ramaker_cortical %>% select(gene_symbol,AnCg_DLPFC_Female_directions), male_ramaker_cortical %>% select(gene_symbol,AnCg_DLPFC_Male_directions ))
cortical_summary %<>% unite("AnCg.F_DLPFC.F_AnCg.M_DLPFC.M", AnCg_DLPFC_Female_directions, AnCg_DLPFC_Male_directions, sep = "")
cortical_summary %<>% left_join(Ramaker_cortical) %>% select(-AnCg_DLPFC_directions) %>% distinct()
#Save full cortical meta-analysis results 
cortical_summary %>% write_csv(here("ProcessedData", "RamakerEtAl", "fullCorticalRamakerTableMagma.csv"))

###########################################################
###### SEX-INTERACTION ANALYSIS (FULL AND CORTICAL) ######
###########################################################
#### Full meta-analysis

#Read in the post-model female data 
full_female_summary_results <- read_csv(here("Processed_Data/RamakerEtAl/CompleteFemaleRamakerTableMagma.csv"))
#add sex for unique identifier 
full_female_summary_results %<>% mutate(sex = "female")

#Read in the post-model female data 
full_male_summary_results <- read_csv(here("ProcessedData", "RamakerEtAl", "CompleteMaleRamakerTableMagma.csv"))
full_male_summary_results_flip <- male_summary_results %>% mutate(t = t*-1)
#add sex for unique identifier
male_summary_results_flip %<>% mutate(sex = "male")

#merge male summary results with female summary results into one table to perform calculations
full_flipped <- rbind(male_summary_results_flip, female_summary_results)
full_flipped %<>% write_csv(here("ProcessedData", "RamakerEtAl", "CompleteRamakerTableMagma_flipped.csv"))

#Calculations from female data and male flipped data combined
full_flipped %<>% RamakerMetaAnalysis(regions)
#flips the directions of the male gene expressions back to original
full_flipped %<>% flipDirections("male")

#Calculate the new male data flipped
male_summary_results_flip %<>% RamakerMetaAnalysis(regions)
male_summary_results_flip %<>% flipDirections("male") #flip the male directions 
male_summary_results_flip %<>% rename(AnCg_nAcc_DLPFC_Male_directions = AnCg_nAcc_DLPFC_directions)

#Read in Female meta-analysis results to get the directions of gene expression across all brain regions
full_female_results <- read_csv(here("ProcessedData", "RamakerEtAl", "FemaleRamakerTableMagma.csv"))
#combine male directions and female directions 
full_flipped_summary <- left_join(full_female_results %>% select(gene_symbol, AnCg_nAcc_DLPFC_Female_directions), male_summary_results_flip %>% select(gene_symbol, AnCg_nAcc_DLPFC_Male_directions))
#concatenate the two directions columns
full_flipped_summary %<>% unite("AnCg.F_nAcc.F_DLPFC.F_AnCg.M_nAcc.M_DLPFC.M",AnCg_nAcc_DLPFC_Female_directions, AnCg_nAcc_DLPFC_Male_directions, sep = "")
#Combine directions with meta values calculated 
full_flipped %<>% select(-sex, -AnCg_nAcc_DLPFC_directions)
#join calculations for each gene
full_flipped_summary %<>% left_join(full_flipped) %>% distinct()

full_flipped_summary%>%write_csv(here("ProcessedData", "RamakerEtAl", "fullRamakerTableMagma_flipped.csv"))


















#---------------------- Genome Ranking MAGMA --------------------------#
magma <- read_csv(here("data", "HowardEtAl", "FullMagmaGenes.csv"))
full_Ramaker <- read_csv( here("ProcessedData", "RamakerEtAl","CompleteRamakerTable.csv"))

full_Ramaker %<>% right_join(magma %>% select(Ramaker_genes) %>% distinct(), by = c('gene_symbol' = 'Ramaker_genes')) %>% na.omit() %>% write_csv(here("ProcessedData","RamakerEtAl", "CompleteRamakerTableMagma.csv"))
full_Ramaker %<>% RamakerMetaAnalysis(regions)


female_Ramaker_magma <- read_csv(here("ProcessedData", "RamakerEtAl", "CompleteFemaleRamakerTable.csv"))
female_Ramaker_magma %<>% mutate(gene_symbol = gsub("C([X0-9]+)ORF([0-9]+)", "C\\1orf\\2", gene_symbol))
female_Ramaker_magma  %<>% right_join(magma %>% select(Ramaker_genes) %>% distinct(), by = c('gene_symbol' = 'Ramaker_genes')) %>% na.omit() %>% write_csv(here("ProcessedData", "RamakerEtAl","CompleteFemaleRamakerTableMagma.csv"))
female_Ramaker_magma %<>% RamakerMetaAnalysis(regions)
female_Ramaker_magma %<>% rename(AnCg_nAcc_DLPFC_Female_directions = AnCg_nAcc_DLPFC_directions)

male_Ramaker_magma <- read_csv(here("ProcessedData", "RamakerEtAl", "CompleteMaleRamakerTable.csv"))
male_Ramaker_magma %<>% mutate(gene_symbol = gsub("C([X0-9]+)ORF([0-9]+)", "C\\1orf\\2", gene_symbol))
male_Ramaker_magma  %<>% right_join(magma %>% select(Ramaker_genes) %>% distinct(), by = c( 'gene_symbol' = 'Ramaker_genes')) %>% na.omit() %>% write_csv(here("ProcessedData","RamakerEtAl", "CompleteMaleRamakerTableMagma.csv"))
male_Ramaker_magma %<>% RamakerMetaAnalysis(regions)
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
