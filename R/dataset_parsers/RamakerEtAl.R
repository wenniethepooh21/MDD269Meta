library(readr)
library(tidyr)
library(GEOquery)
library(magrittr)
library(dplyr)
library(here)
                                                                #This script runs all meta-analyses on the Ramaker transcriptomic dataset 
# detach("package:here", unload=TRUE)
# setwd('../../../school/thesis/')
# library(here)

#based on signature.R by the biojupies team
#https://github.com/MaayanLab/biojupies-plugins/blob/1024a6ed702ad8b0958d4ccdd2afe89cbe493a51/library/core_scripts/signature/signature.R
#from https://amp.pharm.mssm.edu/biojupies/notebook/tlNmgbMxY (just load biojupies with GSE80655)

#Method for getting sequence counts
#"Raw RNA-seq data for GEO dataset GSE80655 was downloaded from the SRA database (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80655) and quantified to gene-level counts using the ARCHS4 pipeline (Lachmann et al., 2017). Gene counts were downloaded from the ARCHS4 gene expression matrix v6. For more information about ARCHS4, as well as free access to the quantified gene expression matrix, visit the project home page at the following URL: http://amp.pharm.mssm.edu/archs4/download.html."
#(Biojupies)

#parsing of their GEO SOFT format file into R object
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

###########################################################
###### REGULAR META-ANALYSIS (FULL, FEMALE AND MALE) ######
###########################################################
#This file holds all functions used for meta-analyses - call this file
source(here("R/transcriptomic_meta/Ramaker_Meta_Analysis.R"))
#get the unique brain regions
regions <- unique(metadata$`brain region`)

#create empty tibble for data to be populated 
female_results <- tibble()
#full female Ramaker data results
R_female_summary_results <- RamakerDEModel(female_metadata, read_counts, rawcount_dataframe, regions, female_results)
R_female_summary_results %<>% mutate(gene_symbol = gsub("C([X0-9]+)ORF([0-9]+)", "C\\1orf\\2", gene_symbol), sex = "female") #change from upper case to lower case for open reading frame genes
female_summary_ramaker <- R_female_summary_results %>% write_csv(here("Processed_Data/RamakerEtAl/CompleteFemaleRamakerTable.csv"))#write out for easier access in other analyses
#perform meta-analysis on female data
R_female_summary_results %<>% RamakerMetaAnalysis(regions)
R_female_summary_results %<>% rename(AnCg_nAcc_DLPFC_Female_directions = AnCg_nAcc_DLPFC_directions)
#Save full female meta-analysis results
R_female_summary_results %>% write_csv(here("Processed_Data/RamakerEtAl/FemaleRamakerTable.csv"))


#create empty tibble for data to be populated 
male_results <- tibble()
#full male Ramaker data results
R_male_summary_results <- RamakerDEModel(male_metadata, read_counts, rawcount_dataframe, regions, male_results)
R_male_summary_results %<>% mutate(gene_symbol = gsub("C([X0-9]+)ORF([0-9]+)", "C\\1orf\\2", gene_symbol), sex = "male") #change from upper case to lower case for open reading frame genes
male_summary_ramaker <- R_male_summary_results %>% write_csv(here("Processed_Data/RamakerEtAl/CompleteMaleRamakerTable.csv"))#write out for easier access in other analyses
#perform meta-analysis on male data
R_male_summary_results %<>% RamakerMetaAnalysis(regions)
R_male_summary_results %<>% rename(AnCg_nAcc_DLPFC_Male_directions = AnCg_nAcc_DLPFC_directions)
#Save full male meta-analysis results
R_male_summary_results %>% write_csv(here("Processed_Data/RamakerEtAl/MaleRamakerTable.csv"))

#combine female and male raw expression stats together
combined_sex_summary <- rbind(male_summary_ramaker, female_summary_ramaker) %>% write_csv("Processed_Data/RamakerEtAl/CombinedCompleteFemaleMaleRamakerTable.csv")

#Perform Ramaker meta-analysis functions in RamakerMetaAnalysis.R
#R_summary_results <- RamakerDEModel(fullmetadata, read_counts, rawcount_dataframe, regions, full_results)
R_summary_results <- combined_sex_summary 
#R_summary_results %>% write_csv(here("Processed_Data/RamakerEtAl/CompleteRamakerTable.csv")) #write out for easier access in other analyses (cortical)
#perform meta-analysis on full (female and male all brain regions) data

R_summary_results %<>% RamakerMetaAnalysis(regions)



#merge all gender directions into summary_results table full meta analysis visualization
Ramaker_summary <- left_join(R_female_summary_results %>% select(gene_symbol,AnCg_nAcc_DLPFC_Female_directions), R_male_summary_results %>% select(gene_symbol,AnCg_nAcc_DLPFC_Male_directions ))
Ramaker_summary %<>% unite(AnCg.F_nAcc.F_DLPFC.F_AnCg.M_nAcc.M_DLPFC.M, AnCg_nAcc_DLPFC_Female_directions, AnCg_nAcc_DLPFC_Male_directions, sep = "")
Ramaker_summary %<>% left_join(R_summary_results%>% select(-sex,-AnCg_nAcc_DLPFC_directions) %>% distinct())
#Save full meta-analysis results
Ramaker_summary %>% write_csv(here("Processed_Data/RamakerEtAl/FullRamakerTable.csv")) #used for combining across transcriptomic studies 

########################################################
###### CORTICAL ANALYSIS (FULL, FEMALE AND MALE) ######
########################################################
#extract cortical data and re-run meta-analysis -- don't have to re-run model creation
Ramaker_cortical<- read_csv(here("Processed_Data/RamakerEtAl/CombinedCompleteFemaleMaleRamakerTable.csv"))
#Extract the cortical region data & run analysis
Ramaker_cortical %<>% filter(target_region != "nAcc")
cortical_regions <- Ramaker_cortical %>% select(target_region) %>% distinct() %>% pull()
#run the meta-analysis on only the cortical regions sampled
Ramaker_cortical %<>% RamakerMetaAnalysis(cortical_regions)

#Extract the cortical region data & run analysis on female data
female_ramaker_cortical <- read_csv(here("Processed_Data/RamakerEtAl/CompleteFemaleRamakerTable.csv"))
female_ramaker_cortical %<>% filter(target_region != "nAcc")
#run the meta-analysis on only the cortical regions sampled in female data
female_ramaker_cortical %<>% RamakerMetaAnalysis(cortical_regions)
female_ramaker_cortical %<>% rename(AnCg_DLPFC_Female_directions = AnCg_DLPFC_directions)
female_ramaker_cortical %>% write_csv(here("Processed_Data/RamakerEtAl/CorticalFemaleRamakerTable.csv"))

#Extract the cortical region data & run analysis on male data
male_ramaker_cortical<- read_csv(here("Processed_Data/RamakerEtAl/CompleteMaleRamakerTable.csv"))
male_ramaker_cortical %<>% filter(target_region != "nAcc")
#run the meta-analysis on only the cortical regions sampled in male data
male_ramaker_cortical %<>% RamakerMetaAnalysis(cortical_regions)
male_ramaker_cortical %<>% rename(AnCg_DLPFC_Male_directions = AnCg_DLPFC_directions)
male_ramaker_cortical %>% write_csv(here("Processed_Data/RamakerEtAl/CorticalMaleRamakerTable.csv"))

#merge all directions from male and female data to visualize cortical directions across sexes
cortical_summary <- left_join(female_ramaker_cortical %>% select(gene_symbol,AnCg_DLPFC_Female_directions), male_ramaker_cortical %>% select(gene_symbol,AnCg_DLPFC_Male_directions ))
cortical_summary %<>% unite("AnCg.F_DLPFC.F_AnCg.M_DLPFC.M", AnCg_DLPFC_Female_directions, AnCg_DLPFC_Male_directions, sep = "")
cortical_summary %<>% left_join(Ramaker_cortical%>% select(-sex, -AnCg_DLPFC_directions) %>% distinct())
#Save full cortical meta-analysis results 
cortical_summary %>% write_csv(here("Processed_Data/RamakerEtAl/CorticalRamakerTable.csv"))


########################################################
###### SUBCORTICAL ANALYSIS (FULL, FEMALE AND MALE) ######
########################################################
#extract cortical data and re-run meta-analysis -- don't have to re-run model creation
Ramaker_subcortical<- read_csv(here("Processed_Data/RamakerEtAl/CombinedCompleteFemaleMaleRamakerTable.csv"))
#Extract the cortical region data & run analysis
Ramaker_subcortical %<>% filter(target_region == "nAcc")
subcortical_regions <- Ramaker_subcortical %>% select(target_region) %>% distinct() %>% pull()
#run the meta-analysis on only the cortical regions sampled
Ramaker_subcortical %<>% RamakerMetaAnalysis(subcortical_regions)

#Extract the cortical region data & run analysis on female data
female_ramaker_subcortical <- read_csv(here("Processed_Data/RamakerEtAl/CompleteFemaleRamakerTable.csv"))
female_ramaker_subcortical %<>% filter(target_region == "nAcc")
#run the meta-analysis on only the cortical regions sampled in female data
female_ramaker_subcortical %<>% RamakerMetaAnalysis(subcortical_regions)
female_ramaker_subcortical %<>% rename(nAcc_Female_directions = nAcc_directions)
female_ramaker_subcortical %>% write_csv(here("Processed_Data/RamakerEtAl/SubcorticalFemaleRamakerTable.csv"))

#Extract the cortical region data & run analysis on male data
male_ramaker_subcortical<- read_csv(here("Processed_Data/RamakerEtAl/CompleteMaleRamakerTable.csv"))
male_ramaker_subcortical %<>% filter(target_region == "nAcc")
#run the meta-analysis on only the cortical regions sampled in male data
male_ramaker_subcortical %<>% RamakerMetaAnalysis(subcortical_regions)
male_ramaker_subcortical %<>% rename(nAcc_Male_directions = nAcc_directions)
male_ramaker_subcortical %>% write_csv(here("Processed_Data/RamakerEtAl/SubcorticalMaleRamakerTable.csv"))

#merge all directions from male and female data to visualize cortical directions across sexes
subcortical_summary <- left_join(female_ramaker_subcortical %>% select(gene_symbol,nAcc_Female_directions), male_ramaker_subcortical %>% select(gene_symbol,nAcc_Male_directions ))
subcortical_summary %<>% unite("nAcc.F_nAcc.M", nAcc_Female_directions, nAcc_Male_directions, sep = "")
subcortical_summary %<>% left_join(Ramaker_subcortical %>% select(-sex, -nAcc_directions) %>% distinct())
#Save full cortical meta-analysis results 
subcortical_summary %>% write_csv(here("Processed_Data/RamakerEtAl/SubcorticalRamakerTable.csv"))

#############################################
###### SEX-INTERACTION FULL ANALYSIS ######
############################################
#### Full meta-analysis on flipped male expression data
#Read in the post-model female data 
full_female_results <- combined_sex_summary %>% filter(sex == "female")
full_male_results <- combined_sex_summary %>% filter(sex == "male")
#flip the expression levels of each gene in each brain region 
full_male_results_flip <- full_male_results %>% mutate(t = t*-1)
#add sex for unique identifier
full_male_results_flip %<>% mutate(sex = "male")
#merge male summary results with female summary results into one table to perform calculations
full_flipped <- rbind(full_male_results_flip, full_female_results)

#Perform the meta-analysis on the female and flipped male data across all brain regions
full_flipped %<>% RamakerMetaAnalysis(regions)

#Read in to get the original directions data from the full meta-analysis
Ramaker_directions <- read_csv(here("Processed_Data/RamakerEtAl/FullRamakerTable.csv")) %>% select(1:2)

#Keep flipped meta-analysis results 
full_flipped %<>% select(-sex, -AnCg_nAcc_DLPFC_directions) 
#join the direction data with the meta-analysis results on the flipped male expression data across all brain regions
full_flipped_summary <- Ramaker_directions %>% left_join(full_flipped) %>% distinct()
full_flipped_summary %>% write_csv(here("Processed_Data/RamakerEtAl/FullRamakerTable_flipped.csv"))

########################################################
###### SEX-INTERACTION CORTICAL ANALYSIS ######
########################################################
#get female post model differential expression data for cortical regions
cortical_female_results <- combined_sex_summary  %>% filter(target_region != "nAcc") %>% filter(sex == "female")
#get male post model differential expression data for cortical regions
cortical_male_results <- combined_sex_summary  %>% filter(target_region != "nAcc") %>% filter(sex == "male")
#flip the expression levels of each gene in each brain region 
cortical_male_results_flip <- cortical_male_results %>% mutate(t = t*-1)
#merge male cortical results with female cortical results into one table to perform meta-analysis calculations
cortical_flipped <- rbind(cortical_male_results_flip, cortical_female_results)

#Perform the meta-analysis on the female and flipped male data across all cortical brain regions
cortical_flipped %<>% RamakerMetaAnalysis(cortical_regions)

#Read in to get the original directions data from the cortical meta-analysis
Ramaker_cortical_directions <- read_csv(here("Processed_Data/RamakerEtAl/CorticalRamakerTable.csv")) %>% select(1:2)

#Keep flipped meta-analysis results 
cortical_flipped %<>% select(-sex, -AnCg_DLPFC_directions)
#join the direction data with the meta-analysis results on the flipped male expression data across all cortical brain regions
cortical_flipped_summary <- Ramaker_cortical_directions %>% left_join(cortical_flipped) %>% distinct()
cortical_flipped_summary %>% write_csv(here("Processed_Data/RamakerEtAl/CorticalRamakerTable_flipped.csv"))

########################################################
###### SEX-INTERACTION SUBORTICAL ANALYSIS ######
########################################################
#get female post model differential expression data for subcortical regions
subcortical_female_results <- combined_sex_summary  %>% filter(target_region == "nAcc" ) %>% filter(sex == "female")
#get male post model differential expression data for subcortical regions
subcortical_male_results <- combined_sex_summary  %>% filter(target_region == "nAcc") %>% filter(sex == "male")
#flip the expression levels of each gene in each brain region 
subcortical_male_results_flip <- subcortical_male_results %>% mutate(t = t*-1)
#merge male cortical results with female cortical results into one table to perform meta-analysis calculations
subcortical_flipped <- rbind(subcortical_male_results_flip, subcortical_female_results)

#Perform the meta-analysis on the female and flipped male data across all cortical brain regions
subcortical_flipped %<>% RamakerMetaAnalysis(subcortical_regions)

#Read in to get the original directions data from the cortical meta-analysis
Ramaker_subcortical_directions <- read_csv(here("Processed_Data/RamakerEtAl/SubcorticalRamakerTable.csv")) %>% select(1:2)

#Keep flipped meta-analysis results 
subcortical_flipped %<>% select(-sex, -nAcc_directions)
#join the direction data with the meta-analysis results on the flipped male expression data across all cortical brain regions
subcortical_flipped_summary <- Ramaker_subcortical_directions %>% left_join(subcortical_flipped) %>% distinct()
subcortical_flipped_summary %>% write_csv(here("Processed_Data/RamakerEtAl/SubcorticalRamakerTable_flipped.csv"))


###############################################################
###### GENOME PERCENTILE RANKING FOR ALL ABOVE ANALYSES ######
##############################################################
#filter Ramaker table for magma genes
ramaker_magma <- read_csv(here("Raw_Data/HowardEtAl/FullMagmaGenes.csv")) %>% select(Ramaker_genes) %>% distinct() %>% na.omit()
#load script that holds the genome percentile ranking function used by all transcriptomic studies
source(here("R/transcriptomic_meta/Percentile_Rank_Analysis.R"))
#Run genome percentile ranking analysis - calculates the percentage of genes that have a smaller meta p-value than the current gene on our full meta-analysis

Ramaker_summary%<>% right_join(ramaker_magma, by = c('gene_symbol' = 'Ramaker_genes'))
Ramaker_summary %<>% getRank()
Ramaker_summary %>% write_csv(here("Processed_Data/RamakerEtAl/FullRamakerTableMagma.csv"))

#Run genome percentile ranking analysis - calculates the percentage of genes that have a smaller meta p-value than the current gene on our female meta-analysis
R_female_summary_results%<>% right_join(ramaker_magma, by = c('gene_symbol' = 'Ramaker_genes'))
R_female_summary_results %<>% getRank()
R_female_summary_results %>% write_csv(path = here("Processed_Data/RamakerEtAl/FemaleRamakerTableMagma.csv"))

#Run genome percentile ranking analysis - calculates the percentage of genes that have a smaller meta p-value than the current gene on our male meta-analysis
R_male_summary_results%<>% right_join(ramaker_magma, by = c('gene_symbol' = 'Ramaker_genes'))
R_male_summary_results %<>% getRank()
R_male_summary_results %>% write_csv(path = here("Processed_Data/RamakerEtAl/MaleRamakerTableMagma.csv"))

#Run genome percentile ranking analysis - calculates the percentage of genes that have a smaller meta p-value than the current gene on our cortical meta-analysis
cortical_summary%<>% right_join(ramaker_magma, by = c('gene_symbol' = 'Ramaker_genes'))
cortical_summary %<>% getRank()
cortical_summary %>% write_csv(here("Processed_Data/RamakerEtAl/CorticalRamakerTableMagma.csv"))

#Run genome percentile ranking analysis - calculates the percentage of genes that have a smaller meta p-value than the current gene on our subcortical meta-analysis
subcortical_summary%<>% right_join(ramaker_magma, by = c('gene_symbol' = 'Ramaker_genes'))
subcortical_summary %<>% getRank()
subcortical_summary %>% write_csv(here("Processed_Data/RamakerEtAl/SubcorticalRamakerTableMagma.csv"))

#Run genome percentile ranking analysis - calculates the percentage of genes that have a smaller meta p-value than the current gene on our sex-interaction full meta-analysis
full_flipped_summary%<>% right_join(ramaker_magma, by = c('gene_symbol' = 'Ramaker_genes'))
full_flipped_summary %<>% getRank()
full_flipped_summary %>% write_csv(here("Processed_Data/RamakerEtAl/FullRamakerTableMagma_flipped.csv"))

#Run genome percentile ranking analysis - calculates the percentage of genes that have a smaller meta p-value than the current gene on our sex-interaction cortical meta-analysis
cortical_flipped_summary%<>% right_join(ramaker_magma, by = c('gene_symbol' = 'Ramaker_genes'))
cortical_flipped_summary %<>% getRank()
cortical_flipped_summary %>% write_csv(here("Processed_Data/RamakerEtAl/CorticalRamakerTableMagma_flipped.csv"))

#Run genome percentile ranking analysis - calculates the percentage of genes that have a smaller meta p-value than the current gene on our sex-interaction subcortical meta-analysis
subcortical_flipped_summary%<>% right_join(ramaker_magma, by = c('gene_symbol' = 'Ramaker_genes'))
subcortical_flipped_summary %<>% getRank()
subcortical_flipped_summary %>% write_csv(here("Processed_Data/RamakerEtAl/SubcorticalRamakerTableMagma_flipped.csv"))
