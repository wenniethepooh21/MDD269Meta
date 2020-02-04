
#### General Script to run all meta-analysis to prioritize the 269 genes and identify which cell-type taxon and brain regions most highly express each gene 
#set the current working directory
setwd('../../../school/thesis/')
library(here)

#1. run the study-specific analyses (pre-processing included)
source(here("R/dataset_parsers/RamakerEtAl.R"))
source(here("R/dataset_parsers/DingEtAl.R"))
source(here("R/dataset_parsers/LabonteEtAl.R"))


#2. run pre=processing of the Howard GWAS gene-symbol matching 
#source(here("R/dataset_parsers/HowardEtAl.R"))

# Identify which brain region maximally express each gene
source(here("R/dataset_parsers/AllenEtAl.R"))

# Identify which cell type taxon maximally express each gene
source(here("R/dataset_parsers/ZeiselEtAl.R"))

# Merge the two results together into one table 
source(here("R/dataset_parsers/Howard_max_celltypes_regions.R"))

#run overall meta-analysis combining the results from the study-specific meta-analyses
source(here("R/transcriptomic_meta/mergeMetaOnSlims.R"))
#view the final product online in google drive! 

#Create heatmaps for visualization! 
source(here("R/visualization/meta_p_heatmap.R")) #generates a heatmap of the corrected meta p-values for each top gene across all meta-analyses (except ranking)
source(here("R/visualization/Labonte_heatmap.R")) #generates a heatmap of each top gene expression direction across each sampled brain region in Labonte data
source(here("R/visualization/Ramaker_heatmap.R")) #generates a heatmap of each top gene expression direction across each sampled brain region in Ramaker data
source(here("R/visualization/Ding_heatmap.R")) #generates a heatmap of each top gene expression direction across each sampled brain region in Ding data
