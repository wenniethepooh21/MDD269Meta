
#### General Script to run all meta-analysis to prioritize the 269 genes and identify which cell-type taxon and brain regions most highly express each gene 

setwd('../../../school/thesis/')
library(here)

#1. run the study-specific analyses (pre-processing included)
source(here("R/dataset_parsers/LabonteEtAl.R"))
source(here("R/dataset_parsers/DingEtAl.R"))
source(here("R/dataset_parsers/RamakerEtAl.R"))

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


#run visualization of expression data from each transcriptomic study
source(here("R/visualization/heatmaps.R"))