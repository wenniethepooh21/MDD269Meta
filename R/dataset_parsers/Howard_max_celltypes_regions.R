library(homologene)
library(here)
library(magrittr)
library(dplyr)

#This script matches the howard 269 genes to the max brain region and cell type taxon that expressed it
howard_regions <- read_csv(here("Processed_Data/HowardEtAl/HowardRegions_four.csv")) %>% distinct()
howard_cell_types <- read_csv(here("Processed_Data/ZeiselEtAl/HowardCellTypes.csv"))
#merge the two datasets together 
howard_region_cell_types <- howard_regions %>% left_join(howard_cell_types)
howard_region_cell_types %>% write_csv(here("Processed_Data/HowardEtAl/HowardRegionsPolygenicCellTypes_four.csv"))
