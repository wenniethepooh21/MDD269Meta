library(homologene)
library(here)
library(magrittr)
library(dplyr)

#This script matches the howard 269 genes to the max brain region and cell type taxon that expressed it
howard_regions <- read_csv(here("Processed_Data/HowardEtAl/HowardRegions_four.csv"))

howard_cell_types_z <- read_csv(here("Processed_Data/ZeiselEtAl/HowardCellTypes_zscore.csv"))
howard_region_cell_types_z <- howard_regions %>% left_join(howard_cell_types_z)
howard_region_cell_types_z %>% write_csv(here("Processed_Data/HowardEtAl/HowardRegionsPolygenicCellTypesZScore_four.csv"))


