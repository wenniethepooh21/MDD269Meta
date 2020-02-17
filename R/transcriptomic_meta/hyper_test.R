library(here)
library(readr)
library(dplyr)
library(magrittr)
library(stringr)

# hyper geometric probability of observing X successes 

#this function calculates the probability of seeing the (exact) expected value within the sample population
hyper_test <- function(sample_success, pop_success, pop_size, sample_size){
  p_greater <- 1-phyper(sample_success, pop_success,(pop_size - pop_success), sample_size)
  p_greater_equal <- phyper(sample_success-1, pop_success,(pop_size - pop_success), sample_size, lower.tail = FALSE)
  
  p_less <- phyper(sample_success-1, pop_success,(pop_size - pop_success), sample_size)
  p_less_equal <- 1- phyper(sample_success, pop_success,(pop_size - pop_success), sample_size, lower.tail = FALSE)
  
  p_equal <-p_greater_equal - p_greater
  
  return(p_greater_equal)
}

# #reset the 'here' package in case we're not in this directory
# detach("package:here", unload=TRUE)
# setwd('../../../Thesis/Datas/MDD44Characterize/')
# library(here)
# 
# #data used from all 6 donors
# Howard_full <- read_csv(here("ProcessedData", "HowardEtAl", "HowardRegionsPolygenicCellTypes.csv"))
# #data used from 4 donors
# Howard_full <- read_csv(here("ProcessedData", "HowardEtAl", "HowardRegionsPolygenicCellTypes_four.csv"))
# 
# #data used from 4 donors
# six_brains <- read_csv(here("ProcessedData", "AllenEtAl", "max_expression_in_four_donors.csv"))
# 
# cell_types <- read_csv(here("ProcessedData", "ZeiselEtAl", "max_cell_type_expression.csv"))
# cns_cell_types <- read_csv(here("ProcessedData", "ZeiselEtAl", "max_cns_cell_type_expression.csv"))
# #
# #


