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
