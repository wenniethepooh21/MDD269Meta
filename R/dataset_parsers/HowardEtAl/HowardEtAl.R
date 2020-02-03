library(mygene)
library(magrittr)
library(mygene)
library(dplyr)
library(here)
#Using Howard 102 gene list 

#Import the Howard data 
# Pre-process Howard gene list to match the genes used in the other studies 
Howard <- read_xlsx(here("Raw_Data/HowardEtAl/433367-22.xlsx"), skip = 1)
Howard %<>% select(`Gene Name`, `P-value`)
Howard %<>% rename(gene_symbol = `Gene Name`, Howard_pvalue = `P-value`)

#Get the full gene name from the gene symbol
gene_equivalent <- read_csv(here("Processed_Data/HowardEtAl/Gene_Reference_list_269.csv"))
full_gene_names <- queryMany(gene_equivalent$Updated_Gene_Names, scope = "symbol", species = "human") %>% as_tibble() %>% distinct(query,name)

full_gene_names %<>% left_join(gene_equivalent %>% select(-gene_symbol), by = c('query' = 'Updated_Gene_Names'))
Howard_Table <- Howard %>% left_join(full_gene_names, by = c('gene_symbol' = 'query'))
Howard_Table %<>% rename(gene_name = name)
#do these manually
Howard_Table %>% filter(is.na(gene_name))
  # Howard_Table %>% mutate(gene_name = if_else(gene_symbol == "C16orf45", "BMERB Domain Containing 1", gene_name))
  # Howard_Table %>% mutate(gene_name = if_else(gene_symbol == "C16orf45", "BMERB Domain Containing 1", gene_name))
  # Howard_Table %>% mutate(gene_name = if_else(gene_symbol == "C16orf45", "BMERB Domain Containing 1", gene_name))
  # Howard_Table %>% mutate(gene_name = if_else(gene_symbol == "C16orf45", "BMERB Domain Containing 1", gene_name))
  # Howard_Table %>% mutate(gene_name = if_else(gene_symbol == "C16orf45", "BMERB Domain Containing 1", gene_name))
  # Howard_Table %>% mutate(gene_name = if_else(gene_symbol == "C16orf45", "BMERB Domain Containing 1", gene_name))
Howard_Table %>% write_csv(here("Processed_Data/HowardEtAl/fullHowardTable.csv"))

