library(dplyr)
library(loomR)
library(here)

#tried accessing file directly from the website - doesn't work, must download on computer
#lfile <- connect(filename = "https://storage.googleapis.com/linnarsson-lab-loom/l6_r4.agg.loom", mode = "r+")
lfile <- connect(filename = here("data", "ZeiselEtAl", "l6_r4.agg.loom"))

zeisel_265 <- read_tsv(here('./l5_all.agg.tab'), col_names = FALSE)


lfile[["col_graphs"]]

exp_matrix <- as_tibble(t(as.matrix(lfile[["matrix"]][,])))
colnames(exp_matrix) <- lfile$col.attrs$TaxonomyRank4[]
exp_matrix %<>% mutate(Gene = lfile$row.attrs$Gene[]) %>% select(Gene, everything())
exp_matrix %>% write_csv(here("ProcessedData", "ZeiselEtAl", "39_types_zeisel.R"))

#deal with duplicate genes - avg the expression
exp_matrix %<>% group_by(Gene) %>% mutate_each(list(mean))  %>% distinct() %>% ungroup()
exp_matrix_long <- gather(exp_matrix,colnames(exp_matrix)[-1], key = "cell_types", value = "expression_levels")

#get only the max expression level for each gene and the associated cell type
max_cell_type_expression<-exp_matrix_long %>% group_by(Gene) %>% slice(which.max(expression_levels)) %>% ungroup()
max_cell_type_expression %>% write_csv(here("ProcessedData", "ZeiselEtAl", "max_cell_type_expression.csv"))

#---------
#remove PNS neurons 
cns_matrix <- as_tibble(t(as.matrix(lfile[["matrix"]][,])))
taxonomy_group_two <- lfile$col.attrs$TaxonomyRank2[] %>% as_tibble() 
taxonomy_group_two %<>% rename(Taxon_2 = value)
taxonomy_group_four <- lfile$col.attrs$TaxonomyRank4[] %>% as_tibble(name = "Taxon_4")
taxonomy_group_four %<>% rename(Taxon_4 = value)

Taxonomy_assignment <- cbind(taxonomy_group_two, taxonomy_group_four)
cns_colnames<- Taxonomy_assignment %>% filter(Taxon_2 != "PNS neurons") %>% select(Taxon_4) %>% pull()

exp_matrix_slim <- exp_matrix %>% select(Gene, cns_colnames)
exp_matrix_slim %>% write_csv(here("ProcessedData", "ZeiselEtAl", "cns_cell_types_zeisel.R"))


