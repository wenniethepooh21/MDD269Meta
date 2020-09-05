library(readr)
library(magrittr)
library(dplyr)
library(tidyr)
library(here)
library(googledrive)
library(googlesheets4)

##############UPLOAD BRAIN REGION TO GOOGLE SHEETS 
tissue_expected_probs_full <- read_csv(here('Results', 'supplementary_tables', 'hypergeometric_brain_regions_tissue_hyper_expected_four_full.csv'))
tissue_expected_probs_slim <- read_csv(here('Results', 'supplementary_tables', 'hypergeometric_brain_regions_tissue_hyper_expected_four.csv'))

gs4_auth(token = drive_token())

region <- drive_get("~/Thesis/Manuscript/gs_tables/Hypergeometric/tissue_hyper_expected_four_full_new")
if(nrow(region) != 0) {
  drive_rm(region)
}
#create the google worksheet
region <- gs4_create("tissue_hyper_expected_four_full_new",sheets = c('hypergeometric_brain_regions'))
sheet_write(tissue_expected_probs_full, region,  sheet = "hypergeometric_brain_regions")

drive_mv(file = "tissue_hyper_expected_four_full_new", path = "~/Thesis/Manuscript/gs_tables/Hypergeometric/")  # move Sheets file

###### SLIMMED
region <- drive_get("~/Thesis/Manuscript/gs_tables/Hypergeometric/tissue_hyper_expected_four")
if(nrow(region) != 0) {
  drive_rm(region)
}
#create the google worksheet
region <- gs4_create("tissue_hyper_expected_four",sheets = c('hypergeometric_brain_regions'))
sheet_write(tissue_expected_probs_slim, region,  sheet = "hypergeometric_brain_regions")

drive_mv(file = "tissue_hyper_expected_four", path = "~/Thesis/Manuscript/gs_tables/Hypergeometric/") 


######## UPLOAD CELL TYPES TO GOOGLE SHEETS" 

cell_expected_probs <- read_csv(here('Results', 'supplementary_tables', 'hypergeometric_cell_type_taxons_cell_expected_probs.csv'))

gs4_auth(token = drive_token())

cells <- drive_get("~/Thesis/Manuscript/gs_tables/Hypergeometric/cell_hyper_expected")
if(nrow(cells) != 0) {
  drive_rm(cells)
}

#create the google worksheet
cells <- gs4_create("cell_hyper_expected",sheets = c('hypergeometric_cell_type_taxons'))
sheet_write(cell_expected_probs, cells,  sheet = "hypergeometric_cell_type_taxons")

drive_mv(file = "cell_hyper_expected", path = "~/Thesis/Manuscript/gs_tables/Hypergeometric/")  # move Sheets file

####### UPLOAD CNS CELL TYPES TO GOOGLE SHEETS

cns_cell_expected_probs <- read_csv(here('Results', 'supplementary_tables', 'hypergeometric_cns_cell_type_taxons_cell_expected_probs.csv'))

gs4_auth(token = drive_token())

cells <- drive_get("~/Thesis/Manuscript/gs_tables/Hypergeometric/cns_cell_hyper_expected")
if(nrow(cells) != 0) {
  drive_rm(cells)
}

#create the google worksheet
cells <- gs4_create("cns_cell_hyper_expected",sheets = c('hypergeometric_cns_cell_type_taxons'))
sheet_write(cns_cell_expected_probs, cells,  sheet = "hypergeometric_cns_cell_type_taxons")

drive_mv(file = "cns_cell_hyper_expected", path = "~/Thesis/Manuscript/gs_tables/Hypergeometric/")  # move Sheets file


###### UPLOAD CHOLINERGIC AND MONOAMINERGIC AND ENTERIC NEURONS ANALYSIS USING 265 
#upload to google drive
cell_expected_probs <- read_csv(here("Results/supplementary_tables/hypergeometric_cell_type_taxon_265.csv"))
cell_expected_probs_ent_chol <- read_csv(here("Results/supplementary_tables/hypergeometric_cell_type_taxon_265_chol_ent.csv"))

gs4_auth(token = drive_token())
cells_265 <- drive_get("~/Thesis/Manuscript/Supplement_Tables/hypergeometric_cell_type_taxa_265")
if(nrow(cells_265) != 0) {
  drive_rm(cells_265)
}
#create the google worksheet
cells_265 <- gs4_create("hypergeometric_cell_type_taxa_265",sheets = c('hypergeometric_cell_type_taxons_265', 'hypergeometric_cell_type_taxons_265_enteric_cholinergic_neurons'))
sheet_write(cell_expected_probs, cells_265,  sheet = "hypergeometric_cell_type_taxons_265")
sheet_write(cell_expected_probs_ent_chol, cells_265,  sheet = "hypergeometric_cell_type_taxons_265_enteric_cholinergic_neurons")

drive_mv(file = "hypergeometric_cell_type_taxa_265", path = "~/Thesis/Manuscript/Supplement_Tables/")  # move Sheets file





######UPLOAD META=ANALYSES TO GOOGLE SHEETS##########

fullTable <- read_csv( here('Results', 'Tables', 'Meta_Analysis', "Full_Meta_Analysis.csv"))
femaleTable <- read_csv(here('Results', 'Tables', 'Meta_Analysis', 'Female_Meta_Analysis.csv'))
maleTable <- read_csv(here('Results', 'Tables', 'Meta_Analysis', 'Male_Meta_Analysis.csv'))
corticalTable <- read_csv( here('Results', 'Tables', 'Meta_Analysis', 'Cortical_Meta_Analysis.csv'))
subcorticalTable <- read_csv( here('Results', 'Tables', 'Meta_Analysis', 'Subcortical_Meta_Analysis.csv'))
fullTable_Flip <- read_csv(here('Results', 'Tables', 'Meta_Analysis','Sex_interaction_Full_Meta_Analysis.csv'))
corticalTable_Flip <- read_csv( here('Results', 'Tables', 'Meta_Analysis', 'Sex_interaction_Cortical_Meta_Analysis.csv'))
subcorticalTable_Flip <- read_csv(here('Results', 'Tables', 'Meta_Analysis', 'Sex_interaction_Subcortical_Meta_Analysis.csv'))

fullTableRank <- read_csv( here('Results', 'Tables', 'Genome_Percentile_Rank', "Full_Meta_Analysis_Rank.csv"))
femaleTableRank <- read_csv(here('Results', 'Tables', 'Genome_Percentile_Rank', 'Female_Meta_Analysis_Rank.csv'))
maleTableRank <- read_csv(here('Results', 'Tables', 'Genome_Percentile_Rank', 'Male_Meta_Analysis_Rank.csv'))
corticalTableRank <- read_csv(here('Results', 'Tables', 'Genome_Percentile_Rank', 'Cortical_Meta_Analysis_Rank.csv'))
subcorticalTableRank <- read_csv( here('Results', 'Tables', 'Genome_Percentile_Rank', 'Subcortical_Meta_Analysis_Rank.csv'))
fullTableRank_Flip <- read_csv( here('Results', 'Tables', 'Genome_Percentile_Rank', 'Sex_interaction_Full_Meta_Analysis_Rank.csv'))
corticalTableRank_Flip <- read_csv(here('Results', 'Tables', 'Genome_Percentile_Rank', 'Sex_interaction_Cortical_Meta_Analysis_Rank.csv'))
subcorticalTableRank_Flip <- read_csv( here('Results', 'Tables', 'Genome_Percentile_Rank', 'Sex_interaction_Subcortical_Meta_Analysis_Rank.csv'))


#FULL DATA SHEETS
#Access googlesheets to upload the tables online for an interactive experience

gs4_auth(token = drive_token())

ma <- drive_get("~/Thesis/Manuscript/gs_tables/Meta_Analysis_Full_Tables/Meta_Analysis")
if(nrow(ma) != 0) {
  drive_rm(ma)
}
#create the google worksheet
ma <- gs4_create("Meta_Analysis", sheets = c('Full_Meta_Analysis', 'Male_Meta_Analysis','Female_Meta_Analysis','Cortical_Meta_Analysis','Subcortical_Meta_Analysis', 'Sex_interaction_Full_Meta_Analysis', 'Sex_interaction_Cortical_Meta_Analysis', 'Sex_interaction_Subcortical_Meta_Analysis'))
sheet_write(fullTable, ma, sheet = "Full_Meta_Analysis")
sheet_write(femaleTable, ma, sheet = 'Female_Meta_Analysis')
sheet_write(maleTable, ma,'Male_Meta_Analysis')
sheet_write(corticalTable, ma, 'Cortical_Meta_Analysis')
sheet_write(subcorticalTable, ma, 'Subcortical_Meta_Analysis')
sheet_write(fullTable_Flip, ma, 'Sex_interaction_Full_Meta_Analysis')
sheet_write(corticalTable_Flip, ma, 'Sex_interaction_Cortical_Meta_Analysis')
sheet_write(subcorticalTable_Flip, ma, 'Sex_interaction_Subcortical_Meta_Analysis')

drive_mv(file = "Meta_Analysis", path = "~/Thesis/Manuscript/gs_tables/Meta_Analysis_Full_Tables/")  # move Sheets file

mar <- drive_get("Thesis/Manuscript/Tables/gs_tables/Meta_Analysis_Full_Tables/Genome_Percentile_Rank")
if(nrow(mar) != 0) {
  drive_rm(mar)
}
#create the google worksheet
mar <- gs4_create("Genome_Percentile_Rank", sheets = c('Full_Meta_Analysis_Rank', 'Male_Meta_Analysis_Rank','Female_Meta_Analysis_Rank', 'Cortical_Meta_Analysis_Rank','Subcortical_Meta_Analysis_Rank', 'Sex_interaction_Full_Meta_Analysis_Rank',  'Sex_interaction_Cortical_Meta_Analysis_Rank',  'Sex_interaction_Subcortical_Meta_Analysis_Rank'))
sheet_write(fullTableRank, mar, sheet = "Full_Meta_Analysis_Rank")
sheet_write(femaleTableRank, mar, sheet = 'Female_Meta_Analysis_Rank')
sheet_write(maleTableRank, mar,'Male_Meta_Analysis_Rank')
sheet_write(corticalTableRank, mar, 'Cortical_Meta_Analysis_Rank')
sheet_write(subcorticalTableRank, mar, 'Subcortical_Meta_Analysis_Rank')
sheet_write(fullTableRank_Flip, mar, 'Sex_interaction_Full_Meta_Analysis_Rank')
sheet_write(corticalTableRank_Flip, mar, 'Sex_interaction_Cortical_Meta_Analysis_Rank')
sheet_write(subcorticalTableRank_Flip, mar, 'Sex_interaction_Subcortical_Meta_Analysis_Rank')

drive_mv(file = "Genome_Percentile_Rank", path = "~/Thesis/Manuscript/gs_tables/Meta_Analysis_Full_Tables/")  # move Sheets file


ma <- drive_get("~/Thesis/Manuscript/gs_tables/Meta_Analysis_Slim_Tables/Official_Meta_Analysis")
if(nrow(ma) != 0) {
  drive_rm(ma)
}
#create the google worksheet
ma <- gs4_create("Official_Meta_Analysis", sheets = c('First_model_Full_Meta_Analysis', 'First_model_Male_Meta_Analysis','First_model_Female_Meta_Analysis','First_model_Cortical_Meta_Analysis','First_model_Subcortical_Meta_Analysis','Sex_interaction_Full_Meta_Analysis', 'Sex_interaction_Cortical_Meta_Analysis', 'Sex_interaction_Subcortical_Meta_Analysis'))
sheet_write(fullTable %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction,min_p_study, meta_p, signed_log_meta, Bonferroni_meta_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma, sheet = "First_model_Full_Meta_Analysis")
sheet_write(maleTable %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir,  meta_direction,min_p_study, meta_p, signed_log_meta, Bonferroni_meta_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma,'First_model_Male_Meta_Analysis')
sheet_write(femaleTable %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction,min_p_study, meta_p, signed_log_meta, Bonferroni_meta_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma, sheet = 'First_model_Female_Meta_Analysis')
sheet_write(corticalTable %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction,min_p_study, meta_p, signed_log_meta, Bonferroni_meta_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma, 'First_model_Cortical_Meta_Analysis')
sheet_write(subcorticalTable %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir,  meta_direction,min_p_study, meta_p, signed_log_meta, Bonferroni_meta_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma, 'First_model_Subcortical_Meta_Analysis')
sheet_write(fullTable_Flip  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir,  meta_direction,min_p_study, meta_p, signed_log_meta, Bonferroni_meta_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma, sheet = "Sex_interaction_Full_Meta_Analysis")
sheet_write(corticalTable_Flip  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir,  meta_direction,min_p_study, meta_p, signed_log_meta, Bonferroni_meta_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma, 'Sex_interaction_Cortical_Meta_Analysis')
sheet_write(subcorticalTable_Flip  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction,min_p_study, meta_p, signed_log_meta, Bonferroni_meta_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma, 'Sex_interaction_Subcortical_Meta_Analysis')

drive_mv(file = "Official_Meta_Analysis", path = "~/Thesis/Manuscript/gs_tables/Meta_Analysis_Slim_Tables/")  # move Sheets file

mar <- drive_get("Thesis/Manuscript/gs_tables/Meta_Analysis_Slim_Tables/Official_Genome_Percentile_Rank")
if(nrow(mar) != 0) {
  drive_rm(mar)
}
#create the google worksheet
mar <- gs4_create("Official_Genome_Percentile_Rank", sheets = c('First_model_Subcortical_Meta_Analysis_Rank', 'Sex_interaction_Cortical_Meta_Analysis_Rank', 'First_model_Male_Meta_Analysis_Rank','First_model_Full_Meta_Analysis_Rank', 'First_model_Female_Meta_Analysis_Rank', 'First_model_Cortical_Meta_Analysis_Rank','Sex_interaction_Full_Meta_Analysis_Rank','Sex_interaction_Subcortical_Meta_Analysis_Rank'))
sheet_write(subcorticalTableRank  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon, cell_type_taxon_zscore, DE_Prior_Rank), mar, 'First_model_Subcortical_Meta_Analysis_Rank')
sheet_write(corticalTableRank_Flip  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), mar, 'Sex_interaction_Cortical_Meta_Analysis_Rank')
sheet_write(maleTableRank  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), mar,'First_model_Male_Meta_Analysis_Rank')
sheet_write(fullTableRank  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon, cell_type_taxon_zscore,DE_Prior_Rank), mar, sheet = "First_model_Full_Meta_Analysis_Rank")
sheet_write(femaleTableRank  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), mar, sheet = 'First_model_Female_Meta_Analysis_Rank')
sheet_write(corticalTableRank  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), mar, 'First_model_Cortical_Meta_Analysis_Rank')
sheet_write(fullTableRank_Flip  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), mar, sheet = "Sex_interaction_Full_Meta_Analysis_Rank")
sheet_write(subcorticalTableRank_Flip  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon, cell_type_taxon_zscore, DE_Prior_Rank), mar, 'Sex_interaction_Subcortical_Meta_Analysis_Rank')

drive_mv(file = "Official_Genome_Percentile_Rank", path = "~/Thesis/Manuscript/gs_tables/Meta_Analysis_Slim_Tables/")  # move Sheets file



ma <- drive_get("~/Thesis/Manuscript/gs_tables/Meta_Analysis_Slim_Tables/Full_Meta_Analysis")
if(nrow(ma) != 0) {
  drive_rm(ma)
}
#create the google worksheet
ma <- gs4_create("Full_Meta_Analysis", sheets = c('Table S1', 'Table S2','Table S3','Table S4','Table S5','Table S6', 'Table S7', 'Table S8',
                                                      'Table S9', 'Table S10', 'Table S11','Table S12', 'Table S13', 'Table S14','Table S15','Table S16'))
sheet_write(fullTable %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction,min_p_study, meta_p, signed_log_meta, Bonferroni_meta_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma, sheet = "Table S1")
sheet_write(maleTable %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir,  meta_direction,min_p_study, meta_p, signed_log_meta, Bonferroni_meta_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma,'Table S2')
sheet_write(femaleTable %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction,min_p_study, meta_p, signed_log_meta, Bonferroni_meta_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma, sheet = 'Table S3')
sheet_write(corticalTable %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction,min_p_study, meta_p, signed_log_meta, Bonferroni_meta_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma, 'Table S4')
sheet_write(subcorticalTable %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir,  meta_direction,min_p_study, meta_p, signed_log_meta, Bonferroni_meta_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma, 'Table S5')
sheet_write(fullTable_Flip  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir,  meta_direction,min_p_study, meta_p, signed_log_meta, Bonferroni_meta_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma, sheet = "Table S6")
sheet_write(corticalTable_Flip  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir,  meta_direction,min_p_study, meta_p, signed_log_meta, Bonferroni_meta_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma, 'Table S7')
sheet_write(subcorticalTable_Flip  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction,min_p_study, meta_p, signed_log_meta, Bonferroni_meta_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma, 'Table S8')
sheet_write(subcorticalTableRank  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon, cell_type_taxon_zscore, DE_Prior_Rank), ma, 'Table S9')
sheet_write(corticalTableRank_Flip  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma, 'Table S10')
sheet_write(maleTableRank  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma,'Table S11')
sheet_write(fullTableRank  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon, cell_type_taxon_zscore,DE_Prior_Rank), ma, sheet = "Table S12")
sheet_write(femaleTableRank  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma, sheet = 'Table S13')
sheet_write(corticalTableRank  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma, 'Table S14')
sheet_write(fullTableRank_Flip  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma, sheet = "Table S15")
sheet_write(subcorticalTableRank_Flip  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon, cell_type_taxon_zscore, DE_Prior_Rank), ma, 'Table S16')

drive_mv(file = "Full_Meta_Analysis", path = "~/Thesis/Manuscript/gs_tables/Meta_Analysis_Slim_Tables/")  # move Sheets file

