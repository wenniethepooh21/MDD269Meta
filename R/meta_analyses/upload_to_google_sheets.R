library(readr)
library(magrittr)
library(dplyr)
library(tidyr)
library(here)
library(googledrive)
library(googlesheets4)


full <- read_csv(fullTable, here('Results', 'Tables', 'Meta_Analysis', "Full_Meta_Analysis.csv"))
female <- read_csv(femaleTable, here('Results', 'Tables', 'Meta_Analysis', 'Female_Meta_Analysis.csv'))
male <- read_csv(maleTable, here('Results', 'Tables', 'Meta_Analysis', 'Male_Meta_Analysis.csv'))
cortical <- read_csv(corticalTable, here('Results', 'Tables', 'Meta_Analysis', 'Cortical_Meta_Analysis.csv'))
subcortical <- read_csv(subcorticalTable, here('Results', 'Tables', 'Meta_Analysis', 'Subcortical_Meta_Analysis.csv'))
si_full <- read_csv(fullTable_Flip, here('Results', 'Tables', 'Meta_Analysis','Sex_interaction_Full_Meta_Analysis.csv'))
si_cortical <- read_csv(corticalTable_Flip, here('Results', 'Tables', 'Meta_Analysis', 'Sex_interaction_Cortical_Meta_Analysis.csv'))
si_subcortical <- read_csv(subcorticalTable_Flip, here('Results', 'Tables', 'Meta_Analysis', 'Sex_interaction_Subcortical_Meta_Analysis.csv'))

full_rank <- read_csv(fullTableRank, here('Results', 'Tables', 'Genome_Percentile_Rank', "Full_Meta_Analysis_Rank.csv"))
female_rank <- read_csv(femaleTableRank, here('Results', 'Tables', 'Genome_Percentile_Rank', 'Female_Meta_Analysis_Rank.csv'))
male_rank <- read_csv(maleTableRank, here('Results', 'Tables', 'Genome_Percentile_Rank', 'Male_Meta_Analysis_Rank.csv'))
cortical_rank <- read_csv(corticalTableRank, here('Results', 'Tables', 'Genome_Percentile_Rank', 'Cortical_Meta_Analysis_Rank.csv'))
subcortical_rank <- read_csv(subcorticalTableRank, here('Results', 'Tables', 'Genome_Percentile_Rank', 'Subcortical_Meta_Analysis_Rank.csv'))
si_full_rank <- read_csv(fullTableRank_Flip, here('Results', 'Tables', 'Genome_Percentile_Rank', 'Sex_interaction_Full_Meta_Analysis_Rank.csv'))
si_cortical_rank <- read_csv(corticalTableRank_Flip, here('Results', 'Tables', 'Genome_Percentile_Rank', 'Sex_interaction_Cortical_Meta_Analysis_Rank.csv'))
si_subcortical_rank <- read_csv(subcorticalTableRank_Flip, here('Results', 'Tables', 'Genome_Percentile_Rank', 'Sex_interaction_Subcortical_Meta_Analysis_Rank.csv'))


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


ma <- drive_get("~/Thesis/Manuscript/gs_tables/Slim_Tables/Official_Meta_Analysis")
if(nrow(ma) != 0) {
  drive_rm(ma)
}
#create the google worksheet
ma <- gs4_create("Official_Meta_Analysis", sheets = c('First_model_Full_Meta_Analysis', 'First_model_Male_Meta_Analysis','First_model_Female_Meta_Analysis','First_model_Cortical_Meta_Analysis','First_model_Subcortical_Meta_Analysis','Sex_interaction_Full_Meta_Analysis', 'Sex_interaction_Cortical_Meta_Analysis', 'Sex_interaction_Subcortical_Meta_Analysis'))
sheet_write(fullTable %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction,min_p_study, meta_p, signed_log_meta, Bonferroni_meta_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma, sheet = "First_model_Full_Meta_Analysis")
sheet_write(femaleTable %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction,min_p_study, meta_p, signed_log_meta, Bonferroni_meta_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma, sheet = 'First_model_Female_Meta_Analysis')
sheet_write(maleTable %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir,  meta_direction,min_p_study, meta_p, signed_log_meta, Bonferroni_meta_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma,'First_model_Male_Meta_Analysis')
sheet_write(corticalTable %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction,min_p_study, meta_p, signed_log_meta, Bonferroni_meta_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma, 'First_model_Cortical_Meta_Analysis')
sheet_write(subcorticalTable %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir,  meta_direction,min_p_study, meta_p, signed_log_meta, Bonferroni_meta_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma, 'First_model_Subcortical_Meta_Analysis')
sheet_write(fullTable_Flip  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir,  meta_direction,min_p_study, meta_p, signed_log_meta, Bonferroni_meta_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma, sheet = "Sex_interaction_Full_Meta_Analysis")
sheet_write(corticalTable_Flip  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir,  meta_direction,min_p_study, meta_p, signed_log_meta, Bonferroni_meta_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma, 'Sex_interaction_Cortical_Meta_Analysis')
sheet_write(subcorticalTable_Flip  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction,min_p_study, meta_p, signed_log_meta, Bonferroni_meta_p, slim_region_location, cell_type_taxon,cell_type_taxon_zscore, DE_Prior_Rank), ma, 'Sex_interaction_Subcortical_Meta_Analysis')

drive_mv(file = "Official_Meta_Analysis", path = "~/Thesis/Manuscript/gs_tables/Slim_Tables/")  # move Sheets file

mar <- drive_get("Thesis/Manuscript/gs_tables/Slim_Tables/Official_Genome_Percentile_Rank")
if(nrow(mar) != 0) {
  drive_rm(mar)
}
#create the google worksheet
mar <- gs4_create("Official_Genome_Percentile_Rank", sheets = c('First_model_Full_Meta_Analysis_Rank', 'First_model_Male_Meta_Analysis_Rank','First_model_Female_Meta_Analysis_Rank', 'First_model_Cortical_Meta_Analysis_Rank', 'First_model_Subcortical_Meta_Analysis_Rank','Sex_interaction_Full_Meta_Analysis_Rank', 'Sex_interaction_Cortical_Meta_Analysis_Rank', 'Sex_interaction_Subcortical_Meta_Analysis_Rank'))
sheet_write(fullTableRank  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), mar, sheet = "First_model_Full_Meta_Analysis_Rank")
sheet_write(femaleTableRank  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), mar, sheet = 'First_model_Female_Meta_Analysis_Rank')
sheet_write(maleTableRank  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), mar,'First_model_Male_Meta_Analysis_Rank')
sheet_write(corticalTableRank  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), mar, 'First_model_Cortical_Meta_Analysis_Rank')
sheet_write(subcorticalTableRank  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), mar, 'First_model_Subcortical_Meta_Analysis_Rank')
sheet_write(fullTableRank_Flip  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), mar, sheet = "Sex_interaction_Full_Meta_Analysis_Rank")
sheet_write(corticalTableRank_Flip  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), mar, 'Sex_interaction_Cortical_Meta_Analysis_Rank')
sheet_write(subcorticalTableRank_Flip  %>% dplyr::select(gene_symbol, Howard_pvalue, gene_name,RamakerDir, LabonteDir, DingDir, meta_direction, meta_empirical_p, Bonferroni_meta_empirical_p, slim_region_location, cell_type_taxon, DE_Prior_Rank), mar, 'Sex_interaction_Subcortical_Meta_Analysis_Rank')

drive_mv(file = "Official_Genome_Percentile_Rank", path = "~/Thesis/Manuscript/gs_tables/Slim_Tables/")  # move Sheets file

