##### This folder contains the files used to generate the heatmaps

- The main file that generates the combined heatmaps is **combine_heatmaps.R** which is called in the main script **Run_meta_analysis_script.R**
- This script can only be successfully completed after all meta-analyses have been performed and their results saved in the specified folders
- The main file makes calls to:
  - **combine_transcriptomic_heatmaps.R**
  	- This script combines the meta analysis results from each transcriptomic study generated in:
  		- **Labonte_heatmap.R**
  			- Heatmap is generated using function definition in **heatmaps.R**
  		- **Ramaker_heatmap.R**
  			- Heatmap is generated using function definition in **heatmaps.R**
  		- **Ding_heatmap.R**
  			- Heatmap is generated using function definition in **heatmaps.R**
  - **meta_p_heatmaps.R**
  	- This script creates a heatmap to display the corrected meta p-values of our top 12 genes and how they compared in each meta-analysis performed