library(here)
library(cowplot)
library(ggplot2)

source(here("R/visualization/combine_transcriptomic_heatmaps.R")) 
source(here("R/visualization/meta_p_heatmap.R")) 


expression_heats <- getExpressionHeat() 
meta_heatmap <- drawMetaHeat()


combined <- plot_grid(
  expression_heats + theme(plot.margin = unit(c(t=0, r=0.5, b=0, l=0), "cm")),
  meta_plot + theme(plot.margin = unit(c(t=0, r=0, b=0, l=0), "cm")),
  nrow = 1,
  align = "h",
  rel_heights = c(0.1,3))
  # rel_widths = c(.5, 1.3))
print(combined)
ggsave(filename = here('Processed_Data/Meta_Analysis_Results/Heatmaps/combined_transcriptomic_heatmaps.png'), dpi=300, width=15.5, height=8)
