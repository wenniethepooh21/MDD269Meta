library(png)
library(grid)
library(gridExtra)
library(here)
library(cowplot)
library(ggplot2)


Labonte_heat <- readPNG(here('Processed_Data/Meta_Analysis_Results/Heatmaps/top_genes_Labonte_expression_heatmap.png'))
Labonte_heat %<>% rasterGrob(interpolate = FALSE)
Ramaker_heat <- readPNG(here('Processed_Data/Meta_Analysis_Results/Heatmaps/top_genes_Ramaker_expression_heatmap.png'))
Ramaker_heat %<>% rasterGrob(interpolate = FALSE)
Ding_heat <- readPNG(here('Processed_Data/Meta_Analysis_Results/Heatmaps/top_genes_Ding_expression_heatmap.png'))
Ding_heat %<>% rasterGrob(interpolate = FALSE)
plots <- 
heat_plots <- list(Labonte_heat,Ramaker_heat,Ding_heat)
#Render a raster object (bitmap image) at the given location, size, and orientation.
heat_plots_raster <- lapply(heat_plots, rasterGrob)
combined_transcriptomic<- do.call(grid.arrange, c(heat_plots_raster,nrow = 1))
g <- arrangeGrob(heat_plots_raster, nrow=1)
ggsave(filename = here('Processed_Data/Meta_Analysis_Results/Heatmaps/combined_transcriptomic_heatmaps.png'),g, dpi=300, width=11, height=8)


plots <- lapply(ll <- list.files(heat_plots),function(x){
                  heatmap_plot <- as.raster(x)
  rasterGrob(heatmap_plot, interpolate = FALSE)
})

ggsave(here('Processed_Data/Meta_Analysis_Results/Heatmaps/combined_transcriptomic_heatmaps.png'), marrangeGrob(grobs=plots, nrow=1, ncol=3))
