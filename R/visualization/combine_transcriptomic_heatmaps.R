library(here)
library(cowplot)
library(ggplot2)


source(here("R/visualization/Labonte_heatmap.R")) #generates a heatmap of each top gene expression direction across each sampled brain region in Labonte data
source(here("R/visualization/Ramaker_heatmap.R")) #generates a heatmap of each top gene expression direction across each sampled brain region in Ramaker data
source(here("R/visualization/Ding_heatmap.R")) #generates a heatmap of each top gene expression direction across each sampled brain region in Ding data

getExpressionHeat <- function() {
  labonte <- drawLabonte()
  ramaker<-drawRamaker()
  ding <-drawDing()
  
  # # add title to the combined plots 
  # title <- ggdraw() + draw_label("Direction of Gene Expressions", fontface = "bold",size = 20, vjust = 1)
  # 
  
  combined <- plot_grid(
            ramaker + theme(plot.margin = unit(c(t=0, r=0.5, b=0, l=0), "cm")),
             labonte + theme(plot.margin = unit(c(t=0, r=0.5, b=0, l=0), "cm")),
             ding + theme(plot.margin = unit(c(t=0, r=0, b=0, l=0), "cm")),
             nrow = 1, 
            align = "h", 
            rel_widths = c(1, 1.3))
  # all_heats<-plot_grid(title, combined, ncol = 1, rel_heights = c(0.1,3))
  
  ggsave(filename = here('Processed_Data/Meta_Analysis_Results/Heatmaps/combined_transcriptomic_heatmaps.png'), dpi=300, width=15.5, height=8)
  
  return(combined)
}