library(cowplot)
library(ggplot2)
### This script draws the heatmaps showing the expression levels of each top gene across all brain regions for both sexes

# This function draws the heatmap of both sexes and combines them into one plot
# arguments: male and female expression data (dataframe), brain regions used in each sex (list), gene_symbols used in each sex (list)
# returns the combined plot
drawExpressionHeat <- function(dataset_male, dataset_female, male_regions, male_symbol, female_regions, female_symbol, male_direction, female_direction, legend_title) {
  #heatmap
  male_heatmap <- ggplot(dataset_male, aes(x=male_regions,y=male_symbol)) + 
    geom_tile(aes(fill = male_direction), colour='black') +
    labs(x=NULL) +
    scale_fill_distiller(palette = "RdBu") +
    ylab('') +
    ggtitle("Male expression levels") +
    theme(
          axis.text.x = element_text(size = 12, angle = 45, hjust=1,vjust=1),
          axis.text.y = element_text(face = "italic"),
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(t=0, r=-0.5, b=0, l=0), "cm"),
          panel.background = element_blank()) 
  #heatmap
  female_heatmap <- ggplot(dataset_female, aes(x=female_regions,y=female_symbol)) + 
    geom_tile(aes(fill = female_direction), colour='black') +
    labs(x=NULL, fill = "Expression\nsigned log(p_value)") +
    scale_fill_distiller(palette = "RdBu") +
    ylab('') +
    ggtitle("Female expression levels")+
    theme(
          axis.text.x = element_text(size = 12, angle = 45, hjust=1,vjust=1),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.margin = unit(c(t=0, r=0, b=0, l=-0.5), "cm"),
          panel.background = element_blank())
  
  combined_plots <- plot_grid(male_heatmap + theme(legend.position = "none"), 
                              female_heatmap + theme(legend.position = "none"),  
                              align = 'vh'
                              )
  
  # extract the legend from one of the plots
  legend <- get_legend(
    # create some space to the left of the legend
    female_heatmap + theme(legend.box.margin = margin(0, 0, 0))
  )
  # add the legend to the space just made. Give it one-third of 
  # the width of one plot (via rel_widths).
  combined_plots_legend <- plot_grid(combined_plots, legend, rel_widths = c(2,0.3))
  # print(combined_plots_legend)
  return(combined_plots_legend)
}