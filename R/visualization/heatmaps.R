library(cowplot)
library(ggplot2)
library(gridExtra)
### This script draws the heatmaps showing the expression levels of each top gene across all brain regions for both sexes

# This function draws the heatmap of both sexes and combines them into one plot
# arguments: male and female expression data (dataframe), brain regions used in each sex (list), gene_symbols used in each sex (list)
# returns the combined plot

setBins <- function(dataframe) {
  dataframe$bin <- cut(dataframe$Bonferroni_Correction, breaks = c())
}
drawExpressionHeat <- function(dataset_male, dataset_female, male_regions, male_symbol, female_regions, female_symbol, male_direction, female_direction, study = "other") {
  #only draw the gene symbols on y-axis for the Ramaker study
  if(study == "other") {
    y_axis = element_blank()
    r_widths = c(.5,.5)
  } else {
    y_axis = element_text(face = "italic")
    r_widths = c(1,.8)
  }
  

  #determine limits
  low_bound <- min(male_direction, female_direction, na.rm = TRUE)
  upper_bound <- max(male_direction, female_direction, na.rm = TRUE)
  abs_bound <- ceiling(max(abs(low_bound), abs(upper_bound)))
  print(abs_bound)

  #heatmap
  male_heatmap <- ggplot(dataset_male, aes(x=male_regions,y=male_symbol)) + 
    geom_tile(aes(fill = male_direction), colour='black') +
    labs(x=NULL) +
    scale_fill_distiller(limits = c(abs_bound*-1, abs_bound),
                                    palette = "RdBu") +
    # scale_fill_gradientn(limits = c(abs_bound*-1, abs_bound),
    # colors = c("darkblue", "blue", "white", "red", "darkred")) +
    ylab('') +
    ggtitle("Male") +
    theme(
          axis.text.x = element_text(size = 10, angle = 20, hjust=0.5,vjust=1),
          axis.ticks.y = element_blank(),
          axis.text.y = y_axis,
          plot.title = element_text(hjust = 0.5, size = 9),
          plot.margin = unit(c(t=0, r=0, b=0, l=-0.3), "cm"),
          panel.background = element_blank()) 
  #heatmap
  female_heatmap <- ggplot(dataset_female, aes(x=female_regions,y=female_symbol)) + 
    geom_tile(aes(fill = female_direction), colour='black') +
    labs(x=NULL, fill = "Expression\nsigned log(p_value)") +
    scale_fill_distiller(limits = c(abs_bound*-1, abs_bound),
                         palette = "RdBu") +
    # scale_fill_gradientn(limits = c(abs_bound*-1, abs_bound),
    # colors = c("darkblue", "blue", "white", "red", "darkred")) +
    ylab('') +
    ggtitle("Female")+
    theme(
          axis.text.x = element_text(size = 10, angle = 20, hjust=0.5,vjust=1),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 9),
          plot.margin = unit(c(t=0, r=0, b=0, l=-0.1), "cm"),
          panel.background = element_blank())


  combined_plots <- plot_grid(male_heatmap + theme(legend.position = "none"), 
                              female_heatmap + theme(legend.position = "none"),
                              nrow = 1,
                              rel_widths = r_widths)
  
  # extract a legend that is laid out horizontally
  legend <- get_legend(
    female_heatmap + 
      guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "top",
            legend.box.margin = unit(c(t=0, r=0, b=0.5, l=0), "cm"))
  )
  
  # add the legend underneath the row we made earlier. Give it 10%
  # of the height of one plot (via rel_heights).
  
  combined_plots_legend <- plot_grid(combined_plots, legend, ncol = 1, rel_heights = c(1, .1))

  return(combined_plots_legend)
}