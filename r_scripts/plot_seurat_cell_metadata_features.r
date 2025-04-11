# Now I will write a function to view all plots pairwise with desired discrete variable plot:

plot_discrete_vars_w_ref <- function(seurat_object, discrete_vars, ref_feature, my_reduction, my_colours, my_margin) { 


discrete_vars <- setdiff(discrete_vars, ref_feature)
    
plot_list_discrete_vars <- list()

 plot_REF <- 
        DimPlot(object = seurat_object,
            shuffle = TRUE,
            cols = my_colours,
            reduction = my_reduction, 
            group.by =  ref_feature) + my_margin
    
for (each_feature in discrete_vars) {

 plot_list_discrete_vars[[each_feature]] <- 
        DimPlot(object = seurat_object,
            shuffle = TRUE,
            cols = my_colours,
            reduction = my_reduction, 
            group.by =  each_feature) + my_margin

 plot_list_discrete_vars[[paste0('ref_', ref_feature, '_plot_for_', each_feature)]] <- plot_REF

}
     
    return(plot_list_discrete_vars)
    
}

#################################################

# Now I will write a function to view all continuous variable plots pairwise 
# with desired discrete variable plot:

plot_continuous_vars_w_ref <- function(seurat_object, continuous_vars, ref_feature, my_reduction, my_colours, my_margin) { 


continuous_vars <- setdiff(continuous_vars, ref_feature)
    
plot_list_continuous_vars <- list()

 plot_REF <- 
        DimPlot(object = seurat_object,
            shuffle = TRUE,
            cols = my_colours,
            reduction = my_reduction, 
            group.by =  ref_feature) + my_margin
    
for (each_feature in continuous_vars) {

 plot_list_continuous_vars[[each_feature]] <- 
        FeaturePlot(object = seurat_object,
            cols = c('lightgray', 'red'), # set two colours for gradients. 
            reduction = my_reduction, 
            features =  each_feature) + my_margin

 plot_list_continuous_vars[[paste0('ref_', ref_feature, '_plot_for_', each_feature)]] <- plot_REF

}
     
    return(plot_list_continuous_vars)
    
}

################################################

# Now I will write a function generate continuous variables plots:

plot_continuous_vars <- function(seurat_object, continuous_vars, my_reduction, my_colours, my_margin) { 

plot_list_continuous_vars <- list()

for (each_feature in continuous_vars) {

 plot_list_continuous_vars[[each_feature]] <- 
        FeaturePlot(object = seurat_object,
            cols = c('lightgray', 'red'), # set two colours for gradients. 
            reduction = my_reduction, 
            features =  each_feature) + my_margin
}

    return(plot_list_continuous_vars)
    
}

###############################################

# write a function to generate discrete variable plot:

plot_discrete_vars <- function(seurat_object, discrete_vars, my_reduction, my_colours, my_margin) { 

plot_list_discrete_vars <- list()

    
for (each_feature in discrete_vars) {

 plot_list_discrete_vars[[each_feature]] <- 
        DimPlot(object = seurat_object,
            shuffle = TRUE,
            cols = my_colours,
            reduction = my_reduction, 
            group.by =  each_feature) + my_margin
}

    return(plot_list_discrete_vars)
    
}

##################################################

plot_read_count_knee_from_BPCells <- function(read_counts, cutoff = NULL, return_data = FALSE, apply_styling = TRUE) {
  
  # Message to inform the user about how the cutoff should be provided
  message('if used cutoff, make sure to provide it in e3 format (1e3 etc.) ')
  
  # Step 1: Prepare ranks for plotting. We want to keep ~1k cells per order of magnitude
  # We calculate unique ranks using a logarithmic scale. This helps in speeding up plotting by reducing the number of points.
  # floor(10^(1:(1000 * ceiling(log10(length(read_counts)))) / 1000)) will generate a sequence of ranks.
  ranks <- unique(floor(10^(1:(1000 * ceiling(log10(length(read_counts)))) / 1000)))
  
  # Step 2: Prepare the data frame for plotting
  data <- list(
    # Create a tibble that contains 'ranks' and 'reads'
    data = tibble::tibble(
      ranks = ranks[ranks < length(read_counts)],  # Ensure ranks do not exceed the number of read counts
      reads = sort(read_counts, decreasing = TRUE)[ranks]  # Sort read counts in descending order and select by ranks
    ),
    cells_passing = NA  # Initialize 'cells_passing' as NA (we'll compute this if cutoff is provided)
  )
  
  # Step 3: If cutoff is provided, calculate the number of cells passing the cutoff
  # This step counts how many cells have a read count greater than the specified cutoff
  if (!is.null(cutoff)) {
    data$cells_passing <- sum(read_counts > cutoff)
  }
  
  # Step 4: If return_data is TRUE, return the prepared data (without plotting)
  if (return_data) {
    return(data)
  }

  # Step 5: Plotting the data using ggplot2
  # ggplot creates a scatter plot where the x-axis is the log10 of the ranks and y-axis is the log10 of the reads
  plot <- ggplot2::ggplot(data$data, ggplot2::aes(log10(ranks), log10(reads))) +
    ggplot2::geom_point() +  # Add points to the plot

    # Formatting the axes: Set breaks on the x and y axes for clarity and use log10 scaling
    ggplot2::scale_x_continuous(labels = scales::label_math(), breaks = scales::breaks_width(1)) +
    ggplot2::scale_y_continuous(labels = scales::label_math(), breaks = scales::breaks_width(1))

  # Step 6: Label the axes with appropriate descriptions
  plot <- plot + ggplot2::labs(x = "Barcode Rank", y = "Reads")

  # Step 7: If apply_styling is FALSE, return the plot without applying further styling
  if (!apply_styling) {
    return(plot)
  }

  # Step 8: If cutoff is provided, enhance the plot with additional elements:
  #  - Add a shaded rectangle indicating cells that pass the cutoff
  #  - Add a dashed horizontal line at the cutoff level
  #  - Add a label displaying the number of cells passing the cutoff
  
  if (!is.null(cutoff)) {
    # Create a label that will display the number of cells passing the cutoff
    cell_label <- tibble::tibble(
      label = sprintf("%s cells", scales::label_comma()(data$cells_passing)),  # Format label for cell count
      x = max(log10(data$data$ranks)),  # Position the label at the max rank value
      y = max(log10(data$data$reads))  # Position the label at the max read count value
    )
    
    # Define the rectangle that will highlight the area above the cutoff threshold
    rectangle_highlight <- tibble::tibble(
      xmin = -Inf, xmax = Inf,  # Set the rectangle to span the entire width of the plot
      ymin = log10(cutoff), ymax = Inf  # The rectangle will start at the log10 of the cutoff value and extend upwards
    )
    
    # Add the cutoff line, label, and rectangle to the plot
    plot <- plot +
      ggplot2::geom_hline(yintercept = log10(cutoff), linetype = "dashed") +  # Dashed line at the cutoff level
      ggplot2::geom_text(data = cell_label, ggplot2::aes(x, y, label = label), hjust = "inward", vjust = "inward") +  # Add label for the number of passing cells
      ggplot2::geom_rect(
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, x = NULL, y = NULL),  # Draw the shaded rectangle
        data = rectangle_highlight,
        alpha = 0.1  # Set transparency of the shaded rectangle
      )
  }

  # Step 9: Final plot styling with log ticks and classic theme
  plot <- plot +
    ggplot2::annotation_logticks() +  # Add log scale ticks to both axes
    ggplot2::theme_classic()  # Use the classic ggplot2 theme for the plot
  
  # Step 10: Return the final plot
  return(plot)
}

