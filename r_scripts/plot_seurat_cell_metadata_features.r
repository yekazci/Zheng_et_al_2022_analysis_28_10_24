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

