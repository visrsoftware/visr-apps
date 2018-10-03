
################################################################
visr.category("Dimensionality reduction",
              info = "Values used to reduce dimensionality of data (Genes)",
              active.condition = "visr.param.output_dir != '' && visr.param.enable_dim_red")
################################################################

visr.param("reduce_num_dim", label = "Number of top principal components to use", default = 50L, min = 2L)

visr.param("reduction_method", label = "Algorithm for dimensionality reduction",
           items = c("DDRTree", "ICA", "tSNE", "SimplePPT", "L1-graph", "SGL-tree"),
           default = 'tSNE',
           info = "Algorithm to use for dimensionality reduction")

################################################################
################################################################
################################################################

#'
#' Performs PCA and dimensionality reduction
#'
perform_dim_reduction <- function(monocle_app_object) {
  my_cds <- monocle_app_object$cds
  stopifnot(is(my_cds, "CellDataSet"))

  if (visr.var.add_title_pages)
    plotTitlePage("Dimensionality Reduction")

  visr.logProgress("Plotting the percentage of variance explained by each component (takes a few minutes)")
  p <-  monocle::plot_pc_variance_explained(my_cds, return_all = FALSE) + #TODO: norm_method= c("log", "vstExprs", "none")
    ggtitle("Percentage of variance explained by each component\nbased on a PCA performed on the normalised expression data") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept = visr.param.reduce_num_dim, linetype = "dotted", color = 'red') +
    labs(x = "PCA Components", y = "Variance explained by each PCA component")
  print(p)

  visr.logProgress("Computing a projection of a CellDataSet object into a lower dimensional space (takes a few minutes)")
  my_cds <- monocle::reduceDimension(my_cds, max_components = 2, num_dim = visr.param.reduce_num_dim,
                                     reduction_method = visr.param.reduction_method, verbose = TRUE) # norm_method = c("log", "vstExprs", "none")

  p <- ggplot(data.frame(t(monocle::reducedDimA(my_cds)))) +
    geom_point(aes(X1, X2)) +
    theme_bw() +
    xlab("Dimension1") + ylab("Dimension2") +
    ggtitle(sprintf("Dimensionality reduction using %s method", visr.param.reduction_method))
  print(p)

  monocle_app_object$cds <- my_cds
  return(monocle_app_object)
}
