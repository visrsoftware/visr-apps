################################################################
visr.category("Construct single-cell trajectories",
              info = "Discover cells transition from one state to another, in development, disease, and throughout life.",
              active.condition = "visr.param.output_dir != '' && visr.param.enable_trajectories")
################################################################

visr.param("trajectory_num_genes", "Number of top DE genes to use", default = 1000L,
           info = "Number of top significantly differentially expressed genes used as the ordering genes for the trajectory reconstruction.")

visr.param("trajectory_max_qval", "Max FDR threshold for selected DE genes", default = 0.1,
           items = "NULL", item.labels = "None",
           info = "Select differentially expressed genes that are significant at an FDR < specified threshold.")

################################################################
################################################################
################################################################

#'
#' Construct single cell trajectories
#'
perform_construct_trajectories <- function(monocle_app_object) {
  my_cds <- monocle_app_object$cds
  stopifnot(is(my_cds, "CellDataSet"))

  if (visr.var.add_title_pages)
    plotTitlePage("Constructing single cell trajectories", 2)

  # Step 1: determine ordering genes based on genes that differ between clusters
  de_genes <- monocle_app_object$de_genes
  if (!is.null(visr.param.trajectory_max_qval)) {
    trajectory_ordering_genes <- de_genes[which(de_genes$qval <= visr.param.trajectory_max_qval),]
  }
  else {
    trajectory_ordering_genes <- de_genes[which(is.na(de_genes$qval)),]
  }
  # pick top 1,000
  trajectory_ordering_genes <- trajectory_ordering_genes[1:min(visr.param.trajectory_num_genes, nrow(trajectory_ordering_genes)),]
  my_cds <- monocle::setOrderingFilter(my_cds, ordering_genes = trajectory_ordering_genes$id)

  p <- monocle::plot_ordering_genes(my_cds) +
    scale_x_log10(labels=axis_plain_format) +
    scale_y_log10(labels=axis_plain_format) +
    labs(x = "average gene expression (log scale)", y = "dispersion (log scale)") +
    ggtitle(sprintf("Step 1: selected %d genes for ordering", nrow(trajectory_ordering_genes)))
  print(p)

  # Step 2: reducing dimension using DDRTree (Reversed Graph Embedding)
  visr.logProgress(sprintf("Reducing the dimensionaly of %d ordering genes using 'DDRTree' method ...",
                           length(trajectory_ordering_genes$id)))
  my_cds <- monocle::reduceDimension(my_cds,
                                     max_components = 2,
                                     reduction_method = 'DDRTree', # c("DDRTree", "ICA", "tSNE", "SimplePPT", "L1-graph", "SGL-tree")
                                     norm_method = 'log',
                                     verbose = T)

  p <- ggplot(data.frame(t(monocle::reducedDimS(my_cds)))) +
    geom_point(aes(X1, X2)) +
    theme_bw() +
    xlab("Component1") + ylab("Component2") +
    ggtitle("Step 2: reducing dimension using DDRTree method")
  print(p)

  # Step 3: ordering the cells in pseudotime (fit the best tree it can to the data)
  visr.logProgress(sprintf("Ordering %d cells over the trajectory...", ncol(my_cds)))
  my_cds <- monocle::orderCells(my_cds)

  ## plot trajectories
  p <- monocle::plot_cell_trajectory(my_cds , color_by = "Pseudotime") +
    ggtitle("Step 3: ordering cells in pseudotime")
  print(p)

  p <- monocle::plot_cell_trajectory(my_cds , color_by = "State")
  print(p)

  if ("Cluster" %in% names(pData(my_cds))) {
    p <- monocle::plot_cell_trajectory(my_cds , color_by = "Cluster")
    print(p)
  }

  monocle_app_object$cds <- my_cds
  return(monocle_app_object)
}

