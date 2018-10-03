################################################################
visr.category("Clustering",
              info = "Identify new (and possibly rare) subtypes of cells in single-cell RNA-seq experiments.",
              active.condition = "visr.param.output_dir != '' && visr.param.enable_clustering")
################################################################

CLUSTER_METHOD_DENSITY_PEAK = "densityPeak"
CLUSTER_METHOD_LOUVAIN = "louvain"
CLUSTER_METHOD_DDRTREE = "DDRTree"
visr.param("cluster_method", items = c(CLUSTER_METHOD_DENSITY_PEAK, CLUSTER_METHOD_LOUVAIN, CLUSTER_METHOD_DDRTREE),
           info = "Method for clustering cells. For big datasets (like data with 50 k cells or so), we recommend using the 'louvain' clustering algorithm.")

visr.param("num_clusters", label = "Number of clusters", type = "integer", min = 1L, default = "NULL", items = c("NULL"), item.labels = c("auto"), debugvalue = NULL,
           info = "Number of clusters. When auto, Use top 95% of the delta and top 95% of the rho as the cutoff for assigning density peaks and clusters",
           active.condition = sprintf("visr.param.cluster_method == '%s'", CLUSTER_METHOD_DENSITY_PEAK))

#visr.param("skip_rho_sigma", default = FALSE, info = "skip the calculation of rho / sigma",
#           active.condition = sprintf("visr.param.cluster_method == '%s'", CLUSTER_METHOD_DENSITY_PEAK))

visr.param("rho_threshold", label = "rho (cell's local density) threshod", type = "double", min = 0,
           items = "NULL", item.labels = "auto (95%)", default = "NULL", debugvalue = NULL,
           info = "The threshold of cell's local density (rho) used to select the density peaks",
           active.condition = sprintf("visr.param.cluster_method == '%s'", CLUSTER_METHOD_DENSITY_PEAK))

visr.param("delta_threshold", label = "delta (local distance) threshod", type = "double", min = 0,
           items = "NULL", item.labels = "auto (95%)", default = "NULL", debugvalue = NULL,
           info = "The threshold of local distance (nearest distance of a cell to another cell with higher distance) used to select the density peaks",
           active.condition = sprintf("visr.param.cluster_method == '%s'", CLUSTER_METHOD_DENSITY_PEAK))

visr.param("gaussian", label = "Use gaussian kernel?", default = T,
           info = "Whether or not Gaussian kernel will be used for calculating the local density in densityClust function",
           active.condition = sprintf("visr.param.cluster_method == '%s'", CLUSTER_METHOD_DENSITY_PEAK))

visr.param("num_centers", label = "Number of centroids", type = "integer", min = 1L, default = 3L,
           info = "Number of number of centroids passed to DDRTree ('Dimensionality Reduction via Graph Structure Learning' by Qi Mao, et al)",
           active.condition = sprintf("visr.param.cluster_method == '%s'", CLUSTER_METHOD_DDRTREE))

visr.param("louvain_k", default = 50L,
           info = "number of kNN used in creating the k nearest neighbor graph for Louvain clustering. The number of kNN is related to the resolution of the clustering result, bigger number of kNN gives low resolution and vice versa.",
           active.condition = sprintf("visr.param.cluster_method == '%s'", CLUSTER_METHOD_LOUVAIN))

visr.param("louvain_iter", default = 1L,
           info = "number of iterations used for Louvain clustering. The clustering result gives the largest modularity score will be used as the final clustering result.",
           active.condition = sprintf("visr.param.cluster_method == '%s'", CLUSTER_METHOD_LOUVAIN))

visr.param("louvain_weight", default = F,
           info = "Use Jaccard coefficent for two nearest neighbors (based on the overlapping of their kNN) as the weight used for Louvain clustering.",
           active.condition = sprintf("visr.param.cluster_method == '%s'", CLUSTER_METHOD_LOUVAIN))

################################################################
################################################################
################################################################

#'
#' performs cell clustering after dimensionality reduction
#'
perform_clustering <- function(monocle_app_object) {
  my_cds <- monocle_app_object$cds
  stopifnot(is(my_cds, "CellDataSet"))

  if (visr.var.add_title_pages)
    plotTitlePage("Clustering")

  num_clusters = NULL
  if (visr.param.cluster_method == CLUSTER_METHOD_DENSITY_PEAK) {
    # for some reason it alwasy produces one less cluster than specified value. so add 1 for now
    num_clusters = if (!is.null(visr.param.num_clusters)) (visr.param.num_clusters + 1) else NULL
  } else if (visr.param.cluster_method == CLUSTER_METHOD_DDRTREE) {
    num_clusters = visr.param.num_centers
  }

  visr.logProgress(paste("Performing unsupervised clustering. method:", visr.param.cluster_method))
  my_cds <- monocle::clusterCells(my_cds,
                                  verbose = T,
                                  skip_rho_sigma = F, #visr.param.skip_rho_sigma, # whether you want to skip the calculation of rho / sigma
                                  num_clusters = num_clusters,
                                  inspect_rho_sigma = F, # whether you want to interactively select the rho and sigma
                                  rho_threshold = visr.param.rho_threshold,
                                  delta_threshold = visr.param.delta_threshold,
                                  gaussian = visr.param.gaussian, # whether use Gaussian kernel for calculating the local density
                                  peaks = NULL, # numeric vector indicating the index of density peaks used for clustering.
                                  cell_type_hierarchy = NULL,
                                  frequency_thresh = NULL,
                                  enrichment_thresh = NULL,
                                  clustering_genes = NULL,
                                  k = visr.param.louvain_k,
                                  louvain_iter = visr.param.louvain_iter,
                                  weight = visr.param.louvain_weight,
                                  method = visr.param.cluster_method)

  head(pData(my_cds))
  print("Cluster sizes:")
  print(table(pData(my_cds)$Cluster))
  if (!is.null(visr.param.rho_threshold) & !is.null(visr.param.delta_threshold)) {
    title_peak <- sprintf("Detected %d peaks using rho=%4.2f and delta = %4.2f ",
                          length(which(pData(my_cds)$peaks)),
                          visr.param.rho_threshold, visr.param.delta_threshold)
  } else {
    title_peak <- sprintf("Detected %d peaks using top 95%% delta and top 95%% rho for threshold.",
                          length(which(pData(my_cds)$peaks)))
  }

  if (visr.param.cluster_method == CLUSTER_METHOD_DENSITY_PEAK) {
    visr.logProgress("Plotting the decision map of density clusters (delta vs. rho)")
    p <- monocle::plot_rho_delta(my_cds, rho_threshold = visr.param.rho_threshold, delta_threshold = visr.param.delta_threshold) +
      ggtitle(paste0("Decision map of density clusters\nPeaks are cells with high local density that are far away from other cells with high local density\n", title_peak)) +
      theme(plot.title = element_text(size = 12)) + # hjust = 0.5
      labs(x = "rho: local density", y = "delta: local distance (to another cell with higher density)")
    print(p)
  }

  p <- monocle::plot_cell_clusters(my_cds) +
    xlab("tSNE1") + ylab("tSNE2") +
    ggtitle(sprintf("Unsupervised clustering of cells\nusing %s method", visr.param.cluster_method))
  print(p)

  monocle_app_object$cds <- my_cds
  return(monocle_app_object)
}
