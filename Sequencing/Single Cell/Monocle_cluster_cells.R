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
