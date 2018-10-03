################################################################
visr.category("Find pseudotime changing genes",
              info = "Find genes that change as a function of pseudotime.",
              active.condition = "visr.param.output_dir != '' && visr.param.enable_trajectories && visr.param.find_pseudotime_genes")
################################################################

visr.param("cluster_genes_by_pseudo", label = "Cluster genes by pseudotime", default = T,
           info = "Hierarchical clustering of genes by pseudotemporal expression pattern")

visr.param("cluster_genes_pseudo_count", label = "Number of genes to cluster",
           default = 50L, min = 2L,
           info = "Number of top genes varying as a function of pseudotime to be used for clustering",
           active.condition = "visr.param.cluster_genes_by_pseudo")

visr.param("num_pseudo_gene_clusters", label = "Number of heatmap clusters", default = 3L, min = 1L,
           info = "Number of clusters for the heatmap of branch genes",
           active.condition = "visr.param.cluster_genes_by_pseudo")
