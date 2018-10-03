################################################################
visr.category("Find pseudotime changing genes",
              info = "Find genes that change as a function of pseudotime.",
              active.condition = "visr.param.output_dir != '' && visr.param.find_pseudotime_genes")
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

################################################################
################################################################
################################################################

#'
#' Finding genes that change as a function of pseudotime
#'
perform_pseudotime_genes <- function(monocle_app_object) {
  my_cds <- monocle_app_object$cds
  stopifnot(is(my_cds, "CellDataSet"))

  if (visr.var.add_title_pages)
    plotTitlePage("DE genes as a function of pseudotime", 2)

  visr.logProgress("Finding genes that change as a function of pseudotime ...")
  my_pseudotime_de <- monocle::differentialGeneTest(
    my_cds, fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = visr.param.num_cores)

  my_pseudotime_de %>% arrange(qval) %>% head()

  my_pseudotime_de %>% arrange(qval) %>% head() %>% select(id) -> my_pseudotime_gene
  my_pseudotime_gene <- my_pseudotime_gene$id

  p <- monocle::plot_genes_in_pseudotime(my_cds[my_pseudotime_gene,])
  print(p)
  # todo: parameters

  if (visr.param.cluster_genes_by_pseudo) {
    visr.logProgress("Clustering genes by pseudotemporal expression pattern")

    # cluster the top genes that vary as a function of pseudotime
    my_pseudotime_de %>% arrange(qval) %>% head(visr.param.cluster_genes_pseudo_count) %>% select(id) -> gene_to_cluster
    gene_to_cluster <- gene_to_cluster$id

    plot.new()
    my_pseudotime_cluster <- monocle::plot_pseudotime_heatmap(
      my_cds[gene_to_cluster,],
      num_clusters = visr.param.num_pseudo_gene_clusters,
      cores = visr.param.num_cores,
      # norm_method = c("log", "vstExprs") # Determines how to transform expression values prior to rendering
      # scale_max = 3, # The maximum value (in standard deviations) to show in the heatmap. Values larger than this are set to the max.
      # scale_min = -3 # The minimum value (in standard deviations) to show in the heatmap. Values smaller than this are set to the min.
      show_rownames = TRUE,
      return_heatmap = TRUE
    )
    mtext(text = "Modules of genes that co-vary across pseudotime",
          adj=0, side=1, line = 4, cex = 1)

    # todo: add the cluster ids to genes.txt
    visr.logProgress("Extracting the genes for each cluster ...")
    my_cluster <- cutree(my_pseudotime_cluster$tree_row, visr.param.num_pseudo_gene_clusters)
    for (cluster_id in seq(visr.param.num_pseudo_gene_clusters)) {
      print(sprintf("Genes for cluster %d", cluster_id))
      print(my_pseudotime_de[names(my_cluster[my_cluster == cluster_id]), "gene_short_name"])
    }
  }

  monocle_app_object$cds <- my_cds
  return(monocle_app_object)
}
