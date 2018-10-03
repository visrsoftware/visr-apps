################################################################
visr.category("Analyze branches in trajectories",
              info = "Analyze branches in single-cell trajectories to identify the genes that differ at a particular branch point.",
              active.condition = "visr.param.output_dir != '' && visr.param.analyze_trajectory_branches")
################################################################

visr.param("trajectory_branch_point", label = "Branch point number", default = 1L,
           info = "The ID of the branch point to analyze")

visr.param("branched_heatmap_num_clusters", label = "Number of clusters for the heatmap",
           default = 6L, min = 2L,
           info = "Number of clusters for the heatmap of branch genes")

visr.param("num_branch_genes_to_plot", label = "Number of branch dependent genes to plot",
           default = 6L,
           "Number of branch dependent genes to plot per cluster")


################################################################
################################################################
################################################################

#'
#' Analyze branches in single-cell trajectories
#'
perform_analyze_branches <- function(monocle_app_object) {
  my_cds <- monocle_app_object$cds
  stopifnot(is(my_cds, "CellDataSet"))

  if (visr.var.add_title_pages)
    plotTitlePage("Analyzing branches in single-cell trajectories", 2)

  trajectory_branch_point <- visr.param.trajectory_branch_point
  # A table of genes is returned with significance values that indicate whether genes have expression patterns that are branch dependent.
  BEAM_res <- monocle::BEAM(my_cds,
                            branch_point = trajectory_branch_point,
                            cores = visr.param.num_cores)
  BEAM_res <- BEAM_res[order(BEAM_res$qval),]
  BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
  # check out the results
  print(head(BEAM_res))

  plot.new()
  my_branched_heatmap <- monocle::plot_genes_branched_heatmap(my_cds[row.names(subset(BEAM_res, qval < 1e-4)),],
                                                              branch_point = trajectory_branch_point,
                                                              num_clusters = visr.param.branched_heatmap_num_clusters,
                                                              cores = visr.param.num_cores,
                                                              # scale_max = 3,
                                                              # scale_min = -3,
                                                              # norm_method = c("log", "vstExprs")
                                                              use_gene_short_name = TRUE,
                                                              show_rownames = TRUE,
                                                              return_heatmap = TRUE)
  mtext(text = sprintf("Bifurcation of gene expressions along two lineages of branch point %d", trajectory_branch_point),
        adj=0, side=1, line = 4, cex = 1)

  print(table(my_branched_heatmap$annotation_row$Cluster))

  my_row <- my_branched_heatmap$annotation_row
  my_row <- data.frame(cluster = my_row$Cluster,
                       gene = row.names(my_row),
                       stringsAsFactors = FALSE)

  for (cluster_id in seq(visr.param.branched_heatmap_num_clusters)) {
    this_plot_title <- sprintf("Top genes of cluster #%d, expressed in a branch dependent manner at branch #%d",
                               cluster_id, trajectory_branch_point)
    visr.logProgress(this_plot_title)
    my_gene <- row.names(subset(fData(my_cds),
                                gene_short_name %in% head(my_row[my_row$cluster == cluster_id, 'gene'],
                                                          visr.param.num_branch_genes_to_plot)))

    # plot genes that are expressed in a branch dependent manner
    p <- monocle::plot_genes_branched_pseudotime(
      my_cds[my_gene,], branch_point = trajectory_branch_point, ncol = 1) +
      ggtitle(this_plot_title) +
      theme(plot.title = element_text(size = 10))
    print(p)
  }

  monocle_app_object$cds <- my_cds
  return(monocle_app_object)
}
