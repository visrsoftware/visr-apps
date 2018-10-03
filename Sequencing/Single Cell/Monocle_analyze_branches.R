################################################################
visr.category("Analyze branches in trajectories",
              info = "Analyze branches in single-cell trajectories to identify the genes that differ at a particular branch point.",
              active.condition = "visr.param.output_dir != '' && visr.param.enable_trajectories && visr.param.analyze_trajectory_branches")
################################################################

visr.param("trajectory_branch_point", label = "Branch point number", default = 1L,
           info = "The ID of the branch point to analyze")

visr.param("branched_heatmap_num_clusters", label = "Number of clusters for the heatmap",
           default = 6L, min = 2L,
           info = "Number of clusters for the heatmap of branch genes")

visr.param("num_branch_genes_to_plot", label = "Number of branch dependent genes to plot",
           default = 6L,
           "Number of branch dependent genes to plot per cluster")
