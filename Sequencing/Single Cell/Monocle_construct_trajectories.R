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
