
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
