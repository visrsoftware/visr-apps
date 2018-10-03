################################################################
visr.category("Analysis steps",
              info = "Different analysis steps",
              active.condition = "visr.param.output_dir != '' && ((visr.param.data_dir_10x != '') || (visr.param.path_to_monocle_object != '') || (visr.param.expression_matrix != '' && visr.param.sample_sheet != '' && visr.param.gene_annotation != ''))")
################################################################

condition_monocle_object <- sprintf("((visr.param.input_type == '%s') && (visr.param.path_to_monocle_object != ''))", INPUT_TYPE_MONOCLE_OBJECT)

visr.param("enable_filtering", label = "Filtering cells and subsetting genes", default = T,
           info = "Remove outlier cells and genes before further processing.",
           debugvalue = F)

visr.param("enable_dim_red", label = "Dimensionality reduction", default = T,
           info = "Reduce dimensionality of data from many genes to fewer number of components.",
           debugvalue = F)

visr.param("enable_clustering", label = "Clustering", default = T,
           active.condition = sprintf("%s || (visr.param.enable_filtering && visr.param.enable_dim_red)", condition_monocle_object),
           info = "Identify subtypes of cells using unsupervised clustering.",
           debugvalue = F)

visr.param("enable_de_analysis", label = "Differential expression analysis", default = T,
           active.condition = sprintf("%s || (visr.param.enable_filtering && visr.param.enable_dim_red && visr.param.enable_clustering)", condition_monocle_object),
           info = "Characterize differentially expressed genes by comparing groups of cells.",
           debugvalue = F)

visr.param("enable_trajectories", label = "Single-cell trajectories", default = T,
           active.condition = sprintf("%s || (visr.param.enable_filtering && visr.param.enable_dim_red && visr.param.enable_clustering && visr.param.enable_de_analysis)", condition_monocle_object),
           info = "Discover cells transition from one state to another.",
           debugvalue = F)

visr.param("find_pseudotime_genes", label = "Find pseudotime changing genes", default = T,
           info = "Find genes that change as a function of pseudotime",
           active.condition = sprintf("%s || (visr.param.enable_filtering && visr.param.enable_dim_red && visr.param.enable_clustering && visr.param.enable_de_analysis && visr.param.enable_trajectories)", condition_monocle_object))

visr.param("analyze_trajectory_branches", label = "Analyze branches in trajectories", default = T,
           info = "Analyze branches in single-cell trajectories to identify the genes that differ at a particular branch point",
           active.condition = sprintf("%s || (visr.param.enable_filtering && visr.param.enable_dim_red && visr.param.enable_clustering && visr.param.enable_de_analysis && visr.param.enable_trajectories)", condition_monocle_object))
