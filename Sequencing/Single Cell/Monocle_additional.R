################################################################
visr.category("Additional parameters", collapsed = T,
              active.condition = "visr.param.output_dir != ''")
################################################################

visr.param("num_cores", label = "Number of cores to use for DE", min = 1L, default = 4L,
           info = "The number of cores to be used while testing each gene for differential expression.")

visr.param("num_histogram_bins", label = "Number of histogram bins", default = 50L,
           info = "Number of histogram bins in the histogram plots")
