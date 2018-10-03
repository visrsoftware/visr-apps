################################################################
visr.category("Filtering (Subsetting)",
              info = "Values used to filter genes and cells",
              active.condition = "visr.param.output_dir != '' && visr.param.enable_filtering")
################################################################

# Remove low (e.g. dead cells or empty wells) and high (e.g. doublets: made from two or more cells accidentally)
visr.param("filter_by_distribution", label = "Filter cells based on average distribution", default = TRUE,
           info = "Remove low (dead cells or empty wells) and high (doublets: made from two or more cells accidentally)")

visr.param("filter_sd_cutoff", label = "Filter range around average (s.d. units)", default = 2,
           info = "Units of standard deviation around average total count per cell to use as cut-off threshold",
           active.condition = "visr.param.filter_by_distribution")

GENE_SUBSET_METHOD_MEAN = "based on min average expression"
GENE_SUBSET_METHOD_PERCENT = "based on min % expressed cells"
visr.param("gene_subset_method", label = "Subset genes", items = c(GENE_SUBSET_METHOD_MEAN, GENE_SUBSET_METHOD_PERCENT),
           info = "Strategy used to subset genes")

visr.param("min_mean_expression", default = 0.1, min = 0,
           label = "Subset genes: minimum average expression", info = "Select genes which have an average expression of specified value",
           active.condition = sprintf("visr.param.gene_subset_method == '%s'", GENE_SUBSET_METHOD_MEAN))

visr.param("min_expressed_cells", default = 0.05, min = 0, max = 1,
           label = "Subset genes: minimum % expressed cells", info = "Select genes which have a minimum percentage of expressed cells",
           active.condition = sprintf("visr.param.gene_subset_method == '%s'", GENE_SUBSET_METHOD_PERCENT))

visr.param("min_expr", label = "Minimum gene expression", default = 1L,
           info = "The minimum expression threshold to be used to tally the number of cells expressing a gene and the number of genes expressed among all cells. A gene is 'expressed' if there is at least the specified count")
