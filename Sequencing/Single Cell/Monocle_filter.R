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

visr.param("min_expr", label = "Minimum gene expression threshold", default = 1L,
           info = "The minimum expression threshold to be used to tally the number of cells expressing a gene and the number of genes expressed among all cells. A gene is 'expressed' if there is at least the specified count",
           active.condition = sprintf("visr.param.gene_subset_method == '%s'", GENE_SUBSET_METHOD_PERCENT))

visr.param("min_expressed_cells", default = 0.05, min = 0, max = 1,
           label = "Subset genes: minimum % expressed cells", info = "Select genes which have a minimum percentage of expressed cells",
           active.condition = sprintf("visr.param.gene_subset_method == '%s'", GENE_SUBSET_METHOD_PERCENT))


#'
#' filter low-quality cells based on the distribution of total counts across the cells
#' based on: http://cole-trapnell-lab.github.io/monocle-release/docs/#filtering-low-quality-cells-recommended
#'
perform_filter_by_distribution <- function(monocle_app_object) {
  my_cds <- monocle_app_object$cds
  stopifnot(is(my_cds, "CellDataSet"))

  if (visr.var.add_title_pages)
    plotTitlePage("Filtering cells")

  pData(my_cds)$TotalCount <- Matrix::colSums(exprs(my_cds))
  my_cds <- my_cds[, pData(my_cds)$TotalCount < 1e6]
  upper_bound <- 10^(mean(log10(pData(my_cds)$TotalCount)) +
                       visr.param.filter_sd_cutoff * sd(log10(pData(my_cds)$TotalCount)))
  lower_bound <- 10^(mean(log10(pData(my_cds)$TotalCount)) -
                       visr.param.filter_sd_cutoff * sd(log10(pData(my_cds)$TotalCount)))

  is_filtered <- pData(my_cds)$TotalCount > lower_bound &
                 pData(my_cds)$TotalCount < upper_bound
  plot_title_for_filter <- sprintf('Filtering %d cells based on distribution of total expression per cell.\nOutput: %d cells in (mean -/+ %4.1f*sd) range',
                                   length(is_filtered), length(which(is_filtered)), visr.param.filter_sd_cutoff)

  my_cds <- monocle::detectGenes(my_cds, min_expr = visr.param.min_expr - 0.0001) # small -0.0001 value is subtracted to allow for >= min_expr
  pData(my_cds)$UMI <- Matrix::colSums(Biobase::exprs(my_cds))

  # scatterplot
  p <- ggplot(pData(my_cds), aes(UMI, num_genes_expressed)) +
    theme_bw() + geom_point(alpha = 0.2) +
    geom_vline(xintercept = c(lower_bound, upper_bound), linetype = "dotted", color = 'red') +
    xlab("total expression per cell") +
    ylab(sprintf("number of expressed genes (>= %4.2f)", visr.param.min_expr)) +
    ggtitle(plot_title_for_filter)
  print(p)

  # histogram
  p <- ggplot(pData(my_cds), aes(TotalCount)) +
    theme_bw() +
    geom_histogram(bins = visr.param.num_histogram_bins) +
    # geom_density() +
    geom_vline(xintercept = c(lower_bound, upper_bound), linetype = "dotted", color = 'red') +
    labs(x = "total expression per cell", y = "number of cells") +
    ggtitle(plot_title_for_filter)
  print(p)

  my_cds <- my_cds[, is_filtered]

  monocle_app_object$cds <- my_cds
  return(monocle_app_object)
}

#'
#' detect genes
#'
perform_detect_genes <- function(monocle_app_object) {
  my_cds <- monocle_app_object$cds
  stopifnot(is(my_cds, "CellDataSet"))

  if (visr.var.add_title_pages)
    plotTitlePage("Subsetting genes")

  my_cds <- monocle::detectGenes(my_cds, min_expr = visr.param.min_expr - 0.0001) # small -0.0001 value is subtracted to allow for >= min_expr
  # head(fData(my_cds)) #  number of cells expressing a particular gene
  summary(fData(my_cds)$num_cells_expressed)
  # head(pData(my_cds)) # The number of genes expressed per cell
  print(summary(pData(my_cds)$num_genes_expressed))

  p <- ggplot(pData(my_cds), aes(num_genes_expressed)) +
    theme_bw() +
    geom_histogram(bins = visr.param.num_histogram_bins) +
    labs(x = "number of expressed genes", y = "number of cells") +
    ggtitle(sprintf('Number of cells with given number of genes expressed per cell\n(min expression = %4.1f) ', visr.param.min_expr))
  print(p)

  if (FALSE) {
    # standardise to Z-distribution
    x <- pData(my_cds)$num_genes_expressed
    x_1 <- (x - mean(x)) / sd(x)
    df <- data.frame(x = x_1)
    p <- ggplot(data.frame(x = x_1), aes(x)) +
      theme_bw() +
      geom_histogram(bins = visr.param.num_histogram_bins) +
      labs(x = "number of expressed genes, standardised to Z-distribution", y = "number of cells") +
      geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = 'red')
    print(p)
  }

  #
  pData(my_cds)$UMI <- Matrix::colSums(Biobase::exprs(my_cds))
  p <- ggplot(pData(my_cds), aes(UMI, num_genes_expressed)) +
    theme_bw() + geom_point(alpha = 0.2) +
    xlab("total expression per cell") +
    ylab("number of expressed genes") +
    ggtitle(sprintf("Number of genes expressed ( >= %4.1f ) vs total expression for %d cells" ,visr.param.min_expr, nrow(pData(my_cds))))
  print(p)

  monocle_app_object$cds <- my_cds
  return(monocle_app_object)
}

#'
#' subset genes based on average expression and variability
#'
perform_subsetting_genes <- function(monocle_app_object) {
  my_cds <- monocle_app_object$cds
  stopifnot(is(my_cds, "CellDataSet"))

  # select gene based on their average expression and variability across cells.
  visr.logProgress("Calculating the mean and dispersion values for genes")
  if (visr.param.data_type != DATA_TYPE_UMI) {
    #TODO:
    stop("dispersionTable() only works for CellDataSet objects containing count-based expression data, either transcripts or reads.")
  }

  visr.logProgress("Deciding which genes to use in clustering the cells ...")

  disp_table <- monocle::dispersionTable(my_cds)
  # head(disp_table)

  is_genes_gt_mean_expression = disp_table$mean_expression >= visr.param.min_mean_expression

  p <-  ggplot(data.frame(x = disp_table$mean_expression), aes(x)) +
    theme_bw() +
    geom_histogram(bins = visr.param.num_histogram_bins) +
    scale_x_log10(labels=axis_plain_format) +
    labs(x = "average expression per gene (log scale)", y = "number of genes") +
    geom_vline(xintercept = visr.param.min_mean_expression, linetype = "dotted", color = 'red') +
    ggtitle(sprintf('Number of genes by average expression.\n%d out of %d have >= %4.1f average expression',
                    length(which(is_genes_gt_mean_expression)), length(is_genes_gt_mean_expression), visr.param.min_mean_expression))
  print(p)

  print(table(disp_table$mean_expression >= visr.param.min_mean_expression))

  # Clustering cells without marker genes: http://cole-trapnell-lab.github.io/monocle-release/docs/#clustering-cells
  # mark genes that will be used for clustering
  if (visr.param.gene_subset_method == GENE_SUBSET_METHOD_MEAN) {
    # subset genes based on minimum average expression
    unsup_clustering_genes <- subset(disp_table, mean_expression >= visr.param.min_mean_expression)
    subset_gene_ids <- unsup_clustering_genes$gene_id
    my_cds <- monocle::setOrderingFilter(my_cds, subset_gene_ids)
    min_threshold_used <- visr.param.min_mean_expression
  }
  else if (visr.param.gene_subset_method == GENE_SUBSET_METHOD_PERCENT) {
    # subset genes based on minimum number of expressed cells
    # num_cells_expressed is previously calculated in detectGenes based on min_expr value
    fData(my_cds)$use_for_ordering <- fData(my_cds)$num_cells_expressed > visr.param.min_expressed_cells * ncol(my_cds)
    min_threshold_used <- visr.param.min_expressed_cells
    subset_gene_ids <- fData(my_cds)$id[fData(my_cds)$use_for_ordering]
  }
  else {
    visr.assert_that(FALSE, sprintf("Invalid value for gene_subset_method: '%s'", visr.param.gene_subset_method))
  }

  # show genes marked for clustering
  p <- monocle::plot_ordering_genes(my_cds) +
    #theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_log10(labels = axis_plain_format) +
    scale_y_log10(labels = axis_plain_format) +
    labs(x = "average gene expression (log scale)", y = "dispersion (log scale)") +
    ggtitle(sprintf(
"Variability in a gene's expression depends on the average expression across cells.
The %d genes in black are marked for use in clustering %s %4.2f.
Red line shows expectation of the dispersion.",
      length(which(fData(my_cds)$use_for_ordering)), visr.param.gene_subset_method, min_threshold_used)) +
    theme(plot.title = element_text(size = 8))
  print(p)

  monocle_app_object$subset_gene_ids <- subset_gene_ids
  monocle_app_object$disp_table <- disp_table
  monocle_app_object$cds <- my_cds
  return(monocle_app_object)
}
