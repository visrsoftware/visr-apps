#TODO:
# - Save and load CellDataSet object
# - save cells.txt and genes.txt, de_genes, BEAM_res
# - subset cells based on barcode
# - add summaries to the report
# - classify cells based on gene-markers: http://cole-trapnell-lab.github.io/monocle-release/docs/#classifying-cells-by-type-recommended
#   Myoblast: [MYF5] >= 1, Fibroblast: [MYF5] < 1 & [ANPEP] > 1
# add filtering / tSNE / clustering results to the cells.txt

source("visrutils.R")

curr_dir <- "Sequencing/Single Cell"

visr.app.start("Monocle")

source(sprintf("%s/Utils.R",curr_dir))
source(sprintf("%s/Monocle_IO.R",curr_dir))
source(sprintf("%s/Monocle_filter.R",curr_dir))
source(sprintf("%s/Monocle_dim_reduction.R",curr_dir))
source(sprintf("%s/Monocle_cluster_cells.R",curr_dir))
source(sprintf("%s/Monocle_DE.R",curr_dir))
source(sprintf("%s/Monocle_construct_trajectories.R",curr_dir))
source(sprintf("%s/Monocle_pseudotime_genes.R",curr_dir))
source(sprintf("%s/Monocle_analyze_branches.R",curr_dir))
source(sprintf("%s/Monocle_additional.R",curr_dir))

visr.app.end()

stop()

visr.applyParameters()

########################################################################################################################
# End of parameter specification
########################################################################################################################

# workflow based on:
# http://cole-trapnell-lab.github.io/monocle-release/docs/
# https://davetang.org/muse/2017/10/01/getting-started-monocle/

visr.library("ggplot2")
visr.biocLite("pheatmap")
visr.biocLite("DDRTree")
# used to order genes
visr.library("dplyr")
try(visr.biocLite("monocle"))
visr.libraryGithub("monocle", "cole-trapnell-lab/monocle-release") # https://github.com/cole-trapnell-lab/monocle-release/issues/118

############################################################
#  Load input data
############################################################
enable_estimations <- TRUE # estimate dispersions
if (visr.param.input_type == INPUT_TYPE_10X) { # if data is 10x data
  visr.assert_file_exists(visr.param.data_dir_10x, '10X data directory')
  visr.assert_that(file.exists(paste0(visr.param.data_dir_10x, '/outs')), msg = paste("Cannot find 'outs' sub-directory inside the specified 10x directory:\n", visr.param.data_dir_10x))
  validateOutputDirectory(visr.param.output_dir, visr.param.create_subdir)

  visr.librarySource("cellrangerRkit", "http://cf.10xgenomics.com/supp/cell-exp/rkit-install-2.0.0.R")

  visr.logProgress("Loading data from the 10x pipeline")
  gbm <- load_cellranger_matrix(visr.param.data_dir_10x)
  # class(gbm)
  # dim(Biobase::exprs(gbm))
  Biobase::exprs(gbm)[1:10, 1:4] # the first 10 genes in the first 4 cells

  # dim(pData(gbm)) # the phenotypic data
  # head(pData(gbm))

  # dim(fData(gbm)) # the feature information
  # head(fData(gbm)) # this table will be useful for matching gene IDs to symbols

  # create a CellDataSet object
  # note: Monocle expects that the gene symbol column in the feature data is called gene_short_name
  my_feat <- fData(gbm)
  names(my_feat) <- c('id', 'gene_short_name') # rename gene symbol column
  familyFunction <- negbinomial.size() # appropriate expression family for UMI data
  if (visr.param.data_type == DATA_TYPE_FPKM) {
    #TODO: http://cole-trapnell-lab.github.io/monocle-release/docs/#converting-tpm-fpkm-values-into-mrna-counts-alternative
    visr.message("DATA_TYPE_FPKM not implemented")
    familyFunction <- tobit() #TODO: parameter for tobit()
    enable_estimations <- FALSE
  } else if (visr.param.data_type == DATA_TYPE_LOG_FPKM) {
    # TODO
    visr.message("DATA_TYPE_LOG_FPKM not implemented")
    familyFunction <- gaussianff()
    enable_estimations <- FALSE
  }

  #TODO:If your data contains relative counts (e.g. FPKM or TPM values), use relative2abs() to convert these measurements into absolute counts: http://cole-trapnell-lab.github.io/monocle-release/docs/#converting-tpm-fpkm-values-into-mrna-counts-alternative
  visr.logProgress("Creating CellDateSet object")
  my_cds <<- newCellDataSet(Biobase::exprs(gbm),
                            phenoData = new("AnnotatedDataFrame", data = pData(gbm)),
                            featureData = new("AnnotatedDataFrame", data = my_feat),
                            lowerDetectionLimit = 0.5,
                            expressionFamily = familyFunction)
  # my_cds
  # slotNames(my_cds)
} else if (visr.param.input_type == INPUT_TYPE_COUNT_MATRIX_TXT) {
  visr.assert_file_exists(visr.param.expression_matrix, 'expression matrix')
  visr.assert_file_exists(visr.param.sample_sheet, 'sample sheet')
  visr.assert_file_exists(visr.param.gene_annotation, 'gene annotation')
  validateOutputDirectory()

  visr.logProgress("Loading data from the expression matrix")

  expr_matrix <- read.table(visr.param.expression_matrix) # e.g. "fpkm_matrix.txt"
  sample_sheet <- read.delim(visr.param.sample_sheet) # e.g. "cell_sample_sheet.txt"
  gene_annotation <- read.delim(visr.param.gene_annotation) # e.g. "gene_annotations.txt"

  visr.assert_that(nrow(expr_matrix) == ncol(sample_sheet), msg = "Expression value matrix does not have the same number of columns as the sample sheet (pheno data) has rows.")
  visr.assert_that(nrow(expr_matrix) == nrow(gene_annotation), msg = "Expression value matrix does not have the same number of rows as the gene_annotation (feature data) data frame has rows.")
  visr.assert_that(all(rownames(expr_matrix) == colnames(sample_sheet)), msg = "row names of the sample sheet (pheno data) does not match the column names of the expression matrix.")
  visr.assert_that(all(rownames(expr_matrix) == rownames(gene_annotation)), msg = "row names of the gene annotation (feature data) should match row names of the expression matrix.")
  visr.assert_that("gene_short_name" %in% rownames(gene_annotation), msg = "one of the columns of the gene annotation (feature data) should be named 'gene_short_name'")


  pd <- new("AnnotatedDataFrame", data = sample_sheet)
  fd <- new("AnnotatedDataFrame", data = gene_annotation)
  my_cds <<- newCellDataSet(as.matrix(expr_matrix), phenoData = pd, featureData = fd)

} else if (visr.param.input_type == INPUT_TYPE_SPARSE_MATRIX) {
  visr.message("input_type INPUT_TYPE_SPARSE_MATRIX Not Implemented Yet")
  validateOutputDirectory()
}

############################################################
#  Prepare output
############################################################

# output current parameters into a file
writeParameters(visr.param.output_dir)

# prepare the pdf to save the plots to
output_plot_file <- startReport(visr.param.output_dir)


############################################################
#  perform normalization and variance estimation steps
############################################################

if (enable_estimations) {
  # perform normalization and variance estimation steps, which will be used in the differential expression analyses later on.
  # estimateSizeFactors() and estimateDispersions() will only work, and are only needed, if you are working with a CellDataSet
  # with a negbinomial() or negbinomial.size() expression family.
  visr.logProgress("Performing normalization and variance estimation (takes a few minutes)")
  my_cds <- estimateSizeFactors(my_cds)
  my_cds <- estimateDispersions(my_cds)
}

################################################################
# Filter low-quality cells
# http://cole-trapnell-lab.github.io/monocle-release/docs/#filtering-low-quality-cells-recommended
################################################################

# Remove low (e.g. dead cells or empty wells) and high (e.g. doublets: made from two or more cells accidentally)
if (visr.param.filter_by_distribution) {
  plotTitlePage("Filtering cells")

  # filter based on the distribution of total counts across the cells
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
}

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

################################################################
# Subsetting genes
################################################################

if (visr.param.enable_clustering && ncol(my_cds) > 5e+05 && visr.param.cluster_method != CLUSTER_METHOD_LOUVAIN) {
  visr.message("Number of cells in your data is larger than 50k. clusterCells with 'densityPeak' or 'DDRTree' may crash. Please try to use the 'Louvain' clustering algorithm!", type = "warning")
}

if (visr.param.enable_filtering) {
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
      ggtitle(sprintf('Number of genes by average expression. %d out of %d have >= %4.1f average expression',
                      length(which(is_genes_gt_mean_expression)), length(is_genes_gt_mean_expression), visr.param.min_mean_expression))
  print(p)

  visr.writeDataTable(disp_table, paste0(visr.param.output_dir, "/genes_dispersion.txt"))

  print(table(disp_table$mean_expression >= visr.param.min_mean_expression))

  # Clustering cells without marker genes: http://cole-trapnell-lab.github.io/monocle-release/docs/#clustering-cells
  # mark genes that will be used for clustering
  if (visr.param.gene_subset_method == GENE_SUBSET_METHOD_MEAN) {
    # subset genes based on minimum average expression
    unsup_clustering_genes <- subset(disp_table, mean_expression >= visr.param.min_mean_expression)
    clustering_gene_ids <- unsup_clustering_genes$gene_id
    my_cds <- monocle::setOrderingFilter(my_cds, clustering_gene_ids)
    min_threshold_used <- visr.param.min_mean_expression

  } else if (visr.param.gene_subset_method == GENE_SUBSET_METHOD_PERCENT) {
    # subset genes based on minimum number of expressed cells
    fData(my_cds)$use_for_ordering <- fData(my_cds)$num_cells_expressed > visr.param.min_expressed_cells * ncol(my_cds)
    min_threshold_used <- visr.param.min_expressed_cells
    clustering_gene_ids <- fData(my_cds)$id[fData(my_cds)$use_for_ordering]
  } else {
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
}

################################################################
# dimensionality reduction
################################################################

if (visr.param.enable_dim_red) {
  plotTitlePage("Dimensionality Reduction")

  visr.logProgress("Plotting the percentage of variance explained by each component (takes a few minutes)")
  p <-  monocle::plot_pc_variance_explained(my_cds, return_all = FALSE) + #TODO: norm_method= c("log", "vstExprs", "none")
    ggtitle("Percentage of variance explained by each component\nbased on a PCA performed on the normalised expression data") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept = visr.param.reduce_num_dim, linetype = "dotted", color = 'red') +
    labs(x = "PCA Components", y = "Variance explained by each PCA component")
  print(p)

  visr.logProgress("Computing a projection of a CellDataSet object into a lower dimensional space (takes a few minutes)")
  my_cds <- monocle::reduceDimension(my_cds, max_components = 2, num_dim = visr.param.reduce_num_dim,
                            reduction_method = visr.param.reduction_method, verbose = TRUE) # norm_method = c("log", "vstExprs", "none")

  p <- ggplot(data.frame(t(monocle::reducedDimA(my_cds)))) +
    geom_point(aes(X1, X2)) +
    theme_bw() +
    xlab("Dimension1") + ylab("Dimension2") +
    ggtitle(sprintf("Dimensionality reduction using %s method", visr.param.reduction_method))
  print(p)
}

################################################################
# Clustering
################################################################

if (visr.param.enable_clustering) {
  plotTitlePage("Clustering")

  num_clusters = NULL
  if (visr.param.cluster_method == CLUSTER_METHOD_DENSITY_PEAK) {
    # for some reason it alwasy produces one less cluster than specified value. so add 1 for now
    num_clusters = if (!is.null(visr.param.num_clusters)) (visr.param.num_clusters + 1) else NULL
  } else if (visr.param.cluster_method == CLUSTER_METHOD_DDRTREE) {
    num_clusters = visr.param.num_centers
  }

  visr.logProgress(paste("Performing unsupervised clustering. method:", visr.param.cluster_method))
  my_cds <- monocle::clusterCells(my_cds, verbose = T,
    skip_rho_sigma = F, #visr.param.skip_rho_sigma, # whether you want to skip the calculation of rho / sigma
    num_clusters = num_clusters,
    inspect_rho_sigma = F, # whether you want to interactively select the rho and sigma
    rho_threshold = visr.param.rho_threshold,
    delta_threshold = visr.param.delta_threshold,
    gaussian = visr.param.gaussian, # whether use Gaussian kernel for calculating the local density
    peaks = NULL, # numeric vector indicating the index of density peaks used for clustering.
    cell_type_hierarchy = NULL,
    frequency_thresh = NULL,
    enrichment_thresh = NULL,
    clustering_genes = NULL,
    k = visr.param.louvain_k,
    louvain_iter = visr.param.louvain_iter,
    weight = visr.param.louvain_weight,
    method = visr.param.cluster_method,
    verbose = T)

  head(pData(my_cds))
  print("Cluster sizes:")
  print(table(pData(my_cds)$Cluster))
  if (!is.null(visr.param.rho_threshold) & !is.null(visr.param.delta_threshold))
    title_peak <- sprintf("Detected %d peaks using rho=%4.2f and delta = %4.2f ",
                          length(which(pData(my_cds)$peaks)),
                          visr.param.rho_threshold, visr.param.delta_threshold)
  else
    title_peak <- sprintf("Detected %d peaks using top 95%% delta and top 95%% rho for threshold.",
                          length(which(pData(my_cds)$peaks)))

  if (visr.param.cluster_method == CLUSTER_METHOD_DENSITY_PEAK) {
    visr.logProgress("Plotting the decision map of density clusters (delta vs. rho)")
    p <- monocle::plot_rho_delta(my_cds, rho_threshold = visr.param.rho_threshold, delta_threshold = visr.param.delta_threshold) +
      ggtitle(paste0("Decision map of density clusters\nPeaks are cells with high local density that are far away from other cells with high local density\n", title_peak)) +
      theme(plot.title = element_text(size = 12)) + # hjust = 0.5
      labs(x = "rho: local density", y = "delta: local distance (to another cell with higher density)")
    print(p)
  }

  p <- monocle::plot_cell_clusters(my_cds) +
        xlab("tSNE1") + ylab("tSNE2") +
    ggtitle(sprintf("Unsupervised clustering of cells\nusing %s method", visr.param.cluster_method))
  print(p)
}

#############################################
# Differential Expression Analysis
#############################################

perform_de <- function(cds, fullModelFormulaStr, de_genes_filename) {
  de_genes <- monocle::differentialGeneTest(cds_de_subset,
                                            fullModelFormulaStr = fullModelFormulaStr,
                                            reducedModelFormulaStr = "~1", # default
                                            relative_expr = TRUE, # default
                                            cores = visr.param.num_cores,
                                            verbose = TRUE)
  dim(de_genes)
  #de_genes_valid <- de_genes[which(de_genes$qval < visr.param.trajectory_max_qval),]
  #de_genes_valid <- de_genes_valid[order(de_genes_valid$qval),]
  de_genes_valid <- de_genes[order(de_genes$qval),]

  visr.logProgress(paste("Writing DE genes to file", de_genes_filename))
  visr.writeDataTable(de_genes_valid, de_genes_filename)

  # the most statistically significant genes
  visr.logProgress(paste("Plotting level of expression for top", visr.param.num_plot_genes_jitter, "most statistically significant gene(s)."))
  p <- monocle::plot_genes_jitter(my_cds[de_genes_valid$id[seq(visr.param.num_plot_genes_jitter)],],
                                  grouping = "Cluster", color_by = "Cluster",
                                  ncol = ceiling(sqrt(visr.param.num_plot_genes_jitter)), nrow = NULL)
  print(p)

  return (de_genes_valid)
}

if (visr.param.enable_de_analysis) {

  plotTitlePage("Differential Expression Analysis", 2)

  if (visr.param.de_subset_by_marker_genes & visr.param.marker_genes_list != '') {
    marker_genes <- row.names(subset(fData(my_cds), gene_short_name %in% strsplit(visr.param.marker_genes_list, split = "[ ]*[,| ]+[ ]*")[[1]]))
    cds_de_subset <- my_cds[marker_genes,]
  } else {
    cds_de_subset <- my_cds[clustering_gene_ids,]
  }

  if (visr.param.de_formula == DE_FORMULA_ALL_CLUSTERS) {
    visr.logProgress("Performing differential expression analysis across all clusters ...")
    de_genes_filename <- paste0(visr.param.output_dir, "/de_genes_all_clusters.txt")
    de_genes <<- perform_de(cds_de_subset, '~Cluster', de_genes_filename)
  }
  else if (visr.param.de_formula == DE_FORMULA_ONE_CLUSTER) {
    # create vector of no's
    de_cluster_ids = if (visr.param.de_cluster_id == Inf) levels(pData(my_cds)$Cluster) else visr.param.de_cluster_id
    i <- 0
    for (de_cluster_id in de_cluster_ids) {
      i <- i + 1
      # change status to yes if the cell was in cluster 1
      my_vector <- rep('no', nrow(pData(my_cds)))
      my_vector[pData(my_cds)$Cluster == de_cluster_id] <- 'yes' #rep('yes', sum(pData(my_cds)$Cluster == de_cluster_id))

      # add vector to phenoData
      pData(my_cds)$test <- my_vector

      head(pData(my_cds))

      # TODO: perform DE based on a selected gene (replace ~test with the specified gene expresion column)
      # question: should we use normalized or unnormalized expression values
      visr.logProgress(paste("(", i, "of", length(de_cluster_ids), ")",
                             "Performing the differential expression analysis on cluster", de_cluster_id,
                             "\nfor the", nrow(unsup_clustering_genes)))

      de_genes_filename <- paste0(visr.param.output_dir, "/de_genes_cluster", de_cluster_id, ".txt")

      de_genes <<- perform_de(cds_de_subset, '~test', de_genes_filename)
    }
  }
  else if (visr.param.de_formula == DE_FORMULA_PHENOTYPE) {
    de_genes_filename <- paste0(visr.param.output_dir, "/de_genes_", visr.param.de_phenotype_name, ".txt")
    de_genes <<- perform_de(cds_de_subset, paste0('~', visr.param.de_phenotype_name), de_genes_filename)
  }
}

#############################################
# Constructing Single Cell Trajectories
#############################################

if (visr.param.enable_trajectories) {
  plotTitlePage("Constructing single cell trajectories", 1)

  # Step 1: determine ordering genes based on genes that differ between clusters
  # get DE genes
  if (!is.null(visr.param.trajectory_max_qval))
    trajectory_ordering_genes <- de_genes[which(de_genes$qval <= visr.param.trajectory_max_qval),]
  else
    trajectory_ordering_genes <- de_genes[which(is.na(de_genes$qval)),]
  # pick top 1,000
  trajectory_ordering_genes <- trajectory_ordering_genes[1:min(visr.param.trajectory_num_genes, nrow(trajectory_ordering_genes)),]
  my_cds <- monocle::setOrderingFilter(my_cds, ordering_genes = trajectory_ordering_genes$id)

  p <- monocle::plot_ordering_genes(my_cds) +
    scale_x_log10(labels=axis_plain_format) +
    scale_y_log10(labels=axis_plain_format) +
    labs(x = "average gene expression (log scale)", y = "dispersion (log scale)") +
    ggtitle(sprintf("Step 1: selected %d genes for ordering", nrow(trajectory_ordering_genes)))
  print(p)

  # Step 2: reducing dimension using DDRTree (Reversed Graph Embedding)
  visr.logProgress(sprintf("Reducing the dimensionaly of %d ordering genes using 'DDRTree' method ...",
                           length(trajectory_ordering_genes$id)))
  my_cds <- monocle::reduceDimension(my_cds,
                                     max_components = 2,
                                     reduction_method = 'DDRTree', # c("DDRTree", "ICA", "tSNE", "SimplePPT", "L1-graph", "SGL-tree")
                                     norm_method = 'log',
                                     verbose = T)

  p <- ggplot(data.frame(t(monocle::reducedDimS(my_cds)))) +
    geom_point(aes(X1, X2)) +
    theme_bw() +
    xlab("Component1") + ylab("Component2") +
    ggtitle("Step 2: reducing dimension using DDRTree method")
  print(p)

  # Step 3: ordering the cells in pseudotime (fit the best tree it can to the data)
  visr.logProgress(sprintf("Ordering %d cells over the trajectory...", ncol(my_cds)))
  my_cds <- monocle::orderCells(my_cds)

  ## plot trajectories
  p <- monocle::plot_cell_trajectory(my_cds , color_by = "Pseudotime") +
    ggtitle("Step 3: ordering cells in pseudotime")
  print(p)

  p <- monocle::plot_cell_trajectory(my_cds , color_by = "State")
  print(p)
  if ("Cluster" %in% names(pData(my_cds))) {
    p <- monocle::plot_cell_trajectory(my_cds , color_by = "Cluster")
    print(p)
  }

  if (visr.param.find_pseudotime_genes) {
    #############################################
    visr.logProgress("Finding genes that change as a function of pseudotime ...")
    #############################################
    my_pseudotime_de <- monocle::differentialGeneTest(
      my_cds, fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = visr.param.num_cores)


    my_pseudotime_de %>% arrange(qval) %>% head()

    my_pseudotime_de %>% arrange(qval) %>% head() %>% select(id) -> my_pseudotime_gene
    my_pseudotime_gene <- my_pseudotime_gene$id

    p <- monocle::plot_genes_in_pseudotime(my_cds[my_pseudotime_gene,])
    print(p)
    # todo: parameters
  }

  if (visr.param.cluster_genes_by_pseudo) {
    #############################################
    visr.logProgress("Clustering genes by pseudotemporal expression pattern")
    #############################################

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

    # todo: add the cluster ids to genes.txt
    visr.logProgress("Extracting the genes for each cluster ...")
    my_cluster <- cutree(my_pseudotime_cluster$tree_row, visr.param.num_pseudo_gene_clusters)
    for (cluster_id in seq(visr.param.num_pseudo_gene_clusters)) {
      print(sprintf("Genes for cluster %d", cluster_id))
      print(my_pseudotime_de[names(my_cluster[my_cluster == cluster_id]), "gene_short_name"])
    }
  }

  # A table of genes is returned with significance values that indicate whether genes have expression patterns that are branch dependent.
  if (visr.param.analyze_trajectory_branches) {
    plotTitlePage("Analyzing branches in single-cell trajectories", 1)

    trajectory_branch_point <- visr.param.trajectory_branch_point
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
  }
}

#############################################
# Output tables
#############################################
finishReport()

browseURL(output_plot_file)

cds_dims <- data.frame(t(monocle::reducedDimA(my_cds)))
colnames(cds_dims) <- c("Component1", "Component2")
cell_table <- cbind(pData(my_cds), cds_dims)
visr.writeDataTable(cell_table, paste0(visr.param.output_dir, "/cells.txt"))

cachedResults <- list(params = visr.getParams(), my_cds = my_cds, disp_table = disp_table)
save(cachedResults, file = paste0(visr.param.output_dir, "/cachedResults.RData"))
