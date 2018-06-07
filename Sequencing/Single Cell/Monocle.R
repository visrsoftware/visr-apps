source("visrutils.R")

visr.app.start("Monocle")

################################ Input

visr.category(label="Input")

INPUT_TYPE_10X = "10X single cell dataset"
INPUT_TYPE_COUNT_MATRIX_TXT = "Count matrix (.txt)"
INPUT_TYPE_SPARSE_MATRIX = "Sparse matrix"
visr.param("input_type",
           info = "Specify the format of your input data. Whether it is a 10X cell ranger dataset or a normal or sparse an expression matrix",
           items = c(INPUT_TYPE_10X, INPUT_TYPE_COUNT_MATRIX_TXT, INPUT_TYPE_SPARSE_MATRIX), debugvalue = INPUT_TYPE_10X)

visr.param("data_dir_10x",
           label="10X dataset directory",
           info="10X Cell Ranger dataset directory. It should contain a sub-directory named 'outs'",
           type="filename", filename.mode = "dir",
           active.condition = sprintf("visr.param.input_type == '%s'", INPUT_TYPE_10X),
           debugvalue= "~/SFU/Datasets/SingleCell/pbmc3k/") # Peripheral Blood Mononuclear Cells (PBMCs) from a healthy donor

visr.param("expression_matrix", type = "filename",
           info = "Numeric matrix of expression values, where rows are genes, and columns are cells",
           active.condition = sprintf("(visr.param.input_type == '%s') || (visr.param.input_type == '%s')", INPUT_TYPE_COUNT_MATRIX_TXT, INPUT_TYPE_SPARSE_MATRIX))

visr.param("sample_sheet", type = "filename",
           info = "Pheno data where rows are cells, and columns are cell attributes (such as cell type, culture condition, day captured, etc.). It should have same number of rows as the number of columns in expression value matrix.",
           active.condition = sprintf("(visr.param.input_type == '%s') || (visr.param.input_type == '%s')", INPUT_TYPE_COUNT_MATRIX_TXT, INPUT_TYPE_SPARSE_MATRIX))

visr.param("gene_annotation", type = "filename",
           info = "Feature data where rows are features (e.g. genes), and columns are gene attributes (such as biotype, gc content, etc.). It should have same number of rows as the number of rows in expression value matrix.",
           active.condition = sprintf("(visr.param.input_type == '%s') || (visr.param.input_type == '%s')", INPUT_TYPE_COUNT_MATRIX_TXT, INPUT_TYPE_SPARSE_MATRIX))

DATA_TYPE_UMI = "UMIs, Raw transcript counts"
DATA_TYPE_FPKM = "FPKM, TPM"
DATA_TYPE_LOG_FPKM = "log FPKM/TPMs, Ct values from SC qPCR"
visr.param("data_type",
           label = "Data type",
           info = "Specify the appropriate data type. Monocle works well with both count-based measures (e.g. UMIs) and relative expression data.",
           items=c(DATA_TYPE_UMI, DATA_TYPE_FPKM, DATA_TYPE_LOG_FPKM),
           active.condition = "(visr.param.data_dir_10x != '') || (visr.param.expression_matrix != '' && visr.param.sample_sheet != '' && visr.param.gene_annotation != '')")

################################ Output tables

visr.category("Output",
              active.condition = "(visr.param.data_dir_10x != '') || (visr.param.expression_matrix != '' && visr.param.sample_sheet != '' && visr.param.gene_annotation != '')")

visr.param("output_dir",
           label="Output Directory to save the results",
           info="Output directory where the analysis results will be saved to",
           type="filename", filename.mode = "dir",
           debugvalue= "~/SFU/Datasets/SingleCell/output/")

visr.param("create_subdir", default = TRUE,
           label = "Create new sub-direcory",
           info = "Create a new sub directory with the name DATE_TIME (YYYYMMDD_hhmmss)")

################################ Analysis Steps

visr.category("Analysis steps",
              info = "Different analysis steps",
              active.condition = "visr.param.output_dir != '' && ((visr.param.data_dir_10x != '') || (visr.param.expression_matrix != '' && visr.param.sample_sheet != '' && visr.param.gene_annotation != ''))")

visr.param("enable_subsetting", label = "Subsetting and dimensionality reduction", default = T)

visr.param("enable_clustering", label = "Clustering", default = T,
           active.condition = "visr.param.enable_subsetting")

visr.param("enable_de_analysis", label = "Differential expression analysis", default = F,
           active.condition = "visr.param.enable_subsetting && visr.param.enable_clustering")

visr.param("enable_trajectories", label = "Single-cell trajectories", default = F,
           active.condition = "visr.param.enable_subsetting && visr.param.enable_clustering && visr.param.enable_de_analysis")

################################ Subsetting and dimensionality reduction

visr.category("Subsetting and dimensionality reduction",
              info = "Values used to reduce dimensionality of data (Genes and Cells)",
              active.condition = "visr.param.output_dir != '' && visr.param.enable_subsetting")

visr.param("min_expr", label = "Expression threshold", default = 1L,
           info = "The expression threshold to be used to tally the number of cells expressing a gene and the number of genes expressed among all cells. A gene is 'expressed' if there is at least the specified count")

# number of histogram bins in the histogram plots
visr.param.num_histogram_bins <- 50L

GENE_SUBSET_METHOD_MEAN = "Based on mean expression"
GENE_SUBSET_METHOD_PERCENT = "Based on number of expressed cells"
visr.param("gene_subset_method", label = "Gene subsetting method", items = c(GENE_SUBSET_METHOD_MEAN, GENE_SUBSET_METHOD_PERCENT),
           info = "Strategy used to subset genes")

visr.param("min_mean_expression", default = 0.1, min = 0,
           label = "Subset genes: Minimum mean expression", info = "Select genes which have a mean expression of specified value",
           active.condition = sprintf("visr.param.gene_subset_method == '%s'", GENE_SUBSET_METHOD_MEAN))

visr.param("min_expressed_cells", default = 0.05, min = 0, max = 1,
           label = "Subset genes: Minimum percentage of expressed cells", info = "Select genes which have a minimum percentage of expressed cells",
           active.condition = sprintf("visr.param.gene_subset_method == '%s'", GENE_SUBSET_METHOD_PERCENT))

visr.param("reduce_num_dim", label = "Number of reduced dimensions", default = 50L, min = 2L)

visr.param("reduction_method", label = "Algorithm for dimensionality reduction", info = "Algorithm to use for dimensionality reduction", items = c("DDRTree", "ICA", "tSNE", "SimplePPT", "L1-graph", "SGL-tree"), default = 'tSNE')

################################ Clustering

visr.category("Clustering",
              info = "Identify new (and possibly rare) subtypes of cells in Single-cell RNA-Seq experiments.",
              active.condition = "visr.param.output_dir != '' && visr.param.enable_clustering")

CLUSTER_METHOD_DENSITY_PEAK = "densityPeak"
CLUSTER_METHOD_LOUVAIN = "louvain"
CLUSTER_METHOD_DDRTREE = "DDRTree"
visr.param("cluster_method", items = c(CLUSTER_METHOD_DENSITY_PEAK, CLUSTER_METHOD_LOUVAIN, CLUSTER_METHOD_DDRTREE),
           info = "Method for clustering cells. For big datasets (like data with 50 k cells or so), we recommend using the 'louvain' clustering algorithm.")

visr.param("num_clusters", label = "Number of clusters", type = "integer", min = 1L, default = "NULL", items = c("NULL"), item.labels = c("auto"), debugvalue = NULL,
           info = "Number of clusters. When auto, Use 0.95 of the delta and 0.95 of the rho as the cutoff for assigning density peaks and clusters",
           active.condition = sprintf("visr.param.cluster_method == '%s'", CLUSTER_METHOD_DENSITY_PEAK))

visr.param("skip_rho_sigma", default = FALSE, info = "skip the calculation of rho / sigma",
           active.condition = sprintf("visr.param.cluster_method == '%s'", CLUSTER_METHOD_DENSITY_PEAK))

visr.param("rho_threshold", label = "rho (cell's local density) threshod", type = "double", min = 0,
           items = "NULL", item.labels = "auto", default = "NULL", debugvalue = NULL,
           info = "The threshold of cell's local density (rho) used to select the density peaks",
           active.condition = sprintf("visr.param.cluster_method == '%s'", CLUSTER_METHOD_DENSITY_PEAK))

visr.param("delta_threshold", label = "delta (local distance) threshod", type = "double", min = 0,
           items = "NULL", item.labels = "auto", default = "NULL", debugvalue = NULL,
           info = "The threshold of local distance (nearest distance of a cell to another cell with higher distance) used to select the density peaks",
           active.condition = sprintf("visr.param.cluster_method == '%s'", CLUSTER_METHOD_DENSITY_PEAK))

visr.param("gaussian", label = "Use gaussian kernel?", default = T,
           info = "Whether or not Gaussian kernel will be used for calculating the local density in densityClust function",
           active.condition = sprintf("visr.param.cluster_method == '%s'", CLUSTER_METHOD_DENSITY_PEAK))

visr.param("num_centers", label = "Number of centroids", type = "integer", min = 1L, default = 3L,
           info = "Number of number of centroids passed to DDRTree ('Dimensionality Reduction via Graph Structure Learning' by Qi Mao, et al)",
           active.condition = sprintf("visr.param.cluster_method == '%s'", CLUSTER_METHOD_DDRTREE))

visr.param("louvain_k", default = 50L,
           info = "number of kNN used in creating the k nearest neighbor graph for Louvain clustering. The number of kNN is related to the resolution of the clustering result, bigger number of kNN gives low resolution and vice versa.",
           active.condition = sprintf("visr.param.cluster_method == '%s'", CLUSTER_METHOD_LOUVAIN))

visr.param("louvain_iter", default = 1L,
           info = "number of iterations used for Louvain clustering. The clustering result gives the largest modularity score will be used as the final clustering result.",
           active.condition = sprintf("visr.param.cluster_method == '%s'", CLUSTER_METHOD_LOUVAIN))

visr.param("louvain_weight", default = F,
           info = "Use Jaccard coefficent for two nearest neighbors (based on the overlapping of their kNN) as the weight used for Louvain clustering.",
           active.condition = sprintf("visr.param.cluster_method == '%s'", CLUSTER_METHOD_LOUVAIN))

################################ Differential expression analysis

visr.category("Differential expression analysis",
              info = "A sophisticated but easy to use system for differential expression, to characterize new cell types and states by comparing them to other",
              active.condition = "visr.param.output_dir != '' && visr.param.enable_de_analysis")


DE_FORMULA_CLUSTER = "Cluster"
DE_FORMULA_PHENOTYPE = "Phenotype (cell attribute)"
visr.param("de_formula", label = "Perform DE by", items = c(DE_FORMULA_CLUSTER, DE_FORMULA_PHENOTYPE))

visr.param("de_cluster_id", label = "Perform DE on which Cluster id?", min = 1L, default = 1L,
           items = c(Inf), item.labels = c("All (slow)"),
           active.condition = sprintf("visr.param.de_formula == '%s'", DE_FORMULA_CLUSTER))

visr.param("de_phenotype_name", label = "Phenotype (cell attribute) name",
           active.condition = sprintf("visr.param.de_formula == '%s'", DE_FORMULA_PHENOTYPE))

visr.param("de_subset_by_marker_genes", label = "Subset using known marker genes", default = F)

visr.param("marker_genes_list", label = "Marker genes (comma separated)", info="Optional comma separated list of short gene names to use as marker genes.",
           debugvalue = "MEF2C, MEF2D, MYF5, ANPEP, PDGFRA, MYOG, TPM1, TPM2, MYH2, MYH3, NCAM1, TNNT1, TNNT2, TNNC1, CDK1, CDK2, CCNB1, CCNB2, CCND1, CCNA1, ID1",
           active.condition = "visr.param.de_subset_by_marker_genes")

visr.param("num_plot_genes_jitter", label = "Draw level of expression for how many top genes?", default = 9L,
           info = "Plots the level of expression for each group of cells per gene,\nfor the specified number of most statistically significant genes.") # min = 0

################################ Single-cell trajectories

visr.category("Single-cell trajectories",
              info = "Discover cells transition from one state to another, in development, disease, and throughout life.",
              active.condition = "visr.param.output_dir != '' && visr.param.enable_trajectories")

visr.app.end()

visr.applyParameters()

########################################################################################################################
# End of parameter specification
########################################################################################################################

try(visr.biocLite("monocle"))
visr.libraryGithub("monocle", "cole-trapnell-lab/monocle-release") # https://github.com/cole-trapnell-lab/monocle-release/issues/118
# used to order genes
visr.library("dplyr")

############################################################
#  Validate output directory
############################################################

validateOutputDirectory <- function() {
  visr.assert_that(!is.null(visr.param.output_dir) && visr.param.output_dir != "", msg = "Output directory not specified")
  if (visr.param.create_subdir)
    visr.param.output_dir <<- paste0(visr.param.output_dir, format(Sys.time(), "/%Y%m%d_%H%M%S"))
  dir.create(visr.param.output_dir, recursive = T, showWarnings = F)
  visr.assert_that(dir.exists(visr.param.output_dir), msg = sprintf("Output directory '%s' is not valid", visr.param.output_dir))
}

############################################################
#  Load input data
############################################################

if (visr.param.input_type == INPUT_TYPE_10X) { # if data is 10x data
  visr.assert_file_exists(visr.param.data_dir_10x, '10X data directory')
  visr.assert_that(file.exists(paste0(visr.param.data_dir_10x, '/outs')), msg = paste("Cannot find 'outs' sub-directory inside the specified 10x directory:\n", visr.param.data_dir_10x))
  validateOutputDirectory()

  visr.librarySource("cellrangerRkit", "http://cf.10xgenomics.com/supp/cell-exp/rkit-install-2.0.0.R")

  visr.logProgress("Loading data from the 10x Genomics Cell Ranger pipeline")
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
  familyFunction <- negbinomial.size()
  if (visr.param.data_type == DATA_TYPE_FPKM) {
    #TODO: http://cole-trapnell-lab.github.io/monocle-release/docs/#converting-tpm-fpkm-values-into-mrna-counts-alternative
    familyFunction <- tobit() #TODO: parameter for tobit()
    visr.message("DATA_TYPE_FPKM not implemented")
  } else if (visr.param.data_type == DATA_TYPE_LOG_FPKM) {
    # TODO
    familyFunction <- gaussianff()
    visr.message("DATA_TYPE_LOG_FPKM not implemented")
  }

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

# output current parameters into a file
write(capture.output(print(visr.getParams())), file = paste0(visr.param.output_dir, "/visr_params.txt"))

# prepare the pdf to save the plots to
visr.dev <- dev.cur()
if (exists("visr_app_pdf_dev") && !is.null(visr_app_pdf_dev)) {
  dev.off(which = visr_app_pdf_dev)
}
output_plot_file <- paste0(visr.param.output_dir, "/plots.pdf")
pdf(file=output_plot_file)
visr_app_pdf_dev <- dev.cur()


# perform normalization and variance estimation steps, which will be used in the differential expression analyses later on.
visr.logProgress("Performing normalization and variance estimation (takes a few minutes)")
my_cds <- estimateSizeFactors(my_cds)
my_cds <- estimateDispersions(my_cds)

my_cds <- detectGenes(my_cds, min_expr = visr.param.min_expr - 0.0001)
# head(fData(my_cds)) #  number of cells expressing a particular gene
summary(fData(my_cds)$num_cells_expressed)
# head(pData(my_cds)) # The number of genes expressed per cell
print(summary(pData(my_cds)$num_genes_expressed))


# standardise to Z-distribution
x <- pData(my_cds)$num_genes_expressed
x_1 <- (x - mean(x)) / sd(x)
visr.library("ggplot2")
df <- data.frame(x = x_1)
p <- ggplot(data.frame(x = x), aes(x)) +
theme_bw() +
geom_histogram(bins = visr.param.num_histogram_bins) +
labs(x = "number of expressed genes", y = "number of cells")
print(p)


p <- ggplot(data.frame(x = x_1), aes(x)) +
  theme_bw() +
  geom_histogram(bins = visr.param.num_histogram_bins) +
  labs(x = "number of expressed genes, standardised to Z-distribution", y = "number of cells") +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = 'red')
print(p)

#
pData(my_cds)$UMI <- Matrix::colSums(Biobase::exprs(my_cds))
p <- ggplot(pData(my_cds), aes(num_genes_expressed, UMI)) +
    theme_bw() + geom_point() +
    ggtitle(paste("Number of genes expressed >=", visr.param.min_expr, "vs UMI per cell\nfor",nrow(pData(my_cds)), "cells"))
print(p)

#############################################
# Clustering
#############################################

if (visr.param.enable_clustering && ncol(my_cds) > 5e+05 && visr.param.cluster_method != CLUSTER_METHOD_LOUVAIN) {
  visr.message("Number of cells in your data is larger than 50k. clusterCells with 'densityPeak' or 'DDRTree' may crash. Please try to use the 'Louvain' clustering algorithm!", type = "warning")
}

if (visr.param.enable_subsetting) {
  # select gene based on their average expression and variability across cells.
  visr.logProgress("Calculating the mean and dispersion values for genes")
  if (visr.param.data_type != DATA_TYPE_UMI) {
    #TODO:
    stop("dispersionTable() only works for CellDataSet objects containing count-based expression data, either transcripts or reads.")
  }
  disp_table <- dispersionTable(my_cds)
  # head(disp_table)

  p <-  ggplot(data.frame(x = disp_table$mean_expression), aes(x)) +
      theme_bw() +
      geom_histogram(bins = 50) +
      scale_x_log10(labels = scales::comma) +
      labs(x = "mean gene expression", y = "number of genes") +
      geom_vline(xintercept = visr.param.min_mean_expression, linetype = "dotted", color = 'red')
  print(p)
  visr.writeDataTable(disp_table, paste0(visr.param.output_dir, "/genes_dispersion.txt"))

  print(table(disp_table$mean_expression >= visr.param.min_mean_expression))


  if (visr.param.gene_subset_method == GENE_SUBSET_METHOD_MEAN) {
    # subset genes based on minimum mean expression
    unsup_clustering_genes <- subset(disp_table, mean_expression >= visr.param.min_mean_expression)
    my_cds <- monocle::setOrderingFilter(my_cds, unsup_clustering_genes$gene_id)
  } else if (visr.param.gene_subset_method == GENE_SUBSET_METHOD_PERCENT) {
    # subset genes based on minimum number of expressed cells
    fData(my_cds)$use_for_ordering <- fData(my_cds)$num_cells_expressed > visr.param.min_expressed_cells * ncol(my_cds_subset)
  } else {
    visr.assert_that(FALSE, "Invalid value for gene_subset_method")
  }

  p <- monocle::plot_ordering_genes(my_cds) +
      ggtitle("Plot of genes by mean expression vs. empirical dispersion\nhighlighted genes (black) are selected for clustering") + theme(plot.title = element_text(hjust = 0.5)) +
      scale_x_log10(labels = scales::comma) +
      scale_y_log10(labels = scales::comma) +
      labs(x = "mean gene expression", y = "dispersion")
  print(p)
  visr.logProgress("Plotting the percentage of variance explained by each component (takes a few minutes)")

  p <-  monocle::plot_pc_variance_explained(my_cds, return_all = FALSE) + #TODO: norm_method= c("log", "vstExprs", "none")
    ggtitle("The percentage of variance explained by each component\nbased on a PCA performed on the normalised expression data") + theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept = visr.param.reduce_num_dim, linetype = "dotted", color = 'red')
  print(p)

  visr.logProgress("Computing a projection of a CellDataSet object into a lower dimensional space (takes a few minutes)")
  my_cds <- monocle::reduceDimension(my_cds, max_components = 2, num_dim = visr.param.reduce_num_dim,
                            reduction_method = visr.param.reduction_method, verbose = TRUE) # norm_method = c("log", "vstExprs", "none")
}

if (visr.param.enable_clustering) {
  num_clusters = NULL
  if (visr.param.cluster_method == CLUSTER_METHOD_DENSITY_PEAK) {
    # for some reason it alwasy produces one less cluster than specified value. so add 1 for now
    num_clusters = if (!is.null(visr.param.num_clusters)) (visr.param.num_clusters + 1) else NULL
  } else if (visr.param.cluster_method == CLUSTER_METHOD_DDRTREE) {
    num_clusters = visr.param.num_centers
  }

  visr.logProgress(paste("Performing unsupervised clustering. method:", visr.param.cluster_method))
  my_cds <- monocle::clusterCells(my_cds, verbose = T,
    skip_rho_sigma = visr.param.skip_rho_sigma,
    num_clusters = num_clusters,
    inspect_rho_sigma = F,
    rho_threshold = visr.param.rho_threshold,
    delta_threshold = visr.param.delta_threshold,
    gaussian = visr.param.gaussian,
    peaks = NULL,
    cell_type_hierarchy = NULL,
    frequency_thresh = NULL,
    enrichment_thresh = NULL,
    clustering_genes = NULL,
    k = visr.param.louvain_k,
    louvain_iter = visr.param.louvain_iter,
    weight = visr.param.louvain_weight,
    method = visr.param.cluster_method)

  visr.logProgress("Plotting the decision map of density clusters (delta vs. rho)")
  plot_rho_delta(my_cds, rho_threshold = visr.param.rho_threshold, delta_threshold = visr.param.delta_threshold) +
    ggtitle("Decision map of density clusters") + theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "rho: local density", y = "delta: local distance (to another cell with higher density)")

  head(pData(my_cds))
  print("Cluster sizes:")
  print(table(pData(my_cds)$Cluster))
  p <- monocle::plot_cell_clusters(my_cds)
  print(p)
}

#############################################
# Differential Expression Analysis
#############################################

if (visr.param.enable_de_analysis) {

  if (visr.param.marker_genes_list != '') {
    marker_genes <- row.names(subset(fData(HSMM_myo), gene_short_name %in% strsplit(visr.param.marker_genes_list, split = "[ ]*[,| ]+[ ]*")[[1]]))
    my_cds_subset <- my_cds[marker_genes,]
  }


  # create vector of no's
  de_cluster_ids = if (visr.param.de_cluster_id == Inf) levels(pData(my_cds)$Cluster) else visr.param.de_cluster_id

  i <- 0
  for (de_cluster_id in de_cluster_ids) {
    i <- i + 1
    # change status to yes if the cell was in cluster 1
    my_vector <- rep('no', nrow(pData(my_cds)))
    my_vector[pData(my_cds)$Cluster == de_cluster_id] <- rep('yes', sum(pData(my_cds)$Cluster == de_cluster_id))

    # add vector to phenoData
    pData(my_cds)$test <- my_vector

    head(pData(my_cds))


    visr.logProgress(paste("(", i, "of", length(de_cluster_ids), ")",
                           "Performing the differential expression analysis on cluster", de_cluster_id,
                           "\nfor the", nrow(unsup_clustering_genes),
                           "genes with the mean expression across cells >=", visr.param.min_mean_expression))

    de_cluster_one <- monocle::differentialGeneTest(my_cds[unsup_clustering_genes$gene_id,],
                                           fullModelFormulaStr = '~test',
                                           cores = 1, verbose = TRUE)

    #TODO:
    #clustering_DEG_genes <- monocle::differentialGeneTest(my_cds_subset,
    #                                                      fullModelFormulaStr = '~Cluster',
    #                                                      cores = 1, verbose = TRUE)

    dim(de_cluster_one)
    #de_cluster_one %>% arrange(qval) %>% head()
    #de_cluster_one_valid <- de_cluster_one %>% arrange(qval)
    de_cluster_one_valid <- de_cluster_one[which(de_cluster_one$status == "OK"),]
    de_cluster_one_valid <- de_cluster_one_valid[order(de_cluster_one_valid$qval),]

    de_genes_filename <- paste0(visr.param.output_dir, "/de_genes_cluster", de_cluster_id, ".txt")
    visr.logProgress(paste("Writing DE genes to file", de_genes_filename))
    visr.writeDataTable(de_cluster_one_valid, de_genes_filename)

    # the most statistically significant genes
    visr.logProgress(paste("Plotting level of expression for top", visr.param.num_plot_genes_jitter, "most statistically significant gene(s)."))
    monocle::plot_genes_jitter(my_cds[de_cluster_one_valid$id[seq(visr.param.num_plot_genes_jitter)],],
                               grouping = "Cluster", color_by = "Cluster",
                               ncol = ceiling(sqrt(visr.param.num_plot_genes_jitter)), nrow = NULL)
  }
}

#############################################
# Constructing Single Cell Trajectories
#############################################

#parameters, todo
trajectory_clusters <- c(1, 4)
trajectory_min_num_cells_expressed <- 10

if (visr.param.enable_trajectories) {
  #1- Choose genes that define progress
  #2- Reduce the dimensionality of the data
  #3- Order cells in pseudotime

  # expressed_genes <- row.names(subset(fData(my_cds), num_cells_expressed >= trajectory_min_num_cells_expressed))
  # my_cds_subset <- my_cds[expressed_genes, pData(my_cds)$Cluster %in% trajectory_clusters]
  # my_cds_subset <- monocle::detectGenes(my_cds_subset, min_expr = visr.param.min_expr)
  # fData(my_cds_subset)$use_for_ordering <- fData(my_cds_subset)$num_cells_expressed > visr.param.min_expressed_cells * ncol(my_cds_subset)
  # plot_pc_variance_explained(my_cds_subset, return_all = FALSE)
  # plot_rho_delta(my_cds_subset, rho_threshold = 2, delta_threshold = 10)
  my_cds_subset <- my_cds

  # TODO:
  # Alternative choices for ordering genes (http://cole-trapnell-lab.github.io/monocle-release/docs/#alternative-choices-for-ordering-genes)
  #   Ordering based on genes that differ between clusters (Recommended)
  #   Selecting genes with high dispersion across cells
  #   Ordering cells using known marker genes

  # how many genes are used?
  table(fData(my_cds_subset)$use_for_ordering)

  clustering_DEG_genes <- monocle::differentialGeneTest(my_cds_subset,
                                               fullModelFormulaStr = '~Cluster',
                                               cores = 1, verbose = TRUE)

  # use the top 1,000 most significantly differentially expressed genes as the set of ordering genes
  dim(clustering_DEG_genes)
  my_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
  my_cds_subset <- setOrderingFilter(my_cds_subset, ordering_genes = my_ordering_genes)
  my_cds_subset <- reduceDimension(my_cds_subset, method = 'DDRTree') # norm_method = c("log", "vstExprs", "none")
  # the warnings were for use of deprecated code
  my_cds_subset <- orderCells(my_cds_subset)
  plot_cell_trajectory(my_cds_subset, color_by = "Cluster")


  ######### Finding Genes that Change as a Function of Pseudotime
  my_pseudotime_de <- differentialGeneTest(my_cds_subset,
                                           fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                           cores = 1, verbose = TRUE)
  my_pseudotime_de %>% arrange(qval) %>% head()
  my_pseudotime_de %>% arrange(qval) %>% head() %>% select(id) -> my_pseudotime_gene
  my_pseudotime_gene <- my_pseudotime_gene$id
  plot_genes_in_pseudotime(my_cds_subset[my_pseudotime_gene,])


  ######### Clustering Genes by Pseudotemporal Expression Pattern
  # cluster the top 50 genes that vary as a function of pseudotime
  my_pseudotime_de %>% arrange(qval) %>% head(50) %>% select(id) -> gene_to_cluster
  gene_to_cluster <- gene_to_cluster$id

  # The columns of the heatmap are pseudotime values binned into 100 bins.
  my_pseudotime_cluster <- plot_pseudotime_heatmap(my_cds_subset[gene_to_cluster,],
                                                   num_clusters = 3,
                                                   cores = 1,
                                                   show_rownames = TRUE,
                                                   return_heatmap = TRUE)

  newdata <- data.frame(Pseudotime = seq(min(pData(my_cds_subset)$Pseudotime),
                                         max(pData(my_cds_subset)$Pseudotime), length.out = 100))

  # hierarchical clustering was used to cluster the genes
  # we can cut the dendrogram to form the same 3 clusters as plot_pseudotime_heatmap
  my_cluster <- cutree(my_pseudotime_cluster$tree_row, 3)
  my_cluster
  my_pseudotime_de[names(my_cluster[my_cluster == 1]),"gene_short_name"]
  my_pseudotime_de[names(my_cluster[my_cluster == 2]),"gene_short_name"]
  my_pseudotime_de[names(my_cluster[my_cluster == 3]),"gene_short_name"]


  ######### Analyzing Branches in Single-Cell Trajectories
  plot_cell_trajectory(my_cds_subset, color_by = "Cluster")

  # warnings not shown
  BEAM_res <- BEAM(my_cds_subset, branch_point = 1, cores = 8)
  BEAM_res <- BEAM_res[order(BEAM_res$qval),]
  BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

  # check out the results
  head(BEAM_res)
  table(BEAM_res$qval < 1e-4)

  # The heatmap shows how some genes are over-expressed or under-expressed depending on the trajectory path.
  my_branched_heatmap <- plot_genes_branched_heatmap(my_cds_subset[row.names(subset(BEAM_res, qval < 1e-4)),],
                                                     branch_point = 1,
                                                     num_clusters = 4,
                                                     cores = 1,
                                                     use_gene_short_name = TRUE,
                                                     show_rownames = TRUE,
                                                     return_heatmap = TRUE)


  head(my_branched_heatmap$annotation_row)
  dim(my_branched_heatmap$annotation_row)
  table(my_branched_heatmap$annotation_row$Cluster)
  my_row <- my_branched_heatmap$annotation_row
  my_row <- data.frame(cluster = my_row$Cluster,
                       gene = row.names(my_row),
                       stringsAsFactors = FALSE)

  head(my_row[my_row$cluster == 3,'gene'])

  my_gene <- row.names(subset(fData(my_cds_subset),
                              gene_short_name %in% head(my_row[my_row$cluster == 3,'gene'])))

  # plot genes that are expressed in a branch dependent manner
  # The trend lines show the expression pattern of genes along the paths formed by branch point 1.
  plot_genes_branched_pseudotime(my_cds_subset[my_gene,],
                                 branch_point = 1,
                                 ncol = 1)
}

#############################################
# Output tables
#############################################
dev.off(which=visr_app_pdf_dev)
rm(visr_app_pdf_dev)
browseURL(output_plot_file)

cds_dims <- data.frame(t(monocle::reducedDimA(my_cds)))
colnames(cds_dims) <- c("Component1", "Component2")
cell_table <- cbind(pData(my_cds), cds_dims)
visr.writeDataTable(cell_table, paste0(visr.param.output_dir, "/cells.txt"))

cachedResults <- list(params = visr.getParams(), my_cds = my_cds, disp_table = disp_table)
save(cachedResults, file = paste0(visr.param.output_dir, "/cachedResults.RData"))
