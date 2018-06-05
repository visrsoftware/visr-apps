source("visrutils.R")

visr.app.start("10X Cell Ranger",
               info = "Simple secondary analysis of the gene-barcode matrix output of 10X Cell Ranger pipeline",
               instructions = " - Specify the Cell Ranger data directory\n - Specify the desired workflow and parameters\n - Click [ Run ]"
               )

################################ Input
# pipeline output directory (pipestance path)
visr.app.category(label="Input",
                  info="10X Cell Ranger pipeline dataset directory. It should contain a sub-directory named 'outs'")

visr.param("data_dir_10x",
           label="10X dataset directory",
           info="10X Cell Ranger dataset directory. It should contain a sub-directory named 'outs'",
           type="filename", filename.mode = "dir",
           debugvalue= "~/Research/Data/SingleCell/pbmc3k/") # Peripheral Blood Mononuclear Cells (PBMCs) from a healthy donor

################################ Output tables

visr.category("Output",
              active.condition = "visr.param.data_dir_10x != ''")

visr.param("output_dir",
           label="Output Directory to save the results",
           info="Output directory where the analysis results will be saved to",
           type="filename", filename.mode = "dir",
           debugvalue= "~/Research/Data/SingleCell/output/")

visr.param("create_subdir", default = TRUE,
           label = "Create new sub-direcory",
           info = "Create a new sub directory with the name DATE_TIME (YYYYMMDD_hhmmss)")

################################ Analysis Steps

visr.category("Analysis steps",
              info = "Different analysis steps",
              active.condition = "visr.param.output_dir != '' && visr.param.data_dir_10x != ''")

visr.param("enable_clustering", label = "Output clustering results", default = T)

visr.param("enable_de", "Enable differential expression analysis", default = F, debugvalue = T)

################################ Precomputed analysis results

visr.app.category(label="Export precomputed clustering results",
                  active.condition = "visr.param.output_dir != '' && visr.param.data_dir_10x != '' && visr.param.enable_clustering",
                  info="The analysis table of cells together with precomputed dimensionality reduction and clustering results.")

visr.param("analysis_table", label="Name for table to output barcodes to", type="output-table",
           info="Name of the output table to appear in VisR",
           options = "importRowNames=FALSE",
           debugvalue = "cell_ranger_analysis")

visr.param("include_umi_counts", label="Include total UMI counts", type="boolean", default=T, debugvalue = T,
           info="Include a column with total number of UMIs for each cell in the output table?")

visr.param("include_tsne", label="Include t-SNE projection", type="boolean", default=T, debugvalue = T,
           info="Include t-SNE projections (2 columns) in the output table?")

visr.param("include_pca", label="Include PCA projection", type="boolean", default=T, debugvalue = T,
           info="Include PCA projections in the output table (10 columns: PCA.1 to PCA.10)?")

visr.param("include_kmeans", label="Include k-means results", type="boolean", default=T, debugvalue = T,
           info="Include precomputed kmeans clusterings in the output table (9 columns: kmeans.2 to kmeans.10)?")

visr.param("include_genecounts", label="Include counts for genes", type="boolean", default=FALSE, debugvalue = TRUE,
           info="Include gene counts for the genes specified below.")

visr.param("analysis_table_gene_list", label="Gene list to export (comma separated)", info="Comma separated list of genes to export to the output table",
           active.condition = "visr.param.include_genecounts",
           debugvalue = "CD79A, NKG7, CD3D, CST3, CD8A, PF4")

visr.param("gene_value_normalization", label="Normalization", info = "Normalization of gene values",
           items = c("norm_none", "norm_median_sum", "norm_log10_median_sum"),
           item.labels = c("None: raw values", "Sum over median sum", "log10 (Sum over median sum)"),
           default = "norm_log10_median_sum",
           active.condition = "visr.param.include_genecounts")

################################ Differential expression analysis

visr.app.category(label="Differential expression analysis",
                  active.condition = "visr.param.data_dir_10x != '' && visr.param.output_dir != '' && visr.param.enable_de",
                  info="The analysis table of genes together with the differential expression analysis results")

visr.param("de_cluster_k", label = "Use k-means clustering of k",
           default = 5L, min = 2L, max = 10L,
           active.condition = "visr.param.enable_de")

visr.param('de_only_significant', label="Export significant genes only", default=TRUE,
           info = "Only export genes that were significant in atleast one of the clusters. Uncheck to export all genes.")

visr.param('all_de_columns', label = "Include all columns in DE output table", default = T,
           info = "Whether to specify columns you want to be exported to the DE analysis output table. By default all columns are selected")

de_genes_fields = list(list('sum_a', 'Include sum of inside cluster values', T),
                       list('sum_b', 'Include sum of outside cluster values', T),
                       list('common_mean', 'Include common mean', T),
                       list('dispersion', 'Include dispersion', T),
                       list('mean_a_sizenorm', 'Include mean A size norm', T),
                       list('mean_b_sizenorm', 'Include mean B size norm', T),
                       list('log2fc', 'Include log2 fold change', T),
                       list('p', 'Include P', T),
                       list('p_adj', 'Include P (adjusted)', T),
                       list('order', 'Include significance order', T),
                       list('significant', 'Include significance condition', T))


for (f in de_genes_fields) {
  visr.param(paste('de_include_', f[[1]], sep=''), label=paste0('    ', f[[2]]), default=f[[3]],
             active.condition = "!visr.param.all_de_columns")
}

################################ Visualization Options

visr.app.category(label="Visualization Options",
                  active.condition = "visr.param.data_dir_10x != '' && visr.param.output_dir != ''",
                  info="Visualization of pre-computed results under a 2D t-SNE projection")

visr.param("visualize_tsne_kmeans", label="Plot t-SNE for ClusterIDs", default = T)

visr.param("visualize_tsne_umi", label="Plot t-SNE for total UMI count", default = T)

visr.param("vis_umi_min_value", type="double", label="min saturate value", info = "min saturate value on the color bar", default = 3,
           active.condition = "visr.param.visualize_tsne_umi")

visr.param("vis_umi_max_value", type="double", label="max saturate value", info = "max saturate value on the color bar", default = 4,
           active.condition = "visr.param.visualize_tsne_umi")

visr.param("vis_colormap_tsne_umi", label = "Color map", type = "multi-color", default="BuPu 7", debugvalue = "gray, red",
           active.condition = "visr.param.visualize_tsne_umi")

visr.param("visualize_tsne_markers", label="Plot t-SNE for gene markers", default = T)

visr.param("vis_gene_list", label="Gene list to visualize (comma separated)", info="Comma separated list of genes to visualize",
           active.condition = "visr.param.visualize_tsne_markers",
           debugvalue = "CD79A, NKG7, CD3D, CST3, CD8A, PF4")

visr.param("vis_min_value", type="double", label="min saturate value", info = "min saturate value on the color bar", default = 0.0,
           active.condition = "visr.param.visualize_tsne_markers")

visr.param("vis_max_value", type="double", label="max saturate value", info = "max saturate value on the color bar", default = 1.5,
           active.condition = "visr.param.visualize_tsne_markers")

visr.param("vis_colormap_tsne_markers", label = "Color map", type = "multi-color", default="BuPu 7", debugvalue = "gray, red",
           active.condition = "visr.param.visualize_tsne_markers")

visr.param("vis_show_de_heatmap", label = "show heatmap of significant genes",
           default = FALSE, debugvalue = TRUE,
           active.condition = "visr.param.enable_de")

visr.param("cluster_heatmap_n_genes", label = "number of genes per cluster", default = 3L, min = 1L, max = 100L,
           info = "Number of top differentially expressed genes to show per cluster in the heatmap",
           active.condition = "visr.param.enable_de && visr.param.vis_show_de_heatmap")

visr.param("cluster_heatmap_min", default = -1.0, label = "heatmap min value",
           info = "min saturate value on the color bar",
           active.condition = "visr.param.enable_de && visr.param.vis_show_de_heatmap")

visr.param("cluster_heatmap_max", default = 2.0, label= "heatmap max value",
           info = "max saturate value on the color bar",
           active.condition = "visr.param.enable_de && visr.param.vis_show_de_heatmap")


############################################################


#visr.param" barcode list
visr.app.end(printjson = TRUE, writefile = TRUE)

visr.applyParameters()

############################################################
############################################################
############################################################


############################################################
#  Validate directories
############################################################

validateOutputDirectory <- function() {
  visr.assert_that(!is.null(visr.param.output_dir) && visr.param.output_dir != "", msg = "Output directory not specified")
  if (visr.param.create_subdir)
    visr.param.output_dir <<- paste0(visr.param.output_dir, format(Sys.time(), "/%Y%m%d_%H%M%S"))
  dir.create(visr.param.output_dir, recursive = T, showWarnings = F)
  visr.assert_that(dir.exists(visr.param.output_dir), msg = sprintf("Output directory '%s' is not valid", visr.param.output_dir))
}

# validate input directory
visr.assert_file_exists(visr.param.data_dir_10x, '10X data directory')
visr.assert_that(file.exists(paste0(visr.param.data_dir_10x, '/outs')), msg = paste("Cannot find 'outs' sub-directory inside the specified 10x directory:\n", visr.param.data_dir_10x))

# validate output directory
validateOutputDirectory()

############################################################
#  Prepare output
############################################################

# output current parameters into a file
write(capture.output(print(visr.getParams())), file = paste0(visr.param.output_dir, "/visr_params.txt"))

# prepare the pdf to save the plots to
visr.dev <- dev.cur()
if (exists("visr_app_pdf_dev") && !is.null(visr_app_pdf_dev)) {
  dev.off(which = visr_app_pdf_dev) # close the previous one
}
output_plot_file <- paste0(visr.param.output_dir, "/plots.pdf")
pdf(file=output_plot_file)
visr_app_pdf_dev <- dev.cur()

############################################################
#  Load data
############################################################

# load cell-ranger library
visr.librarySource("cellrangerRkit", "http://cf.10xgenomics.com/supp/cell-exp/rkit-install-2.0.0.R")

# download sample data
# data_dir_10x <- "~/Research/Data/SingleCell/pbmc3k"
# download_sample(sample_name="pbmc3k",sample_dir=data_dir_10x, host="https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/")

visr.logProgress("Loading Cell Ranger matrix")
# gbm is based on the ExpressionSet class and stores the barcode filtered gene expression matrix and metadata (gene symbols and cell barcode IDs)
gbm <- cellrangerRkit::load_cellranger_matrix(visr.param.data_dir_10x)

gbm_summary <- cbind(rownames(t(gbm@summary)), t(gbm@summary))
colnames(gbm_summary) <- c("name", "value")
visr.writeDataTable(gbm_summary, file = paste0(visr.param.output_dir, "/gbm_summary.txt"))

visr.logProgress("Loading Cell Ranger analysis results")
analysis_results <- cellrangerRkit::load_cellranger_analysis_results(visr.param.data_dir_10x)

#
# Function to extract the list of valid gene names from a comma sepearted string
#
get_genes_from_string <- function(comma_string) {
    split_gene_list <- strsplit(comma_string, split = "[ ]*[,| ]+[ ]*")[[1]]
    if (length(which(split_gene_list != "")) == 0) {
        visr.message("No genes are specified for parameter 'Gene list'")
    }
    else
    {
        genes <- split_gene_list[which(split_gene_list != "")]
        suppressWarnings(indices <- sapply(genes, function(x) {return (cellrangerRkit:::get_gene_index(gbm, x))}))
        if (length(which(is.na(indices))) > 0) {
            visr.message(paste(c("Did not find and will ignore gene symbol(s):", genes[which(is.na(indices))]), collapse = " "), type = "warning")
        }
        genes <- genes[which(!is.na(indices))]
        if (length(genes) == 0) {
            visr.message("Could not find any of the specified gene names in the input data gene symbols")
        }
        return(genes)
    }
    return (c())
}

############################################################
#  Output table
############################################################

enable_analysis_table <- TRUE

if (enable_analysis_table) {
  analysis_table <- data.frame(analysis_results$tsne$Barcode)
  colnames(analysis_table) <- "Barcode"
}

## total UMI counts
if (enable_analysis_table && visr.param.include_umi_counts) {
  umi_counts <- data.frame(colSums(exprs(gbm)))
  umi_counts <- cbind(umi_counts, log(umi_counts[,1], base = 10))
  colnames(umi_counts) <- c("UMI.counts", "Log10.UMI.counts")
  analysis_table <- cbind(analysis_table, umi_counts)
}

## t-SNE results
if (enable_analysis_table && visr.param.include_tsne) {
  analysis_table <- cbind(analysis_table, analysis_results$tsne[,-1])
}

## PCA results
if (enable_analysis_table && visr.param.include_pca) {
  analysis_table <- cbind(analysis_table, analysis_results$pca[,-1])
}

## k-means results
if (enable_analysis_table && visr.param.include_kmeans) {
  n_clu <- 2:10
  km_res <- analysis_results$clustering # load pre-computed kmeans results
  clu_res <- sapply(n_clu, function(x) as.character(km_res[[paste("kmeans",x,"clusters",sep="_")]]$Cluster))
  colnames(clu_res) <- sapply(n_clu, function(x) paste("kmeans",x,sep="."))
  analysis_table <- cbind(analysis_table, clu_res)
}

## gene counts
if (enable_analysis_table && visr.param.include_genecounts) {
    genes <- get_genes_from_string(visr.param.analysis_table_gene_list)
    if (length(genes) > 0) {
        gbm_norm <- gbm
        if (visr.param.gene_value_normalization == "norm_median_sum"  || visr.param.gene_value_normalization == "norm_log10_median_sum") {
            use_genes <- cellrangerRkit::get_nonzero_genes(gbm)
            gbm_norm <- cellrangerRkit::normalize_barcode_sums_to_median(gbm[use_genes,]) # Normalize barcodes by sum to the median barcode by sum
            if (visr.param.gene_value_normalization == "norm_log10_median_sum") {
                gbm_norm <- log_gene_bc_matrix(gbm_norm, base=10)
            }
        }
        gbm_trunc <- cellrangerRkit::trunc_gbm_by_genes(gbm_norm, genes)
        gene_values <- t(as.matrix(exprs(gbm_trunc)))
        colnames(gene_values) <- genes
        analysis_table <- cbind(analysis_table, gene_values)
    }
}

if (enable_analysis_table) {
  row.names(analysis_table) <- NULL
  visr.writeDataTable(analysis_table, paste0(visr.param.output_dir, "/cells.txt"))
}

############################################################
#  Visualization
############################################################

## visualize umi counts
if (visr.param.visualize_tsne_umi) {
  p <- visualize_umi_counts(gbm, analysis_results$tsne[c("TSNE.1","TSNE.2")],
                            limits=c(visr.param.vis_umi_min_value, visr.param.vis_umi_max_value),
                            marker_size=0.05)
  p <- p + visr.util.scale_color_gradient(visr.param.vis_colormap_tsne_umi, label="log10(UMI_counts)")
  print(p)
}

## visualize clusters
if (visr.param.visualize_tsne_kmeans) {
  n_clu <- 2:10
  km_res <- analysis_results$clustering # load pre-computed kmeans results
  clu_res <- sapply(n_clu, function(x) km_res[[paste("kmeans",x,"clusters",sep="_")]]$Cluster)
  colnames(clu_res) <- sapply(n_clu, function(x) paste("kmeans",x,sep="."))
  p <- cellrangerRkit::visualize_clusters(clu_res, analysis_results$tsne[c("TSNE.1","TSNE.2")])
  print(p)
}

## visualize gene markers
if (visr.param.visualize_tsne_markers) {
  genes <- get_genes_from_string(visr.param.vis_gene_list)
  if (length(genes) > 0) {
    use_genes <- get_nonzero_genes(gbm)
    gbm_bcnorm <- cellrangerRkit::normalize_barcode_sums_to_median(gbm[use_genes,])
    gbm_log <- cellrangerRkit::log_gene_bc_matrix(gbm_bcnorm, base=10)
    print(dim(gbm_log))
    p <- cellrangerRkit::visualize_gene_markers(gbm_log, genes,
                                analysis_results$tsne[c("TSNE.1","TSNE.2")],
                                limits = c(visr.param.vis_min_value, visr.param.vis_max_value))
    p <- p + visr.util.scale_color_gradient(visr.param.vis_colormap_tsne_markers, label="normalized\nsums")
    print(p)
  }
}



############################################################
#  Differential Expression Analysis
############################################################
if (visr.param.enable_de) {

  visr.logProgress("Computing differential expression parameters")

  cluster_result <- analysis_results$clustering[[paste("kmeans", visr.param.de_cluster_k,"clusters",sep="_")]]

  # sort the cells by the cluster labels
  cells_to_plot <- cellrangerRkit::order_cell_by_clusters(gbm, cluster_result$Cluster)

  # order the genes from most up-regulated to most down-regulated in each cluster
  prioritized_genes <- cellrangerRkit::prioritize_top_genes(gbm, cluster_result$Cluster, method = "sseq",
                                                            p_cutoff = 0.05,
                                                            min_log2fc = 0,
                                                            min_mean = 1)
  # TODO: expose these parameters?
  # \item{method}{Method to prioritize genes: currently only supports mean difference ('mean_diff', 'unique', or 'sseq')}
  # \item{logscale}{Logical if the input is log-scale or not (default: TRUE)}
  # \item{p_cutoff}{If method is "sseq," only consider genes w/ adjusted p-value below this value}
  # \item{min_mean}{If method is "sseq," only consider genes w/ mean normalized UMI counts/cell exceeding this value}
  # \item{min_log2fc}{If method is "sseq," only consider genes with at least this log2 fold-change}
  # \item{order_by}{If method is "sseq," sort genes by 'pvalue' or by 'log2fc'}
  # \item{n_top_genes}{Number of top genes for each cluster group to consider}

  if (visr.param.vis_show_de_heatmap) {
    visr.logProgress("Rendering heatmap of DE genes")

    brewer_pal_name = "Set2"
    if (visr.param.de_cluster_k > 8) {
      brewer_pal_name = "Set3"
    }
    example_col <- rev(brewer.pal(visr.param.de_cluster_k, brewer_pal_name)) # customize plotting colors
    # scatter plot
    # visualize_clusters(cluster_result$Cluster, tsne_proj[c("TSNE.1","TSNE.2")], colour=example_col)

    # create values and axis annotations for pheatmap
    cellrangerRkit::gbm_pheatmap(log_gene_bc_matrix(gbm),
                 prioritized_genes,
                 cells_to_plot,
                 n_genes=visr.param.cluster_heatmap_n_genes,
                 colour=example_col,
                 limits=c(visr.param.cluster_heatmap_min, visr.param.cluster_heatmap_max))
  }

  visr.logProgress("Exporting DE table")

  data_table_genes <- prioritized_genes[[1]][order(prioritized_genes[[1]]$ix),c("gene_id", "gene_name", "tested")]
  data_table_genes$significant_in_any = FALSE
  data_table_genes$significant_in_cluster = ''
  data_table_genes_colnames <- c("gene_id", "gene_name", "tested", "significant_in_any", 'significant_in_cluster')

  for (c in seq(visr.param.de_cluster_k)) {
    genes_order <- order(prioritized_genes[[c]]$ix)
    genes_df = prioritized_genes[[c]][genes_order,]
    genes_df$order <- genes_order
    stopifnot(all.equal(genes_df$gene_id, data_table_genes$gene_id))

    sig_label <- data.frame(new_label = rep("", nrow(genes_df)), comma='', stringsAsFactors=FALSE)
    sig_label$new_label[which(genes_df$significant)] = as.character(c)
    sig_label$comma[which(genes_df$significant & data_table_genes$significant_in_any)] = ','
    data_table_genes$significant_in_cluster = paste(data_table_genes$significant_in_cluster, sig_label$comma, sig_label$new_label, sep = '')

    data_table_genes$significant_in_any = data_table_genes$significant_in_any | genes_df$significant

    for (field in de_genes_fields) {
      field.name <- field[[1]]
      if (eval(parse(t=paste('visr.param.de_include_', field.name, sep = '')))) {
        data_table_genes <- cbind(data_table_genes, genes_df[, field.name])
        data_table_genes_colnames <- c(data_table_genes_colnames, paste("c", c, "_", field.name, sep=""))
      }
    }
  }
  colnames(data_table_genes) <- data_table_genes_colnames
  rownames(data_table_genes) <- NULL
  de_analysis_table <- data_table_genes
  if (visr.param.de_only_significant) {
    de_analysis_table <- data_table_genes[which(data_table_genes$significant_in_any), ]
  }
  visr.writeDataTable(de_analysis_table, paste0(visr.param.output_dir, "/de_analysis.txt"))
}

#############################################
# Output tables
#############################################
dev.off(which=visr_app_pdf_dev)
rm(visr_app_pdf_dev)
browseURL(output_plot_file)

cachedResults <- list(params = visr.getParams())
save(cachedResults, file = paste0(visr.param.output_dir, "/cachedResults.RData"))

rm(gbm)
