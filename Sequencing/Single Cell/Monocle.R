#TODO:
# - Save and load CellDataSet object
# - save cells.txt and genes.txt, de_genes, BEAM_res
# - subset cells based on barcode
# - add summaries to the report
# - classify cells based on gene-markers: http://cole-trapnell-lab.github.io/monocle-release/docs/#classifying-cells-by-type-recommended
# - cluster using marker genes: http://cole-trapnell-lab.github.io/monocle-release/docs/#clustering-cells-using-marker-genes-recommended
#   Myoblast: [MYF5] >= 1, Fibroblast: [MYF5] < 1 & [ANPEP] > 1
# - ordering genes using marker genes: http://cole-trapnell-lab.github.io/monocle-release/docs/#alternative-choices-for-ordering-genes
# add filtering / tSNE / clustering results to the cells.txt

source("visrutils.R")

curr_dir <- "Sequencing/Single Cell"

visr.app.start("Monocle")

source(sprintf("%s/Utils.R",curr_dir))
source(sprintf("%s/Monocle_IO.R",curr_dir))
source(sprintf("%s/Monocle_steps.R",curr_dir))
source(sprintf("%s/Monocle_filter.R",curr_dir))
source(sprintf("%s/Monocle_dim_reduction.R",curr_dir))
source(sprintf("%s/Monocle_cluster_cells.R",curr_dir))
source(sprintf("%s/Monocle_DE.R",curr_dir))
source(sprintf("%s/Monocle_construct_trajectories.R",curr_dir))
source(sprintf("%s/Monocle_pseudotime_genes.R",curr_dir))
source(sprintf("%s/Monocle_analyze_branches.R",curr_dir))
source(sprintf("%s/Monocle_additional.R",curr_dir))

visr.app.end()

visr.applyParameters()

visr.var.add_title_pages <- TRUE

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

#  Preparing output directory
output_dir <- prepare_output()

# Loading input
monocle_app_object <- load_input()

# Filtering: remove low (e.g. dead cells or empty wells) and high (e.g. doublets: made from two or more cells accidentally)
if (visr.param.filter_by_distribution) {
  monocle_app_object <- perform_filter_by_distribution(monocle_app_object)
  monocle_app_object <- perform_detect_genes(monocle_app_object)
}

# Here is the earliest place the clustering parameters could be validated
num_cells <- ncol(monocle_app_object$cds)
if (visr.param.enable_clustering && num_cells > 5e+05 && visr.param.cluster_method != CLUSTER_METHOD_LOUVAIN) {
  visr.message("Number of cells in your data is larger than 50k. clusterCells with 'densityPeak' or 'DDRTree' may crash. Please try to use the 'Louvain' clustering algorithm!", type = "warning")
}

# Subsetting genes
if (visr.param.enable_filtering) {
  monocle_app_object <- perform_subsetting_genes(monocle_app_object)
}

# dimensionality reduction
if (visr.param.enable_dim_red) {
  monocle_app_object <- perform_dim_reduction(monocle_app_object)
}

# Clustering
if (visr.param.enable_clustering) {
  monocle_app_object <- perform_clustering(monocle_app_object)
}

# Differential Expression Analysis
if (visr.param.enable_de_analysis) {
  monocle_app_object <- perform_de_analysis(monocle_app_object, output_dir)
}

# Construct Single Cell Trajectories
if (visr.param.enable_trajectories) {
  monocle_app_object <- perform_construct_trajectories(monocle_app_object)
}

# Cluster genes along pseudotime
if (visr.param.find_pseudotime_genes) {
  monocle_app_object <- perform_pseudotime_genes(monocle_app_object)
}

# Analyze trajectory branches
if (visr.param.analyze_trajectory_branches) {
  monocle_app_object <- perform_analyze_branches(monocle_app_object)
}

# Output tables
finalize_output(output_dir, monocle_app_object)
