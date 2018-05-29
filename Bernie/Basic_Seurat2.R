source("visrutils.R")

# Parameters --------------------------------------------------------------

visr.app.start("Basic_Seurat2", debugdata=mtcars, input.type = "none")

op <- "("
cp <- ")"

#input
visr.app.category("Input")
visr.param("Import_method", label = "Choose Import Method",
           items = c("load_raw","load_seurat"),
           item.labels = c("Load CellRanger Output","Load Seurat Object"),default = "load_raw", 
           debugvalue = "load_raw")
visr.param("Path_to_outs", type="filename",filename.mode = "dir", 
           info="Cell Ranger pipeline output directory. It should contain another directory named \"outs\"",
           debugvalue = "C:/Users/Yiwei Zhao/Desktop/BRC/tools/meta/single/",
           active.condition = "visr.param.Import_method == 'load_raw'")

visr.param("Path_to_seurat_object", type="filename",filename.mode = "load",
           info = "Path to the Seurat object file.",
           #debugvalue = "C:/Users/Yiwei Zhao/Desktop/BRC/tools/meta/single/sample.Robj",
           debugvalue = "C:/Users/Yiwei Zhao/Desktop/DN88DF.raw.Robj",
           active.condition = "visr.param.Import_method == 'load_seurat'")

load_raw_valid <- "(visr.param.Import_method == 'load_raw' && visr.param.Path_to_outs != '')"
load_seurat_valid <- "(visr.param.Import_method == 'load_seurat' && visr.param.Path_to_seurat_object != '')"
load_data_valid <- paste(c(op,load_raw_valid, "||", load_seurat_valid,cp), collapse = "")

#output
visr.app.category("Output", active.condition = load_data_valid)
visr.param("Output_directory", type = "filename", filename.mode="dir",
           label="Output directory to save the results",
           info="Output directory where the analysis results will be saved to",
           debugvalue="C:/Users/Yiwei Zhao/Desktop/")
visr.param("create_subdir", default = T,
           label = "Create new sub-direcory",
           info = "Create a new sub directory with the name DATE_TIME (YYYYMMDD_hhmmss)")

#Create Seurat Object and filter
create_seurat_cond <- paste(c(op,load_raw_valid,"&&",op,"visr.param.Output_directory != ''",cp,cp),collapse = "")
visr.app.category("Create Seurat Object", active.condition = create_seurat_cond)
visr.param("max_nGenes", type = "int", min = 0, debugvalue = 2500,items = c("Inf"), default = "Inf",
           info = "High cutoff for the number of genes expressed in a cell")
visr.param("min_nGenes",min=0,default = 0, debugvalue = 200,
           info = "Low cutoff for the number of genes expressed in a cell")
visr.param("max_percent_mito", label = "Max fraction Mito", min = 0, max = 1,default = 0.05,debugvalue = 0.05,
           info = "High cutoff for the proportion of UMIs from mitochondrial genes; Enter a value between 0 and 1")
visr.param("min_fraction_cells", min = 0, max = 1, default = 0.001, debugvalue = 0.001,
           info = "Low cutoff for the fraction of cells that a gene is expressed")

#options
analysis_cond <-  paste(c(op,load_data_valid,"&&",op,"visr.param.Output_directory != ''",cp,cp),collapse = "")
visr.app.category("Analysis steps",info = "Select the analysis steps to run",
                  active.condition = analysis_cond)
# visr.param("Filter_Cells", default = F, debugvalue = F)
visr.param("Find_Variable_Genes", default = F, debugvalue = F)
visr.param("Dim_Reduction", label = "Dimensionality Reduction", default = F, debugvalue = F)
visr.param("Cluster_Cells", default = F, debugvalue = F)
visr.param("Find_Marker_Genes", default = F, debug = F)

#Find variable genes
find_var_cond <- paste(c(op,analysis_cond,"&&",op,"visr.param.Find_Variable_Genes == T",cp,cp),collapse = "")
visr.app.category("Variable gene detection", active.condition = find_var_cond)
visr.param("Normalization_scale_factor", min = 1, default = 10000, debugvalue = 10000,
           info = "Sets the scale factor for cell-level normalization")
visr.param("Mean_exp_low",min = 0, default = 0.0125,debugvalue = 0.0125,
           info = "Low cutoff on x-axis (mean expression) for identifying variable genes")
visr.param("Mean_exp_high", type = "int", min = 0, default = 3.0, debugvalue = 3.0, items = c("Inf"),
           info = "High cutoff on x-axis (mean expression) for identifying variable genes")
visr.param("Dispersion_low", min = 0, default = 0.5, debugvalue = 0.5,
           info = "Low cutoff on y-axis (standard deviation) for identifying variable genes")
visr.param("Dispersion_high", type = "int", min = 0, debugvalue = Inf, items = c("Inf"),default="Inf",
           info = "High cutoff on y-axis (standard deviation) for identifying variable genes")

#Dim Reduction
dim_red_cond <- paste(c(op,analysis_cond,"&&",op,"visr.param.Dim_Reduction == T",cp,cp),collapse = "")
visr.app.category("Dimensionality Reduction",active.condition = dim_red_cond)
visr.param("Run_PCA", default = F, debugvalue = F,
           info = "Run a PCA dimensionality reduction")
visr.param("nPC_compute", min = 1, default = 20, debugvalue = 20,label = "Number of PCs to compute",
           active.condition = "visr.param.Run_PCA == T",
           info = "Total Number of PCs to compute and store")
visr.param("jackstraw", label = "Run Jackstraw", default = F, debugvalue = F,
           info = "Randomly permutes a subset of data, and calculates projected PCA scores for these 'random' genes. Then compares the PCA scores for the 'random' genes with the observed PCA scores to determine statistical signifance. End result is a p-value for each gene's association with each principal component.")
visr.param("jackstrawRep",label = "Number of replicates", min = 10, default = 100,
           active.condition = "visr.param.jackstraw == T",
           info = "Number of replicate samplings to perform")
visr.param("nPC_jackstrawPlot", label = "Number of PCs to plot", default = 12,
           active.condition = "visr.param.jackstraw == T",
           info = "Number of PCs to include on jackstraw plot")
visr.param("elbow", label = "Draw Elbow Plot", default = F, debugvalue = F,
           info = "Plots the standard deviation of each principle component")

visr.param("PC_heatmap","Plot PC Heatmap", default = F,
           info = "Draws a heatmap focusing on each principal component. Both cells and genes are sorted by their principal component scores. Allows for nice visualization of sources of heterogeneity in the dataset.")
visr.param("nPC_PCheatmap", label = "Number of PCs to plot", default = 12,
           active.condition = "visr.param.PC_heatmap == T",
           info = "Number of PCs to include on PC heatmap")
visr.param("nGene_PCHeatmap", label = "Number of Genes to plot", default = 30,
           info = "Number of top genes to plot for each PC")
visr.param("nCell_PCHeatmap", label = "Number of Cells to plot", default = 500,
           active.condition = "visr.param.PC_heatmap == T",
           info = "Number of top cells to plot for each PC")

visr.param("Run_tSNE", default = F, debug = F, active.condition = "visr.param.Dim_Reduction == T",
           info = "Run t-SNE dimensionality reduction on selected PCs")
visr.param("calculate_tsne_nPC", label = "Automatically calculate number of PCs", default = T, debug = F,
           active.condition = "visr.param.Run_tSNE == T",
           info = "Automatically calculate number of PCs selected for tSNE. Uncheck to specify the number of PCs to use.")
visr.param("tsne_nPC", label = "Number of PCs for calculating tSNE", min = 1, default = 10, 
           active.condition = "visr.param.calculate_tsne_nPC == F")


#cluster cells
cluster_cond <- paste(c(op,analysis_cond,"&&",op,"visr.param.Cluster_Cells == T",cp,cp),collapse = "")
visr.app.category("Cluster Cells", active.condition = cluster_cond)
visr.param("calculate_cluster_nPC", label = "Automatically calculate number of PCs", default = T, debug = F,
           info = "Automatically calculate number of PCs selected for clustring cells. Uncheck to specify the number of PCs to use.")
visr.param("cluster_nPC", label = "Number of PCs for clustering", min = 1, default = 10,
           active.condition = "visr.param.calculate_cluster_nPC == F")
visr.param("cluster_resolution", label = "resolution", min = 0.1, default = 0.6, debugvalue = 0.6,
           info = "Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.") #increase for large dataset

#Differential expression
DE_cond <- paste(c(op,analysis_cond,"&&",op,"visr.param.Find_Marker_Genes == T",cp,cp),collapse = "")
visr.app.category("Differential Expression",active.condition = DE_cond)
visr.param("Choose_Clusters",
           items=c("all","selected"), 
           item.labels = c("Compare each cluster to the rest of the cells","Select sepcific group of clusters"),
           default = "all")
visr.param("group_1", type = "character", default = "0,1", active.condition = "visr.param.Choose_Clusters == 'selected'",
           info = "Group of clusters to define markers for")
visr.param("group_2", type = "character", default = "2,3", active.condition = "visr.param.Choose_Clusters == 'selected'",
           info = "Group of clusters for comparison. If empty, assume all clusters that are not in group 1.")
visr.param("DE_test_method", items = c("wilcox","bimod","roc","t","tobit","poisson","negbinom"), default = "wilcox",
           info = "Denotes which test to use")
visr.param("min_pct",label = "Min percent of cells", default = 0.1, 
           info = "Only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations")
visr.param("min_logfc", label = "LogFC threshold", default = 0.25,
           info = "Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells.")

visr.param("Top_gene_n",label = "Number of top genes to export", default = 2, min = 1, items = c(Inf), 
           item.labels = c("All genes"),info = "Number of top genes for each target group in the output file")


#Additional Analysis results
visr.app.category("Additional Output", active.condition = analysis_cond)
visr.param("export_results",label = "Export analysis results", default = F)
visr.param("include_id",label="Include cluster id", default = T, active.condition = "visr.param.export_results == T")
visr.param("include_umi",label = "Include total UMI", default = F, active.condition = "visr.param.export_results == T")
visr.param("include_pc",label="Include PCA projection", default = F,  active.condition = "visr.param.export_results == T")
visr.param("nPC_export",label="Number of PCs to export", default = 2, min = 1, 
           active.condition = "visr.param.export_results == T && visr.param.include_pc")
visr.param("include_tsne",label="Include t-SNE projection", default = F,  active.condition = "visr.param.export_results == T")
visr.param("plot_genes",label = "Visualize selected genes", default = F)
visr.param("gene_name", label = "Gene list ot visualize (comma separated)",
           debugvalue = "MS4A1", active.condition = "visr.param.plot_genes == T")
visr.param("gene_plot_color", label = "Color map", type = "multi-color", default="BuPu 7", debugvalue = "gray, red",
           active.condition = "visr.param.plot_genes == T")

visr.app.end(printjson=T, writefile=T)
visr.applyParameters()


# Check Path --------------------------------------------------------------
if (visr.param.Import_method == "load_raw"){
  if (!dir.exists(visr.param.Path_to_outs)){visr.message("Path to outs directory is not found")}
} else {
  if (!file.exists(visr.param.Path_to_seurat_object)){visr.message("Path to seurat object is not found")}
}

if (!dir.exists(visr.param.Output_directory)){visr.message("Path to output directory not found")}

# Load Packages -----------------------------------------------------------

library(Seurat)
library(dplyr)
library(reshape2)

# Functions ---------------------------------------------------------------

filter_Cells <- function(gbmData){
  print(paste(project_name,":","Filtering cells"))
  max_number_of_genes <- visr.param.max_nGenes
  max_fraction_of_mito <- visr.param.max_percent_mito
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = gbmData@data), value = T, ignore.case = T)
  percent.mito <- Matrix::colSums(gbmData@raw.data[mito.genes, ])/Matrix::colSums(gbmData@raw.data)
  gbmData <- AddMetaData(object = gbmData, metadata = percent.mito, col.name = "percent.mito")
  
  print(VlnPlot(object = gbmData, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3))
  
  par(mfrow = c(1, 2))
  GenePlot(object = gbmData, gene1 = "nUMI", gene2 = "percent.mito")
  GenePlot(object = gbmData, gene1 = "nUMI", gene2 = "nGene")
  
  gbmData <- FilterCells(object = gbmData, subset.names = c("nGene", "percent.mito"), low.thresholds = c(min_number_of_genes, -Inf), high.thresholds = c(max_number_of_genes, max_fraction_of_mito))
  return(gbmData)
}

find_variable_genes <- function(gbmData){
  if (is.null(gbmData@meta.data$percent.mito)){
    gbmData <- filter_Cells(gbmData)
  }
  # Normalization
  print(paste(project_name,":","Performing normalization"))
  normalization_scale_factor <- visr.param.Normalization_scale_factor
  gbmData <- NormalizeData(object = gbmData, normalization.method = "LogNormalize", scale.factor = normalization_scale_factor)
  
  # Find variable genes
  print(paste(project_name,":","Finding variable genes"))
  x.low.cutoff <- visr.param.Mean_exp_low
  x.high.cutoff <- visr.param.Mean_exp_high
  y.cutoff <- visr.param.Dispersion_low
  y.high.cutoff <- visr.param.Dispersion_high
  gbmData <- FindVariableGenes(object= gbmData, mean.function = ExpMean, dispersion.function = LogVMR, 
                               x.low.cutoff = x.low.cutoff, x.high.cutoff = x.high.cutoff, y.cutoff = y.cutoff, y.high.cutoff = y.high.cutoff)
  print(paste(project_name,":",paste("Number of variable genes: ", toString(length(x = gbmData@var.genes)), sep = "")))
  if (length(gbmData@var.genes) == 0){
    visr.message("No variable genes found")
  }
  return(gbmData)
}

run_PCA <- function(gbmData){
  if (length(gbmData@var.genes) == 0){
    gbmData <- find_variable_genes(gbmData)
  }
  
  # Scale data (regression)
  print(paste(project_name,":","Scaling data"))
  gbmData <- ScaleData(object = gbmData, vars.to.regress = c("nUMI", "percent.mito"))
  
  print(paste(project_name,":","Running PCA"))
  gbmData <- RunPCA(object = gbmData, pc.genes = gbmData@var.genes, do.print=F, 
                    pcs.compute = visr.param.nPC_compute)
  return(gbmData)
}

draw_elbow <- function(gbmData){
  if (is.null(gbmData@dr$pca)){
    gbmData <- run_PCA(gbmData)
  }
  print(paste(project_name,":","Drawing elbow plot"))
  print(PCElbowPlot(object = gbmData))
  return(gbmData)
}

run_jackstraw <- function(gbmData){
  if (is.null(gbmData@dr$pca)){
    gbmData <- run_PCA(gbmData)
  }
  print(paste(project_name,":","Running JackStraw"))
  gbmData <- JackStraw(object = gbmData, num.replicate = visr.param.jackstrawRep, do.print = T)
  nPC <- min(visr.param.nPC_jackstrawPlot,length(gbmData@dr$pca@sdev))
  print(JackStrawPlot(object = gbmData, PCs = 1:nPC))
  return(gbmData)
}

draw_PCHeatmap <- function(gbmData){
  if (is.null(gbmData@dr$pca)){
    gbmData <- run_PCA(gbmData)
  }
  nPC <- min(visr.param.nPC_PCheatmap,length(gbmData@dr$pca@sdev))
  nCell <- min(visr.param.nCell_PCHeatmap,nrow(gbmData@dr$pca@cell.embeddings))
  nGene <- min(visr.param.nGene_PCHeatmap,nrow(gbmData@dr$pca@gene.loadings))
  print(PCHeatmap(object = gbmData, pc.use = 1:nPC, 
                  cells.use = nCell,
                  num.genes = nGene,
                  do.balanced = T, label.columns = F, use.full = F))
  return(gbmData)
}

calculate_nPC <- function(gbmData){
  variance <- gbmData@dr$pca@sdev^2
  mean_var <- sum(variance)/length(variance)
  percent.var <- variance/sum(variance)
  min_change <- 0.001
  for (i in 2:length(percent.var)){
    #print(paste(project_name,":",c(i,variance[i],sum(percent.var[1:i])))
    if ((variance[i]<mean_var) & ((percent.var[i-1]-percent.var[i]) < min_change)){
      break
    }
  }
  auto_num_pc_to_use <- i
  print(paste(project_name,": auto_num_pc = ",auto_num_pc_to_use))
  return(auto_num_pc_to_use)
}

cluster_cells <- function(gbmData){
  if (is.null(gbmData@dr$pca)){
    gbmData <- run_PCA(gbmData)
  }
  if (visr.param.calculate_cluster_nPC){
    num_pc_to_use <- calculate_nPC(gbmData)
  }else{
    num_pc_to_use <- min(visr.param.cluster_nPC,length(gbmData@dr$pca@gene.loadings.full))
  }
  print(num_pc_to_use)
  print(paste(project_name,":","Clustering cells"))
  gbmData <- FindClusters(object = gbmData, reduction.type = "pca", dims.use = 1:num_pc_to_use, 
                          resolution = visr.param.cluster_resolution, print.output = 0, save.SNN = T, 
                          temp.file.location = visr.param.Output_directory)
  # PrintFindClustersParams(object = gbmData)
  if (length(levels(gbmData@ident)) == 1){
    visr.message("Cells cannot be clustered")
  }
  return(gbmData)
}

run_tSNE <- function(gbmData){
  if (is.null(gbmData@dr$pca)){
    gbmData <- run_PCA(gbmData)
  }
  
  if (visr.param.calculate_tsne_nPC){
    num_pc_to_use <- calculate_nPC(gbmData)
  }else{
    num_pc_to_use <- visr.param.tsne_nPC
  }
  print(num_pc_to_use)
  print(paste(project_name,":","Running tSNE"))
  gbmData <- RunTSNE(object = gbmData, dims.use = 1:num_pc_to_use, do.fast = T)
  TSNEPlot(object = gbmData)
  return(gbmData)
}

plot_genes <- function (gbm, gene_probes, projection, limits = c(0, 10), marker_size = 0.1, title = NULL) 
{
  gene_values <- t(as.matrix(rbind(numeric(0),gbm@data[gene_probes,])))
  gene_values[gene_values < limits[1]] <- limits[1]
  gene_values[gene_values > limits[2]] <- limits[2]
  colnames(gene_values) <- gene_probes
  projection_names <- colnames(projection)
  colnames(projection) <- c("Component.1", "Component.2")
  proj_gene <- data.frame(cbind(projection, gene_values))
  proj_gene_melt <- melt(proj_gene, id.vars = c("Component.1", 
                                                "Component.2"))
  p <- ggplot(proj_gene_melt, aes(Component.1, Component.2)) + 
    geom_point(aes(colour = value), size = marker_size) + 
    facet_wrap(~variable) + scale_colour_gradient(low = "grey", 
                                                  high = "red", name = "val") + labs(x = projection_names[1], 
                                                                                     y = projection_names[2])
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  p <- p + theme_bw() + theme(plot.title = element_text(hjust = 0.5), 
                              panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(p)
}


# Analysis result file names ----------------------------------------------

plot_output <- "plots.pdf"
cell_info_output <- "Cells.tsv"
DE_output <- "Differential_Expression_Analysis_result.tsv"
selected_gene_output <- "Selected_Marker_Genes.tsv"

# Create PDF --------------------------------------------------------------

# shut off exisiting device
visr.dev <- dev.cur()
if (exists("seurat_app_pdf_dev") && !is.null(seurat_app_pdf_dev)){dev.off(which = seurat_app_pdf_dev)}

# save plots to this location
if (visr.param.create_subdir){
  subfolder <- Sys.time()
  subfolder <- gsub(" ","_",subfolder)
  subfolder <- gsub(":","-",subfolder)
  subfolder <- gsub("-","",subfolder)
  output_folder <- paste(visr.param.Output_directory,subfolder,sep = "/")
} else {
  output_folder <- visr.param.Output_directory
}
dir.create(output_folder)

pdf(file = paste(output_folder,plot_output,sep="/"))
seurat_app_pdf_dev <- dev.cur()
plot.new()

# Load Data ---------------------------------------------------------------

min_fraction_of_cells <- visr.param.min_fraction_cells
min_number_of_genes <- visr.param.min_nGenes

if (visr.param.Import_method == "load_raw"){
  #Input data path
  path_to_outs <- visr.param.Path_to_outs
  path <- paste(path_to_outs,"outs/filtered_gene_bc_matrices",sep="/")
  #search for genome
  genome <- list.files(path=path)
  path <- paste(path,genome,sep = "/")
  project_name <- tail(strsplit(path_to_outs,split="/")[[1]],n=1) #10x_pbmc
  print(project_name)
  
  # load from outs folder
  print(paste(project_name,":","Loading data"))
  gbmData <- Read10X(path)
  print(paste(project_name,":","Creating Seurat Object"))
  num_cells <- dim(gbmData)[2]
  gbmData <- CreateSeuratObject(raw.data = gbmData, min.cells = max(1,ceiling(num_cells*min_fraction_of_cells)), min.genes = visr.param.min_nGenes, 
                                project = project_name)
  gbmData <- filter_Cells(gbmData)
  
}else{
  # load Seurate object directly
  print("Loading Seurat Object")
  gbmData <- readRDS(visr.param.Path_to_seurat_object)
  project_name <- gbmData@project.name
}

# Analysis ----------------------------------------------------------------

# Find variable genes
if (visr.param.Find_Variable_Genes){
  gbmData <- find_variable_genes(gbmData)
}

if (visr.param.Dim_Reduction){
  #PCA
  if (visr.param.Run_PCA){
    gbmData <- run_PCA(gbmData)
  }
  
  #elbow plot
  if (visr.param.elbow){
    gbmData <- draw_elbow(gbmData)
    
    #display plot on screen
    dev.set(which = visr.dev)
    gbmData <- draw_elbow(gbmData)
    dev.set(which = seurat_app_pdf_dev)
  }
  
  # Jackstraw
  if (visr.param.jackstraw){
    gbmData <- run_jackstraw(gbmData)
  }
  
  # PC heatmap
  if (visr.param.PC_heatmap){
    gbmData <- draw_PCHeatmap(gbmData)
  }
}

# Cluster cells
if (visr.param.Cluster_Cells){
  gbmData <- cluster_cells(gbmData)
}

# tSNE
if (visr.param.Run_tSNE && visr.param.Dim_Reduction) {
  gbmData <- run_tSNE(gbmData) 
}

#DE analysis
if (visr.param.Find_Marker_Genes){
  if (length(levels(gbmData@ident)) == 1){
    gbmData <- cluster_cells(gbmData)
  }
  top_genes_n <- min(visr.param.Top_gene_n,nrow(gbmData@raw.data))
  print(paste(project_name,":","Running DE analysis"))
  if (visr.param.Choose_Clusters == "selected"){
    all.clusters <- as.numeric(levels(gbmData@ident))
    group_1 <- eval(parse(text=paste("c(",visr.param.group_1,")")))
    if (is.null(group_1)){visr.message("Enter clusters of interest")}
    
    group_2 <- eval(parse(text=paste("c(",visr.param.group_2,")")))
    group_1.markers <- FindMarkers(object = gbmData, ident.1 = group_1, ident.2 = group_2, 
                                   test.use = visr.param.DE_test_method,min.pct = visr.param.min_pct,
                                   logfc.threshold = visr.param.min_logfc)
    print(head(group_1.markers,n=5))
    table <- data.frame(head(group_1.markers,n=top_genes_n))
    gene <- rownames(table)
    gene <- data.frame(gene)
    table <- cbind(gene,table)
      
  }else{
    DE.markers <- FindAllMarkers(object = gbmData, test.use = visr.param.DE_test_method,
                                 min.pct = visr.param.min_pct, logfc.threshold = visr.param.min_logfc)
    table <- data.frame(DE.markers %>% group_by(cluster) %>% top_n(top_genes_n, avg_logFC))
    table <- table[,c("gene","cluster","avg_logFC","pct.1","pct.2","p_val","p_val_adj")]
  }
  write.table(x = table, file = paste(output_folder, DE_output,sep = "/"), row.names = F, quote = F, sep = "\t")
}

# Additional output -------------------------------------------------------

# Cell Info
if (visr.param.export_results){
  barcodes <- data.frame(gbmData@cell.names)
  colnames(barcodes) <- "Barcode"
  table <- barcodes
  if (visr.param.include_id){
    if (length(levels(gbmData@ident)) == 1){
      gbmData <- cluster_cells(gbmData)
    }
    clusterID <- data.frame(gbmData@ident)
    colnames(clusterID) <- "Cluster ID"
    table <- cbind(table,clusterID)
  }
  
  if (visr.param.include_umi){
    umi_counts <- data.frame(colSums(gbmData@data))
    colnames(umi_counts) <- "Total UMI"
    table <- cbind(table,umi_counts)
  }
  
  if (visr.param.include_pc){
    if (is.null(gbmData@dr$pca)){
      gbmData <- run_PCA(gbmData)
    }
    PCs <- data.frame(gbmData@dr$pca@cell.embeddings)
    nPC <- min(visr.param.nPC_export,ncol(PCs))
    colnames(PCs) <- colnames(gbmData@dr$pca@cell.embeddings)
    PCs <- PCs[,1:nPC]
    table <- cbind(table,PCs)
  }
  
  if (visr.param.include_tsne){
    if (is.null(gbmData@dr$tsne)){
      gbmData <- run_tSNE(gbmData)
    }
    tSNEs <- gbmData@dr$tsne@cell.embeddings
    colnames(tSNEs) <- colnames(gbmData@dr$tsne@cell.embeddings)
    table <- cbind(table,tSNEs)
  }
  
  row.names(table) <- NULL
  write.table(x = table, file = paste(output_folder,cell_info_output,sep = "/"),
              quote = F, row.names = F, sep = "\t")
}

# Plot genes

if (visr.param.plot_genes){
  gene_probes <- strsplit(gsub(visr.param.gene_name, pattern = ' ', replacement = ''),split = ',')[[1]]
  if (length(gene_probes) > 0){
    p <- plot_genes(gbm = gbmData,gene_probes = gene_probes,projection = gbmData@dr$tsne@cell.embeddings)
    p <- p + visr.util.scale_color_gradient(visr.param.gene_plot_color, label="log10(UMI_counts)")
    dev.set(which = visr.dev)
    print(p)
    dev.set(which = seurat_app_pdf_dev)
    print(p)
    
    gene_values <- data.frame(t(as.matrix(rbind(numeric(0),gbmData@data[gene_probes,]))))
    colnames(gene_values) <- gene_probes
    gene_values$cluster <- gbmData@ident
    table <- gene_values %>% group_by(cluster) %>% summarise_at(.vars = gene_probes,.funs = mean)
    table <- t(table)
    write.table(table,file = paste(output_folder,selected_gene_output,sep = "/"), row.names = T, col.names = F, quote = F, sep = "\t")
  }
}

# Save Object -------------------------------------------------------------

# close current device
dev.off(which=seurat_app_pdf_dev)
rm(seurat_app_pdf_dev)

print(paste(project_name,":","Saving object"))
saveRDS(gbmData, file = paste(output_folder,"Seurat.Robj",sep = "/"))

#rm(gbmData)
#gc()

#export parameters

#visr.getParams()
#visr.assertthat()
#visualize marker -> output marker/additional output (plot and export table)
        


