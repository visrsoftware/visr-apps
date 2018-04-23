source("visrutils.R")

visr.app.start("Basic_Seurat2", debugdata=mtcars, input.type = "none")

###parameters

#input
visr.app.category("Input")
visr.param("Import_method",items = c("load_raw","load_seurat"),
           item.labels = c("Load from cellRanger outs","Load Seurat Object"),default = "load_raw", 
           debugvalue = "load_raw")
visr.param("Path_to_outs", type="filename",filename.mode = "dir", 
           debugvalue = "C:/Users/Yiwei Zhao/Desktop/BRC/tools/meta/single/",
           active.condition = "visr.param.Import_method == 'load_raw'")

visr.param("Path_to_seurat_object", type="filename",filename.mode = "load",
           #debugvalue = "C:/Users/Yiwei Zhao/Desktop/BRC/tools/meta/single/sample.Robj",
           debugvalue = "C:/Users/Yiwei Zhao/Desktop/DN88DF.raw.Robj",
           active.condition = "visr.param.Import_method == 'load_seurat'")

visr.app.category("Output")
visr.param("Output_plots", type="filename",filename.mode = "save", 
           debugvalue = "C:/Users/Yiwei Zhao/Desktop/visr_seurat.pdf")

visr.param("Save_output_object", default = T, debugvalue = F)
visr.param("Output_object", type="filename",filename.mode = "save", 
           #debugvalue = "C:/Users/Yiwei Zhao/Desktop/BRC/tools/meta/single/sample.Robj",
           debugvalue = "C:/Users/Yiwei Zhao/Desktop/DN88DF.delete.Robj",
           active.condition = "visr.param.Save_output_object == true")

#options
visr.app.category("Options")
visr.param("Filter_Cells", default = F, debugvalue = F)
visr.param("Find_Variable_Genes", default = F, debugvalue = F)
visr.param("Run_PCA", default = F, debugvalue = F)
visr.param("nPC_compute", min = 1, default = 20, debugvalue = 20,label = "Number of PC to compute",
           active.condition = "visr.param.Run_PCA == true")
visr.param("elbow", label = "Draw Elbow Plot", default = F, debugvalue = F)
visr.param("jackstraw", label = "Run Jackstraw", default = F, debugvalue = F)
visr.param("jackstrawRep",label = "number of replicates", min = 1, default = 100,
            active.condition = "visr.param.jackstraw == true")
visr.param("Cluster_Cells", default = F, debugvalue = F)
visr.param("Run_tSNE", default = F, debugvalue = F)

#filter cells
visr.app.category("Filter cells", active.condition = "visr.param.Filter_Cells == true")
visr.param("max_nGenes",min = 0, default = 2500, debugvalue = 2500,items = c("Inf"))
visr.param("min_nGenes",min=0,default = 200, debugvalue = 200)
visr.param("max_percent_mito",min = 0, max = 1,default = 0.05,debugvalue = 0.05)

#Find variable genes
visr.app.category("Variable gene detection", active.condition = "visr.param.Find_Variable_Genes == true")
visr.param("Normalization_scale_factor", min = 1, default = 10000, debugvalue = 10000)
visr.param("Mean_exp_low",min = 0, default = 0.0125,debugvalue = 0.0125)
visr.param("Mean_exp_high",min = 0, default = 3.0, debugvalue = 3.0)
visr.param("Dispersion_low", min = 0, default = 0.5, debugvalue = 0.5)
visr.param("Dispersion_high",min = 0, default = 100, debugvalue = Inf, items = c(Inf))

#cluster cells
visr.app.category("Cluster Cells", active.condition = "visr.param.Cluster_Cells == true")
visr.param("calculate_cluster_nPC", label = "Automatically calculate number of PCs", default = T, debug = F)
visr.param("cluster_nPC", label = "Number of PCs", min = 1, default = 10,
           active.condition = "visr.param.calculate_cluster_nPC == false")
visr.param("cluster_resolution", label = "resolution", min = 0.1, default = 0.6, debugvalue = 0.6) #increase for large dataset

#tsne
visr.app.category("Run tSNE", active.condition = "visr.param.Run_tSNE == true")
visr.param("calculate_tsne_nPC", label = "Automatically calculate number of PCs", default = T, debug = F)
visr.param("tsne_nPC", label = "Number of PCs", min = 1, default = 10, 
           active.condition = "visr.param.calculate_tsne_nPC == false")

#Differential expression
visr.app.category("Differential Expression")
visr.param("Find_Marker_Genes", default = F, debug = F)
visr.param("Choose_Clusters",
           items=c("all","selected"), 
           item.labels = c("Compare each group to the rest of the cells","Select sepcific group of clusters"),
           default = "all", active.condition = "visr.param.Find_Marker_Genes == true")
visr.param("group_1", type = "character", default = "0,1",
           active.condition = "visr.param.Choose_Clusters == 'selected'")
visr.param("group_2", type = "character", default = "2,3",
           active.condition = "visr.param.Choose_Clusters == 'selected'")

visr.app.end(printjson=TRUE, writefile=TRUE)
visr.applyParameters()


### check path here


plot_output_dir <- paste( head(strsplit(visr.param.Output_plots,split="/")[[1]],-1), collapse = "/")
object_output_dir <- paste( head(strsplit(visr.param.Output_object,split="/")[[1]],-1), collapse = "/")

if (visr.param.Import_method == "load_raw"){
  if (!dir.exists(visr.param.Path_to_outs)){stop("Path to outs directory is not found")}
} else {
  if (!file.exists(visr.param.Path_to_seurat_object)){stop("Path to seurat object is not found")}
}
if (!dir.exists(plot_output_dir)){print(plot_output_dir);stop("Path to plot output directory is not found")}
if (visr.param.Save_output_object){
  if (!dir.exists(object_output_dir)){stop("Path to object output directory is not found")}
}

### 

library(Seurat)
library(dplyr)
library(Matrix)

### functions

filter_Cells <- function(gbmData){
  print(paste(project_name,":","Filtering cells"))
  max_number_of_genes <- visr.param.max_nGenes
  max_fraction_of_mito <- visr.param.max_percent_mito
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = gbmData@data), value = TRUE)
  percent.mito <- Matrix::colSums(gbmData@raw.data[mito.genes, ])/Matrix::colSums(gbmData@raw.data)
  gbmData <- AddMetaData(object = gbmData, metadata = percent.mito, col.name = "percent.mito")
  
  print(VlnPlot(object = gbmData, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3))
  
  par(mfrow = c(1, 2))
  GenePlot(object = gbmData, gene1 = "nUMI", gene2 = "percent.mito")
  GenePlot(object = gbmData, gene1 = "nUMI", gene2 = "nGene")
  
  gbmData <- FilterCells(object = gbmData, subset.names = c("nGene", "percent.mito"), low.thresholds = c(min_number_of_genes, -Inf), high.thresholds = c(max_number_of_genes, max_fraction_of_mito))
  
  # Normalization
  print(paste(project_name,":","Performing normalization"))
  normalization_scale_factor <- visr.param.Normalization_scale_factor
  gbmData <- NormalizeData(object = gbmData, normalization.method = "LogNormalize", scale.factor = normalization_scale_factor)
  return(gbmData)
}

find_variable_genes <- function(gbmData){
  if (is.null(gbmData@meta.data$percent.mito)){
    gbmData <- filter_Cells(gbmData)
  }
  print(paste(project_name,":","Finding variable genes"))
  x.low.cutoff <- visr.param.Mean_exp_low
  x.high.cutoff <- visr.param.Mean_exp_high
  y.cutoff <- visr.param.Dispersion_low
  
  gbmData <- FindVariableGenes(object= gbmData, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = x.low.cutoff, x.high.cutoff = x.high.cutoff, y.cutoff = y.cutoff)
  print(paste(project_name,":",paste("Number of variable genes: ", toString(length(x = gbmData@var.genes)), sep = "")))
  if (length(gbmData@var.genes) == 0){
    stop("No variable genes found")
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
    num_pc_to_use <- visr.param.cluster_nPC
  }
  print(num_pc_to_use)
  print(paste(project_name,":","Clustering cells"))
  gbmData <- FindClusters(object = gbmData, reduction.type = "pca", dims.use = 1:num_pc_to_use, 
                          resolution = visr.param.cluster_resolution, print.output = 0, save.SNN = TRUE, 
                          temp.file.location = visr.param.Output_object)
  # PrintFindClustersParams(object = gbmData)
  if (length(levels(gbmData@ident)) == 1){
    stop("Cells cannot be clustered")
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
  gbmData <- RunTSNE(object = gbmData, dims.use = 1:num_pc_to_use, do.fast = TRUE)
  TSNEPlot(object = gbmData)
  return(gbmData)
}

### load data

min_fraction_of_cells <- 0.001
min_number_of_genes <- visr.param.min_nGenes

if (visr.param.Import_method == "load_raw"){
  #Input data path
  path_to_outs <- visr.param.Path_to_outs #single
  path <- paste(path_to_outs,"outs/filtered_gene_bc_matrices",sep="/")
  #search for genome
  genome <- list.files(path=path)
  path <- paste(path,genome,sep = "/")
  print(path)
  project_name <- tail(strsplit(path_to_outs,split="/")[[1]],n=1) #10x_pbmc
  print(project_name)
  
  # load from outs folder
  print(paste(project_name,":","Loading data"))
  gbmData <- Read10X(path)
  
  # create Seurat object
  print(paste(project_name,":","Creating Seurat Object"))
  num_cells <- dim(gbmData)[2]
  gbmData <- CreateSeuratObject(raw.data = gbmData, min.cells = max(1,ceiling(num_cells*min_fraction_of_cells)), min.genes = min_number_of_genes, project = project_name)

  }else{
  # load Seurate object directly
  print("Loading Seurat Object")
  gbmData <- readRDS(visr.param.Path_to_seurat_object)
  project_name <- gbmData@project.name
}

# shut off exisiting device
if (exists("seurat_app_pdf_dev") && !is.null(seurat_app_pdf_dev)){dev.off(which = seurat_app_pdf_dev)}

# save plots to this location
pdf(file=visr.param.Output_plots)
seurat_app_pdf_dev <- dev.cur()

### Analysis

# filter out cells with excessive mitochondrial genes
if (visr.param.Filter_Cells){
  gbmData <- filter_Cells(gbmData)
}

# Find variable genes
if (visr.param.Find_Variable_Genes){
  gbmData <- find_variable_genes(gbmData)
}

#PCA
if (visr.param.Run_PCA){
  gbmData <- run_PCA(gbmData)
}

#elbow plot
if (visr.param.elbow){
  gbmData <- draw_elbow()
}

# Jackstraw
if (visr.param.jackstraw){
  print(paste(project_name,":","Running JackStraw"))
  gbmData <- JackStraw(object = gbmData, num.replicate = visr.param.jackstrawRep, do.print = T)
  print(JackStrawPlot(object = gbmData, PCs = 1:20))
}

# Cluster cells
if (visr.param.Cluster_Cells){
  gbmData <- cluster_cells(gbmData)
}

# tSNE
if (visr.param.Run_tSNE) {
  gbmData <- run_tSNE(gbmData)
}

#DE analysis
#visr.param.group_1
if (visr.param.Find_Marker_Genes){
  if (length(levels(gbmData@ident)) == 1){
    gbmData <- cluster_cells(gbmData)
  }
  if (visr.param.Choose_Clusters == "selected"){
    all.clusters <- as.numeric(levels(gbmData@ident))
    group_1 <- eval(parse(text=paste("c(",visr.param.group_1,")")))
    if (is.null(visr.param.group_2)){
      group_2 <- (all.clusters+1)[group_1+1]-1
    }else{
      group_2 <- eval(parse(text=paste("c(",visr.param.group_2,")")))
    }
    group_1.markers <- FindMarkers(object = gbmData, ident.1 = group_1, ident.2 = group_2)
    print(head(group_1.markers,n=5))
  }else{
    DE.markers <- FindAllMarkers(object = gbmData, min.pct = 0.25)
    print(DE.markers %>% group_by(cluster) %>% top_n(2, avg_logFC))
  }
}

# close current device
dev.off(which=seurat_app_pdf_dev)
rm(seurat_app_pdf_dev)


# save object
if (visr.param.Save_output_object){
  print(paste(project_name,":","Saving object"))
  saveRDS(gbmData, file = visr.param.Output_object)
}
