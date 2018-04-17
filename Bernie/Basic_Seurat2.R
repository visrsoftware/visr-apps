source("visrutils.R")

visr.app.start("Basic_Seurat2", debugdata=mtcars)

#parameters

#input
visr.app.category("Input")
visr.param("Import_method",items = c("load_raw","load_seurat"),
           item.labels = c("Load from cellRanger outs","Load Seurat Object"),default = "load_raw", 
           debugvalue = "load_seurat")
visr.param("Path_to_outs", type="filename",filename.mode = "dir", 
           debugvalue = "C:/Users/Yiwei Zhao/Desktop/BRC/tools/meta/single/",
           active.condition = "visr.param.Import_method == 'load_raw'")

visr.param("Path_to_seurat_object", type="filename",filename.mode = "load",
           #debugvalue = "C:/Users/Yiwei Zhao/Desktop/BRC/tools/meta/single/sample.Robj",
           debugvalue = "C:/Users/Yiwei Zhao/Desktop/DN88DF.Robj",
           active.condition = "visr.param.Import_method == 'load_seurat'")

visr.param("Output_plots", type="filename",filename.mode = "save", 
           debugvalue = "C:/Users/Yiwei Zhao/Desktop/BRC/tools/meta/single/visr_seurat.pdf")

visr.param("Save_output_object", default = T, debugvalue = F)
visr.param("Output_object", type="filename",filename.mode = "save", 
           #debugvalue = "C:/Users/Yiwei Zhao/Desktop/BRC/tools/meta/single/sample.Robj",
           debugvalue = "C:/Users/Yiwei Zhao/Desktop/DN88DF.Robj",
           active.condition = "visr.param.Save_output_object == true")

#options
visr.param("Filter_Cells", default = F, debugvalue = F)
visr.param("Find_Variable_Genes", default = F, debugvalue = F)
visr.param("Run_PCA", default = F, debugvalue = F)
visr.param("elbow", label = "Draw Elbow Plot", default = F, debugvalue = F)
visr.param("jackstraw", label = "Run Jackstraw", default = F, debugvalue = F)
visr.param("jackstrawRep",label = "number of replicates", default = 100,
            active.condition = "visr.param.jackstraw == true")
visr.param("Run_tSNE", default = F, debugvalue = F)
visr.param("Cluster_Cells", default = F, debugvalue = T)

#filter cells
visr.app.category("Filter cells", active.condition = "visr.param.Filter_Cells == true")
visr.param("max_nGenes",default = 2500, debugvalue = 2500)
visr.param("min_nGenes",default = 200, debugvalue = 200)
visr.param("max_percent_mito",default = 0.05,debugvalue = 0.05)

#Find variable genes
visr.app.category("Variable gene detection", active.condition = "visr.param.Find_Variable_Genes == true")
visr.param("Normalization_scale_factor", default = 10000,
           debugvalue = 10000)
visr.param("Mean_exp_low",default = 0.0125,debugvalue = 0.0125)
visr.param("Mean_exp_high",default = 3.0, debugvalue = 3.0)
visr.param("Dispersion_low",default = 0.5, debugvalue = 0.5)
#visr.param("Dispersion_high")

#tsne
visr.app.category("Run tSNE", active.condition = "visr.param.Run_tSNE == true")
visr.param("specify_tsne_nPC", label = "Manually specify PC number", default = F, debug = T)
visr.param("tsne_nPC", label = "Number of PCs", default = 10, 
           active.condition = "visr.param.specify_tsne_nPC == true")

#cluster cells
visr.app.category("Cluster Cells", active.condition = "visr.param.Cluster_Cells == true")
visr.param("specify_cluster_nPC", label = "Manually specify PC number", default = F, debug = T)
visr.param("cluster_nPC", label = "Number of PCs", default = 10, 
           active.condition = "visr.param.specify_cluster_nPC == true")
visr.param("cluster_resolution", label = "resolution", default = 0.6, debugvalue = 0.6) #increase for large dataset

visr.app.end(printjson=TRUE, writefile=TRUE)
visr.applyParameters()

###
### check path here

library(Seurat)
library(dplyr)
library(Matrix)


###

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

# save plots to this location
pdf(file=visr.param.Output_plots)

# filter out cells with excessive mitochondrial genes
if (visr.param.Filter_Cells){
  print(paste(project_name,":","Filtering cells"))
  max_number_of_genes <- visr.param.max_nGenes
  max_fraction_of_mito <- visr.param.max_percent_mito
  
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = gbmData@data), value = TRUE)
  percent.mito <- Matrix::colSums(gbmData@raw.data[mito.genes, ])/Matrix::colSums(gbmData@raw.data)
  gbmData <- AddMetaData(object = gbmData, metadata = percent.mito, col.name = "percent.mito")
  print(VlnPlot(object = gbmData, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3))
  # par(mfrow = c(1, 2))
  # GenePlot(object = gbmData, gene1 = "nUMI", gene2 = "percent.mito")
  # GenePlot(object = gbmData, gene1 = "nUMI", gene2 = "nGene")
  
  gbmData <- FilterCells(object = gbmData, subset.names = c("nGene", "percent.mito"), low.thresholds = c(min_number_of_genes, -Inf), high.thresholds = c(max_number_of_genes, max_fraction_of_mito))
  
  # Normalization
  print(paste(project_name,":","Performing normalization"))
  normalization_scale_factor <- visr.param.Normalization_scale_factor
  gbmData <- NormalizeData(object = gbmData, normalization.method = "LogNormalize", scale.factor = normalization_scale_factor)
}

# Find variable genes
if (visr.param.Find_Variable_Genes){
  print(paste(project_name,":","Finding variable genes"))
  x.low.cutoff <- visr.param.Mean_exp_low
  x.high.cutoff <- visr.param.Mean_exp_high
  y.cutoff <- visr.param.Dispersion_low
  
  gbmData <- FindVariableGenes(object = gbmData, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = x.low.cutoff, x.high.cutoff = x.high.cutoff, y.cutoff = y.cutoff)
  print(paste(project_name,":",paste("Number of variable genes: ", toString(length(x = gbmData@var.genes)), sep = "")))
}

# PCA
if (visr.param.Run_PCA) {
  # Scale data (regression)
  print(paste(project_name,":","Scaling data"))
  gbmData <- ScaleData(object = gbmData, vars.to.regress = c("nUMI", "percent.mito"))
  print(paste(project_name,":","Running PCA"))
  gbmData <- RunPCA(object = gbmData, pc.genes = gbmData@var.genes,do.print=F)
}

#elbow plot
if (visr.param.elbow){
  print(paste(project_name,":","Drawing elbow plot"))
  print(PCElbowPlot(object = gbmData)) 
}

# Jackstraw
if (visr.param.jackstraw){
  print(paste(project_name,":","Running JackStraw"))
  gbmData <- JackStraw(object = gbmData, num.replicate = visr.param.jackstrawRep, do.print = T)
  print(JackStrawPlot(object = gbmData, PCs = 1:20))
}

# Select cutoff for elbow automatically
if (!(is.null(gbmData@dr$pca))){
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
}

# Cluster cells
if (visr.param.Cluster_Cells){
  if (visr.param.specify_cluster_nPC){
    num_pc_to_use <- visr.param.cluster_nPC
  } else {
    num_pc_to_use <- auto_num_pc_to_use
  }
  print(paste(project_name,":","Clustering cells"))
  gbmData <- FindClusters(object = gbmData, reduction.type = "pca", dims.use = 1:num_pc_to_use, resolution = visr.param.cluster_resolution, print.output = 0, save.SNN = TRUE)
  # PrintFindClustersParams(object = gbmData)
}

# tSNE
if (visr.param.Run_tSNE) {
  if (visr.param.specify_tsne_nPC){
    num_pc_to_use <- visr.param.tsne_nPC
  } else {
    num_pc_to_use <- auto_num_pc_to_use
  }
  print(paste(project_name,":","Running tSNE"))
  gbmData <- RunTSNE(object = gbmData, dims.use = 1:num_pc_to_use, do.fast = TRUE)
  TSNEPlot(object = gbmData)
}


#
dev.off()

if (visr.param.Save_output_object){
  print(paste(project_name,":","Saving object"))
  saveRDS(gbmData, file = visr.param.Output_object) 
}