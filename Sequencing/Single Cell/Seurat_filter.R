create_seurat_cond <- sprintf("visr.param.Output_directory != '' && %s", load_raw_valid)

visr.app.category("Create Seurat Object", active.condition = create_seurat_cond)
visr.param("max_nGenes", label = "Maximum number of genes per cell", type = "int", min = 0, debugvalue = 2500,items = c("Inf"), default = "Inf",
           info = "High cutoff for the number of genes expressed in a cell")
visr.param("min_nGenes", label = "Minimum number of genes per cell", type = "int", min=0,default = 200, debugvalue = 200,
           info = "Low cutoff for the number of genes expressed in a cell")
visr.param("max_percent_mito", label = "Maximum fraction of mitocondrial genes", min = 0, max = 1,default = 0.05,debugvalue = 0.05,
           info = "High cutoff for the proportion of UMIs from mitochondrial genes; Enter a value between 0 and 1")
visr.param("min_fraction_cells", label = "Minimum fraction of cells per gene", min = 0, max = 1, default = 0.001, debugvalue = 0.001,
           info = "Low cutoff for the fraction of cells that a gene is expressed")

merge_obj <- function(gbmData1,gbmData2){
  if (length(gbmData1@var.genes) == 0){
    gbmData1 <- find_variable_genes(gbmData1)
  }
  if (length(gbmData2@var.genes) == 0){
    gbmData2 <- find_variable_genes(gbmData2)
  }
  print("Merging objects")
  gbmData1@meta.data$group <- visr.param.dataset1_name
  #gbmData1 <- ScaleData(gbmData1)
  
  gbmData2@meta.data$group <- visr.param.dataset2_name
  #gbmData2 <- ScaleData(gbmData2)
  
  genes.use <- unique(c(gbmData1@var.genes, gbmData2@var.genes))
  genes.use <- intersect(genes.use, rownames(gbmData1@data))
  genes.use <- intersect(genes.use, rownames(gbmData2@data))
  
  combined.object <- MergeSeurat(object1 = gbmData1, object2 = gbmData2, 
                                 do.scale = FALSE, do.center = FALSE)
  #combined.object <- ScaleData(object = combined.object)
  #combined.object@scale.data[is.na(x = combined.object@scale.data)] <- 0
  combined.object@var.genes <- genes.use
  
  return(combined.object)
}

plot_meta_data <- function(gbmData, filtered){
  if (!filtered){
    title1 <- "Cell Meta Data (before filtering)"
    title2 <- "Meta Data Comparison (before filtering)"
  }else{
    title1 <- "Cell Meta Data (after filtering)"
    title2 <- "Meta Data Comparison (after filtering)"
  }
  p <- VlnPlot(object = gbmData, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3,size.title.use = 15)
  p <- p + ggtitle(title1) + theme(plot.title = element_text(lineheight=2,size = 20), plot.margin = margin(20, 10, 10, 10))
  print(p)
  
  par(mfrow = c(1, 2), oma = c(0, 0, 3, 0))
  GenePlot(object = gbmData, gene1 = "nUMI", gene2 = "percent.mito")
  GenePlot(object = gbmData, gene1 = "nUMI", gene2 = "nGene")
  mtext(title2, outer = TRUE, cex = 1.5)
}

filter_Cells <- function(gbmData){
  print(paste("Filtering cells"))
  max_number_of_genes <- visr.param.max_nGenes
  max_fraction_of_mito <- visr.param.max_percent_mito
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = gbmData@data), value = T, ignore.case = T)
  percent.mito <- Matrix::colSums(gbmData@raw.data[mito.genes, ])/Matrix::colSums(gbmData@raw.data)
  
  gbmData <- AddMetaData(object = gbmData, metadata = percent.mito, col.name = "percent.mito")
  
  
  plot_meta_data(gbmData,filtered = F)
  
  gbmData <- FilterCells(object = gbmData, subset.names = c("nGene", "percent.mito"), low.thresholds = c(min_number_of_genes, -Inf), high.thresholds = c(max_number_of_genes, max_fraction_of_mito))
  
  plot_meta_data(gbmData, filtered = T)
  switchPlotToScreen()
  plot_meta_data(gbmData, filtered = T)
  switchPlotToReport()
  
  return(gbmData)
}

load_data <- function(){
  print(paste("Loading data"))
  if (visr.param.workflow == "single"){
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
      gbmData <- Read10X(path)
      print(paste("Creating Seurat Object"))
      num_cells <- dim(gbmData)[2]
      gbmData <- CreateSeuratObject(raw.data = gbmData, min.cells = max(1,ceiling(num_cells*min_fraction_of_cells)), min.genes = visr.param.min_nGenes, 
                                    project = project_name)
      gbmData <- filter_Cells(gbmData)
      
    }else{
      # load Seurat object directly
      gbmData <- readRDS(visr.param.Path_to_seurat_object)
      if (is.null(gbmData@meta.data$percent.mito)){
        # mito.genes <- grep(pattern = "^MT-", x = rownames(x = gbmData@data), value = T, ignore.case = T)
        # percent.mito <- Matrix::colSums(gbmData@raw.data[mito.genes, ])/Matrix::colSums(gbmData@raw.data)
        # gbmData <- AddMetaData(object = gbmData, metadata = percent.mito, col.name = "percent.mito")
        gbmData <- filter_Cells(gbmData)
      }
    }
  }else{
    if (visr.param.Import_method2 == "load_one"){
      gbmData <- readRDS(visr.param.seurat_integrated)
      visr.assert_that((length(gbmData@var.genes) > 0), msg = "Run find variable genes on original datasets")
    }else{
      gbmData1 <- readRDS(visr.param.seurat_obj_1)
      gbmData2 <- readRDS(visr.param.seurat_obj_2)
      gbmData <- merge_obj(gbmData1,gbmData2)
      rm(gbmData1,gbmData2)
    }
  }
  return(gbmData)
}