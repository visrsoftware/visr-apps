#Dim Reduction
dim_red_cond <- sprintf("visr.param.Dim_Reduction == T && %s",analysis_cond)
dim_red_cond2 <- sprintf("visr.param.Dim_Reduction == T && %s",analysis_cond2)
visr.app.category("Dimensionality Reduction",active.condition = paste(dim_red_cond,dim_red_cond2,sep = "||"))
visr.param("Run_PCA", default = F, debugvalue = F, active.condition = dim_red_cond,
           info = "Run a PCA dimensionality reduction")
visr.param("nPC_compute", min = 2, default = 20, debugvalue = 20,label = "Number of PCs to compute", type = "int",
           active.condition = paste("visr.param.Run_PCA == T",dim_red_cond, sep = "&&"),
           info = "Total Number of PCs to compute and store")
visr.param("jackstraw", label = "Run Jackstraw", default = F,active.condition = dim_red_cond,
           info = "Randomly permutes a subset of data, and calculates projected PCA scores for these 'random' genes. Then compares the PCA scores for the 'random' genes with the observed PCA scores to determine statistical signifance. End result is a p-value for each gene's association with each principal component.")
visr.param("jackstrawRep",label = "Number of replicates", min = 10, default = 100, type = "int",
           active.condition = paste(dim_red_cond,"visr.param.jackstraw == T",sep = "&&"),
           info = "Number of replicate samplings to perform")
visr.param("nPC_jackstraw", label = "Number of PCs to compute signigicance for", default = 12, min = 1, type = "int",
           active.condition = paste("visr.param.jackstraw == T", dim_red_cond, sep = "&&"),
           info = "Number of PCs to include on jackstraw plot")
visr.param("elbow", label = "Plot Elbow Plot", default = F, debugvalue = T, active.condition = dim_red_cond,
           info = "Plots the standard deviation of each principle component")

visr.param("PC_heatmap","Plot PC Heatmap", default = F, active.condition = dim_red_cond,
           info = "Draws a heatmap focusing on each principal component. Both cells and genes are sorted by their principal component scores. Allows for nice visualization of sources of heterogeneity in the dataset.")
visr.param("nPC_PCheatmap", label = "Number of PCs to plot", default = 12, min = 1, type = "int",
           active.condition = paste("visr.param.PC_heatmap == T",dim_red_cond,sep = "&&"),
           info = "Number of PCs to include on PC heatmap")
visr.param("nGene_PCHeatmap", label = "Number of Genes to plot", default = 30, min = 1, type = "int",
           active.condition = paste("visr.param.PC_heatmap == T",dim_red_cond, sep = "&&"),
           info = "Number of top genes to plot for each PC")
visr.param("nCell_PCHeatmap", label = "Number of Cells to plot", default = 500, min = 1, type = "int",
           active.condition = paste("visr.param.PC_heatmap == T",dim_red_cond, sep = "&&"),
           info = "Number of top cells to plot for each PC")

visr.param("Run_tSNE", default = F, active.condition = dim_red_cond,
           info = "Run t-SNE dimensionality reduction on selected PCs")
visr.param("calculate_tsne_nPC", label = "Automatically calculate number of PCs", default = T,
           active.condition = paste("visr.param.Run_tSNE == T",dim_red_cond, sep = "&&"),
           info = "Automatically calculate number of PCs selected for tSNE. Uncheck to specify the number of PCs to use.")
visr.param("tsne_nPC", label = "Number of PCs for calculating tSNE", min = 2, default = 10, type = "int",
           active.condition = paste("visr.param.calculate_tsne_nPC == F && visr.param.Run_tSNE == T",dim_red_cond,sep = "&&"))

visr.param("Run_CCA", default = F, debugvalue = F, active.condition = dim_red_cond2,
           info = "Runs a canonical correlation analysis using a diagonal implementation of CCA.")
visr.param("nCC_compute", min = 2, default = 20, debugvalue = 20, type = "int", label = "Number of CCs to compute",
           active.condition = paste("visr.param.Run_CCA == T",dim_red_cond2, sep = "&&"),
           info = "Number of canonical vectors to calculate")
visr.param("bicor", label = "Plot CC bicor saturation plot", default = F, debugvalue = F, active.condition = dim_red_cond2,
           info = "Plots biweight midcorrelation (bicor) of the Xth gene ranked by minimum bicor across the specified CCs for each group in the grouping.var")
visr.param("align_CCA", default = F, debugvalue = F, active.condition = dim_red_cond2, label = "Align CCA subspaces",
           info = "Aligns subspaces across a given grouping variable")
visr.param("nCC_align", min = 2, default = 10, type = "int", label = "Number of CCs to align",
           active.condition = paste("visr.param.align_CCA == T",dim_red_cond2, sep = "&&"),
           info = "Number of CCs to align for downstream analysis")
visr.param("Run_tSNE2", label = "Run tSNE", default = F, active.condition = dim_red_cond2,
           info = "Run t-SNE dimensionality reduction on selected CCs")
visr.param("tsne_nCC", label = "Number of aligned CCs for calculating tSNE", min = 2, default = 10, type = "int",
           active.condition = paste("visr.param.Run_tSNE2 == T",dim_red_cond2,sep = "&&"))


run_PCA <- function(gbmData){
  if (length(gbmData@var.genes) == 0){
    gbmData <- find_variable_genes(gbmData)
  }
  
  # Scale data (regression)
  print(paste("Scaling data"))
  gbmData <- ScaleData(object = gbmData, vars.to.regress = c("nUMI", "percent.mito"))
  
  print(paste("Running PCA"))
  gbmData <- RunPCA(object = gbmData, pc.genes = gbmData@var.genes, do.print=F, 
                    pcs.compute = visr.param.nPC_compute)
  return(gbmData)
}

draw_elbow <- function(gbmData){
  if (is.null(gbmData@dr$pca)){
    gbmData <- run_PCA(gbmData)
  }
  print(paste("Drawing elbow plot"))
  p <- PCElbowPlot(object = gbmData)
  p <- p + ggtitle("Elbow Plot") + theme(plot.title = element_text(lineheight=2,size = 20,face = "plain",hjust = 0.5), plot.margin = margin(20, 10, 10, 10))
  
  print(p)
  switchPlotToScreen()
  print(p)
  switchPlotToReport()
  return(gbmData)
}

run_jackstraw <- function(gbmData){
  if (is.null(gbmData@dr$pca)){
    gbmData <- run_PCA(gbmData)
  }
  print(paste("Running JackStraw"))
  nPC <- min(visr.param.nPC_jackstraw,length(gbmData@dr$pca@sdev))
  gbmData <- JackStraw(object = gbmData, num.replicate = visr.param.jackstrawRep, num.pc = nPC)
  JackStrawPlot(object = gbmData, PCs = 1:nPC)
  switchPlotToScreen()
  JackStrawPlot(object = gbmData, PCs = 1:nPC)
  switchPlotToReport()
  return(gbmData)
}

draw_PCHeatmap <- function(gbmData){
  if (is.null(gbmData@dr$pca)){
    gbmData <- run_PCA(gbmData)
  }
  nPC <- min(visr.param.nPC_PCheatmap,length(gbmData@dr$pca@sdev))
  nCell <- min(visr.param.nCell_PCHeatmap,nrow(gbmData@dr$pca@cell.embeddings))
  nGene <- min(visr.param.nGene_PCHeatmap,nrow(gbmData@dr$pca@gene.loadings))
  
  par(oma = c(0, 0, 3, 0))
  PCHeatmap(object = gbmData, pc.use = 1:nPC, cells.use = nCell, num.genes = nGene, do.balanced = T, label.columns = F, use.full = F)
  mtext(text = sprintf("PC Heatmap for top %d cells and %d genes", nCell, nGene), outer = T, cex = 1.5,line = 1)
  
  switchPlotToScreen()
  par(oma = c(0, 0, 3, 0))
  PCHeatmap(object = gbmData, pc.use = 1:nPC, cells.use = nCell, num.genes = nGene, do.balanced = T, label.columns = F, use.full = F)
  mtext(text = sprintf("PC Heatmap for top %d cells and %d genes", nCell, nGene), outer = T, cex = 1.5,line = 1)
  switchPlotToReport()
  
  return(gbmData)
}

calculate_nPC <- function(gbmData){
  variance <- gbmData@dr$pca@sdev^2
  mean_var <- sum(variance)/length(variance)
  percent.var <- variance/sum(variance)
  min_change <- 0.001
  for (i in 2:length(percent.var)){
    #print(paste(c(i,variance[i],sum(percent.var[1:i])))
    if ((variance[i]<mean_var) & ((percent.var[i-1]-percent.var[i]) < min_change)){
      break
    }
  }
  auto_num_pc_to_use <- i
  print(paste("auto_num_pc = ",auto_num_pc_to_use))
  return(auto_num_pc_to_use)
}

run_tSNE <- function(gbmData){
  if (is.null(gbmData@dr$pca)){
    gbmData <- run_PCA(gbmData)
  }
  
  if (visr.param.calculate_tsne_nPC){
    num_pc_to_use <- calculate_nPC(gbmData)
  }else{
    num_pc_to_use <- min(visr.param.tsne_nPC,length(gbmData@dr$pca@sdev))
  }
  print(sprintf("Number of PCs: %s",num_pc_to_use))
  print(paste("Running tSNE"))
  gbmData <- RunTSNE(object = gbmData, dims.use = 1:num_pc_to_use, do.fast = T)
  p <- plot_tsne(gbmData)
  
  print(p)
  switchPlotToScreen()
  print(p)
  switchPlotToReport()
  return(gbmData)
}

plot_dim <- function(gbmData, reduction){
  title <- if (reduction == "cca"){"CCA"}else{"Aligned CCA"}
  dim <- if (reduction == "cca"){"CC1"}else{"ACC1"}
  p1 <- DimPlot(object = gbmData, reduction.use = reduction, group.by = "group", 
                pt.size = 0.5, do.return = TRUE)
  p2 <- VlnPlot(object = gbmData, features.plot = dim, group.by = "group", 
                do.return = TRUE)
  p <- plot_grid(p1, p2)
  p <- p + ggtitle(sprintf("%s Plot", title)) + theme(plot.title = element_text(lineheight=2,size = 20,face = "plain",hjust = 0.5), plot.margin = margin(20, 10, 10, 10))
  return(p)
}

run_CCA <- function(gbmData){
  groups <- levels(as.factor(gbmData@meta.data$group))
  g1 <- rownames(subset(gbmData@meta.data,group == groups[1]))
  g2 <- rownames(subset(gbmData@meta.data,group == groups[2]))
  genes.use <- gbmData@var.genes
  
  gbmData1 <- SubsetData(object = gbmData, cells.use = g1, do.clean = T)
  gbmData2 <- SubsetData(object = gbmData, cells.use = g2, do.clean = T)
  rm(gbmData)
  gbmData1 <- ScaleData(gbmData1)
  gbmData2 <- ScaleData(gbmData2)
  
  gbmData <- RunCCA(object = gbmData1, object2 = gbmData2, genes.use = genes.use, num.cc = visr.param.nCC_compute)
  p <- plot_dim(gbmData, "cca")
  print(p)
  switchPlotToScreen()
  print(p)
  switchPlotToReport()
  return(gbmData)
}

plot_bicor <- function(gbmData){
  if (is.null(gbmData@dr$cca)){
    gbmData <- run_CCA(gbmData)
  }
  print("Plotting CC bicor saturation plot")
  p <- MetageneBicorPlot(gbmData, grouping.var = "group", dims.eval = 1:visr.param.nCC_compute, display.progress = T)

  switchPlotToScreen()
  print(p)
  switchPlotToReport()
  return(gbmData)
}

align_subspace <- function(gbmData){
  if (is.null(gbmData@dr$cca)){
    gbmData <- run_CCA(gbmData)
  }
  nCC_align <- min(visr.param.nCC_align,ncol(gbmData@dr$cca@cell.embeddings))
  gbmData <- AlignSubspace(gbmData, reduction.type = "cca", grouping.var = "group", dims.align = 1:nCC_align)
  p <- plot_dim(gbmData, "cca.aligned")
  print(p)
  switchPlotToScreen()
  print(p)
  switchPlotToReport()
  return(gbmData)
}

run_tSNE2 <- function(gbmData){
  if (is.null(gbmData@dr$cca.aligned)){
    gbmData <- align_subspace(gbmData)
  }
  num_CC_to_use <- min(visr.param.tsne_nCC,ncol(gbmData@dr$cca.aligned@cell.embeddings))
  print(paste("Running tSNE"))
  gbmData <- RunTSNE(object = gbmData, reduction.use = "cca.aligned",dims.use = 1:num_CC_to_use, do.fast = T)
  p <- plot_tsne(gbmData)
  
  print(p)
  switchPlotToScreen()
  print(p)
  switchPlotToReport()
  return(gbmData)
}

plot_tsne <- function(gbmData){
  group <- if (is.null(gbmData@meta.data$group) || (length(levels(gbmData@ident)))>1){"ident"}else{"group"}
  
  p <- TSNEPlot(object = gbmData,do.return = T,group.by = group)
  p <- p + ggtitle("t-SNE Plot") + theme(plot.title = element_text(lineheight=2,size = 20,face = "plain",hjust = 0.5), plot.margin = margin(20, 10, 10, 10))
  return(p)
}