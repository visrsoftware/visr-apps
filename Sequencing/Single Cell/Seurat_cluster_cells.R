cluster_cond <- sprintf("((%s || %s) && visr.param.Cluster_Cells == T)",analysis_cond,analysis_cond2)
cluster_cond1 <- sprintf("(%s && visr.param.Cluster_Cells == T)",analysis_cond)
cluster_cond2 <- sprintf("(%s && visr.param.Cluster_Cells == T)",analysis_cond2)

visr.app.category("Cluster Cells", active.condition = cluster_cond)
visr.param("calculate_cluster_nPC", label = "Automatically calculate number of PCs", default = T,
           active.condition = cluster_cond1,
           info = "Automatically calculate number of PCs selected for clustring cells. Uncheck to specify the number of PCs to use.")
visr.param("cluster_nPC", label = "Number of PCs for clustering", min = 1, default = 10, type = "int",
           active.condition = sprintf("visr.param.calculate_cluster_nPC == F && %s", cluster_cond1))

visr.param("cluster_nCC", label = "Number of aligned CCs for clustering", min = 1, default = 10, type = "int",
           active.condition = cluster_cond2)
visr.param("cluster_resolution", label = "resolution", min = 0.1, default = 0.6, debugvalue = 0.6,
           info = "Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.") #increase for large dataset

plot_clusters <- function(gbmData,reduction){
  if (!is.null(gbmData@dr$tsne)){
    p <- plot_tsne(gbmData)
  } else {
    p <- DimPlot(object = gbmData, reduction.use = reduction, do.return = T)
    projection <- if (reduction == "pca") {"PC"} else {"Aligned CC"}
    p <- p + ggtitle(sprintf("Clusters on %s Projections", projection)) + theme(plot.title = element_text(lineheight=2,size = 20,face = "plain",hjust = 0.5), plot.margin = margin(20, 10, 10, 10))
  }
  return(p)
}

cluster_cells <- function(gbmData, reduction){
  if (reduction == "pca"){
    if (is.null(gbmData@dr$pca)){
      gbmData <- run_PCA(gbmData)
    }
    if (visr.param.calculate_cluster_nPC){
      num_dim_to_use <- calculate_nPC(gbmData)
    }else{
      num_dim_to_use <- min(visr.param.cluster_nPC,length(gbmData@dr$pca@sdev))
    }
    print(sprintf("Number of PCs: %s",num_dim_to_use))
  }else{
    if (is.null(gbmData@dr$cca.aligned)){
      gbmData <- align_subspace(gbmData)
    }
    num_dim_to_use <- min(visr.param.cluster_nCC,ncol(gbmData@dr$cca.aligned@cell.embeddings))
  }

  print(paste("Clustering cells"))
  gbmData <- FindClusters(object = gbmData, reduction.type = reduction, dims.use = 1:num_dim_to_use, 
                          resolution = visr.param.cluster_resolution, print.output = 0, save.SNN = T, 
                          temp.file.location = output_folder)
  # PrintFindClustersParams(object = gbmData)
  if (length(levels(gbmData@ident)) == 1){
    visr.message("Cells cannot be clustered. Try adjusting the parameters.",type = 'warning')
  }
  
  p <- plot_clusters(gbmData,reduction)
  
  print(p)
  switchPlotToScreen()
  print(p)
  switchPlotToReport()
  return(gbmData)
}

