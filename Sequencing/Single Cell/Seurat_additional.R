visr.app.category("Additional Options", active.condition = sprintf("%s || %s", analysis_cond, analysis_cond2))
visr.param("rename_id", label = "Rename Clusters", default = F)
visr.param("curr_id",label = "Current cluster ids (comma separated)", default = "", active.condition = "visr.param.rename_id == T")
visr.param("new_id", label = "New cluster ids (comma separated)", active.condition = "visr.param.rename_id == T")

visr.param("export_gene_list",label = "List of genes to export/visualize", default = "")
visr.param("plot_genes_scatter", label = "Scatter plot on t-SNE projections", default = F,
           active.condition = "visr.param.export_gene_list  != ''")
visr.param("gene_plot_color", label = "Color map of gene expression", type = "multi-color", default="BuPu 7", debugvalue = "gray, red",
           active.condition = "visr.param.export_gene_list  != '' && visr.param.plot_genes_scatter == T")
visr.param("plot_genes_vln", label = "Violin plot of expression distribution in each cluster", default = F,
           active.condition = "visr.param.export_gene_list  != ''")
visr.param("plot_genes_dot", label = "Dot plot of average expression in each cluster", default = F,
           active.condition = "visr.param.export_gene_list  != ''")
visr.param("plot_genes_ridge", label = "Ridge plot of expression distribution in each clusters", default = F,
           active.condition = "visr.param.export_gene_list  != ''")
visr.param("plot_genes_hm", label = "Expression heatmap grouped by clusters", default = F,
           active.condition = "visr.param.export_gene_list  != ''")


check_gene_list <- function(gbmData,gene_probes){
  valid_gene_probes <- c()
  gene_probes <- strsplit(gsub(gene_probes, pattern = ' ', replacement = ''),split = ',')[[1]]
  for (gene_name in gene_probes){
    if (!(gene_name %in% rownames(gbmData@data))){
      visr.message(sprintf("Gene '%s' is not in this dataset",gene_name),type = "warning")
    }
    else{
      valid_gene_probes <- c(valid_gene_probes, gene_name)
    }
  }
  return(valid_gene_probes)
}

rename_cluster <- function(gbmData,reduction){
  if (length(levels(gbmData@ident)) == 1){
    visr.message("Cells have not been clustered.",type = 'warning')
    return(gbmData)
  }
  curr_id <- trimws(strsplit(visr.param.curr_id,split = ',')[[1]])
  new_id <- trimws(strsplit(visr.param.new_id,split = ',')[[1]])
  if (length(curr_id) == 0  || length(new_id) == 0){
    visr.message(sprintf("Enter cluster IDs. Current cluster ids are: %s", paste(levels(gbmData@ident),collapse = ",")),type = 'warning')
    return(gbmData)
  }
  if (length(curr_id) != length(new_id)){
    visr.message("The number of current ids and new ids must be equal.", type = 'warning')
    return(gbmData)
  }
  for (cluster in curr_id){
    if (!(cluster %in% levels(gbmData@ident))){
      visr.message(sprintf("Cluster %s doesn't exist. Current cluster ids are: %s", cluster, paste(levels(gbmData@ident),collapse = ",")),type = 'warning')
      return(gbmData)
    }
  }
    
  for (i in 1:length(curr_id)){
    if (curr_id[i] != new_id[i]){
      gbmData <- RenameIdent(object = gbmData, old.ident.name = curr_id[i], new.ident.name = new_id[i])
    }
  }
  p <- plot_clusters(gbmData,reduction)
  print(p)
  switchPlotToScreen()
  print(p)
  switchPlotToReport()
  
  return(gbmData)
}

plot_genes <- function (gbm, gene_probes, projection, limits = c(0, 10), marker_size = 0.1, title = NULL){
  # Adapted from cellrangerRkit
  gene_values <- t(as.matrix(rbind(numeric(0),gbm@data[gene_probes,])))
  gene_values[gene_values < limits[1]] <- limits[1]
  gene_values[gene_values > limits[2]] <- limits[2]
  colnames(gene_values) <- gene_probes
  projection_names <- colnames(projection)
  colnames(projection) <- c("Component.1", "Component.2")
  proj_gene <- data.frame(cbind(projection, gene_values))
  proj_gene_melt <- melt(proj_gene, id.vars = c("Component.1","Component.2"))
  p <- ggplot(proj_gene_melt, aes(Component.1, Component.2)) +
    geom_point(aes(colour = value), size = marker_size) +
    facet_wrap(~variable) + scale_colour_gradient(low = "grey", high = "red", name = "val") + labs(x = projection_names[1], y = projection_names[2])
  
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  p <- p + theme_bw() + theme(plot.title = element_text(hjust = 0.5), 
                              panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p <- p + labs(size = "test2")
  return(p)
}

plot_genes2 <- function(gbmData){
  p <- FeatureHeatmap(gbmData, features.plot = valid_gene_probes, group.by = "group", pt.size = 0.25, key.position = "right", max.exp = 3,do.return = T)
  return(p)
}

visualize_gene_scatter <- function(gbmData){
  if (is.null(gbmData@dr$tsne)){
    visr.message("No t-SNE information",type = "warning")
    return()
  }
  if (length(valid_gene_probes) > 0 && !is.null(gbmData@dr$tsne)){
    # plot expression
    if (visr.param.workflow == "single"){
      p <- plot_genes(gbm = gbmData,gene_probes = valid_gene_probes,projection = gbmData@dr$tsne@cell.embeddings)
    }else{
      p <- plot_genes2(gbmData)
    }
    p <- p + visr.util.scale_color_gradient(visr.param.gene_plot_color, label = "Expression")
    p <- p + ggtitle("Expression of Selected Genes") + theme(plot.title = element_text(lineheight=2,size = 20, hjust = 0.5), plot.margin = margin(20, 10, 10, 10),legend.title = element_text(size = 8)) 
    print(p)
    switchPlotToScreen()
    print(p)
    switchPlotToReport()
  }
}

visualize_gene_vln <- function(gbmData){
  if (length(levels(gbmData@ident)) == 1){
    visr.message("No clustering information",type = 'warning')
    return()
  }
  
  p <- VlnPlot(object = gbmData, features.plot = valid_gene_probes, x.lab.rot = TRUE, nCol = ceiling(sqrt(length(valid_gene_probes))))
  p <- p + ggtitle("Expression of Selected Genes") + theme(plot.title = element_text(lineheight=2,size = 20, hjust = 0.5), plot.margin = margin(20, 10, 10, 10),legend.title = element_text(size = 8)) 
  print(p)
  switchPlotToScreen()
  print(p)
  switchPlotToReport()
}

visualize_gene_dot <- function(gbmData){
  if (length(levels(gbmData@ident)) == 1){
    visr.message("No clustering information",type = 'warning')
    return()
  }
  
  switchPlotToScreen() # dot plot automatically print plot when called
  if (is.null(gbmData@meta.data$group)){
    p <- DotPlot(object = gbmData, genes.plot = valid_gene_probes, plot.legend = T, x.lab.rot = T,do.return = T)
  }else{
    # SplitDotPlotGG cannot run if cluster or grouping.var contains '_'
    # p <- SplitDotPlotGG(object = gbmData, genes.plot = rev(valid_gene_probes), cols.use = c("blue","red"), x.lab.rot = T, plot.legend = T, do.return = T, grouping.var = "group")
    p <- DotPlot(object = gbmData, genes.plot = valid_gene_probes, plot.legend = T, x.lab.rot = T,do.return = T)
  }
  switchPlotToReport()
  
  p <- p + ggtitle("Expression of Selected Genes") + theme(plot.title = element_text(lineheight=2,size = 20, hjust = 0.5), plot.margin = margin(20, 10, 10, 10),legend.title = element_text(size = 8)) 
  print(p)
  switchPlotToScreen()
  print(p)
  switchPlotToReport()
}

visualize_gene_ridge <- function(gbmData){
  if (length(levels(gbmData@ident)) == 1){
    visr.message("No clustering information",type = 'warning')
    return()
  }
  p <- RidgePlot(object = gbmData, features.plot = valid_gene_probes, nCol = 2)
  p <- p + ggtitle("Expression of Selected Genes") + theme(plot.title = element_text(lineheight=2,size = 20,face = "plain",hjust = 0.5), plot.margin = margin(20, 10, 10, 10))
  
  print(p)
  switchPlotToScreen()
  print(p)
  switchPlotToReport()
}

visualize_gene_hm <- function(gbmData){
  if (length(levels(gbmData@ident)) == 1){
    visr.message("No clustering information",type = 'warning')
    return()
  }
  p <- DoHeatmap(object = gbmData, genes.use = valid_gene_probes, slim.col.label = TRUE, remove.key = TRUE, do.plot = F, group.label.rot = T)
  p <- p + ggtitle("Expression of Selected Genes") + theme(plot.title = element_text(lineheight=2,size = 20,face = "plain",hjust = 0.5), plot.margin = margin(20, 10, 10, 10))
  
  print(p)
  switchPlotToScreen()
  print(p)
  switchPlotToReport()
}

visualize_gene <- function(gbmData){
  if (length(valid_gene_probes) == 0){
    visr.message("No valid genes specified",type = "warning")
    return()
  }
  
  if (visr.param.plot_genes_scatter){
    visualize_gene_scatter(gbmData)  
  }
  if (visr.param.plot_genes_vln){
    visualize_gene_vln(gbmData)
  }
  if (visr.param.plot_genes_dot){
    visualize_gene_dot(gbmData)
  }
  if(visr.param.plot_genes_ridge){
    visualize_gene_ridge(gbmData)
  }
  if(visr.param.plot_genes_hm){
    visualize_gene_hm(gbmData)
  }
}

visr.param.include_id <- T
visr.param.include_umi <- T
visr.param.include_pc <- T
visr.param.nPC_export <- 2
visr.param.include_cc <- T
visr.param.nCC_export <- 2

visr.param.include_tsne <- T
visr.param.include_gene <- T

export_cells <- function(gbmData){
  barcodes <- data.frame(gbmData@cell.names)
  colnames(barcodes) <- "Barcode"
  table <- barcodes
  if (!is.null(gbmData@meta.data$group)){
    table <- cbind(table, Dataset = gbmData@meta.data$group)
  }

  # get library id
  lib_id <- get_lib_id(gbmData@cell.names)
  table <- cbind(table, Library_id = lib_id)

  # cluster id
  if (visr.param.include_id){
    if (length(levels(gbmData@ident)) > 1){
      clusterID <- data.frame(gbmData@ident)
      colnames(clusterID) <- "Cluster_ID"
      table <- cbind(table,clusterID)
    }
  }

  if (visr.param.include_umi){
    umi_counts <- data.frame(colSums(gbmData@raw.data[rownames(gbmData@data),colnames(gbmData@data)]))
    colnames(umi_counts) <- "Total UMI"
    table <- cbind(table,umi_counts)
  }

  if (visr.param.workflow == "single"){
    if (visr.param.include_pc){
      if (!is.null(gbmData@dr$pca)){
        PCs <- data.frame(gbmData@dr$pca@cell.embeddings)
        nPC <- min(visr.param.nPC_export,length(gbmData@dr$pca@sdev))
        colnames(PCs) <- colnames(gbmData@dr$pca@cell.embeddings)
        PCs <- PCs[,1:nPC]
        table <- cbind(table,PCs)
      }
    }
  }else{
    if (visr.param.include_cc){
      if(!is.null(gbmData@dr$cca.aligned)){
        ACCs <- data.frame(gbmData@dr$cca.aligned@cell.embeddings)
        nACC <- min(visr.param.nCC_export, ncol(gbmData@dr$cca.aligned@cell.embeddings))
        colnames(ACCs) <- colnames(gbmData@dr$cca.aligned@cell.embeddings)
        ACCs <- ACCs[1:nACC]
        table <- cbind(table, ACCs)
      }
    }
  }

  if (visr.param.include_tsne){
    if (!is.null(gbmData@dr$tsne)){
      tSNEs <- gbmData@dr$tsne@cell.embeddings
      colnames(tSNEs) <- colnames(gbmData@dr$tsne@cell.embeddings)
      table <- cbind(table,tSNEs)
    }
  }

  if (visr.param.include_gene){
    if (length(valid_gene_probes) == 1){
      E <- as.matrix(gbmData@data[valid_gene_probes,])
      colnames(E) <- valid_gene_probes
      table <- cbind(table,E)
    }else if (length(valid_gene_probes) > 1){
      E <- t(as.matrix(gbmData@data[valid_gene_probes,]))
      table <- cbind(table,E)
    }
  }

  row.names(table) <- NULL
  write.table(x = table, file = paste(output_folder,cell_info_output,sep = "/"),
              quote = F, row.names = F, sep = "\t")
}
