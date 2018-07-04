DE_cond1 <- sprintf("(visr.param.Find_Marker_Genes == T && %s)",analysis_cond)
DE_cond2 <- sprintf("(visr.param.Find_Marker_Genes == T && %s)",analysis_cond2)
DE_cond <- sprintf("(%s || %s)", DE_cond1, DE_cond2)

visr.app.category("Differential Expression Anlysis",active.condition = DE_cond)
visr.param("choose_de", items = c("conserved","diff"), item.labels = c("Find conserved markers","Find DE genes across conditions"),
           label = "Choose DE analysis", active.condition = DE_cond2)
choose_cluster_cond <- sprintf("(%s || (visr.param.choose_de == 'conserved' && %s))",DE_cond1,DE_cond2)
visr.param("Choose_Clusters", label = "Choose groups", default = "all", items=c("all","selected"),
           item.labels = c("Compare each cluster to the rest of the cells","Select sepcific group of clusters"),
           active.condition = choose_cluster_cond,
           info = "Choose the target group and the group for comparison")
visr.param("group_1", label = "Group 1 cluster ids (comma separated)",type = "character", default = "0,1", active.condition = sprintf("visr.param.Choose_Clusters == 'selected' && %s",choose_cluster_cond),
           info = "Group of clusters to define markers for")
visr.param("group_2", label ="Group 2 cluster ids (comma separated)", type = "character", default = "2,3", active.condition = sprintf("visr.param.Choose_Clusters == 'selected' && %s",choose_cluster_cond),
           info = "Group of clusters for comparison. If empty, assume all clusters that are not in group 1.")

visr.param("de_cluster", "Cluster ID", type = "character", active.condition = sprintf("visr.param.choose_de == 'diff' && %s", DE_cond2),
           info = "Enter the cluster of interes")

de_param_cond <- sprintf("%s || %s",DE_cond1,DE_cond2)
visr.param("DE_test_method", items = c("wilcox","bimod","roc","t","tobit","poisson","negbinom"), default = "wilcox",info = "Denotes which test to use",
           active.condition = de_param_cond)
visr.param("min_pct",label = "Min percent of cells", default = 0.1, min = 0, max = 1, active.condition = de_param_cond,
           info = "Only test genes that are detected in more than the specified fraction of cells in either of the two populations")
visr.param("min_logfc", label = "LogFC threshold", default = 0.25, min = 0, active.condition = de_param_cond,
           info = "Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells.")

get_groups <-  function(gbmData){
  #check and return the group of clusters for DE analysis
  curr_ids <- paste(levels(gbmData@ident),collapse = ",")
  group_1 <- trimws(strsplit(visr.param.group_1,split = ',')[[1]])
  if (length(group_1) == 0){
    visr.message(sprintf("Enter clusters of interest.If ignore, DE analysis will be skipped. Current cluster ids are: %s", curr_ids),type = 'warning')
    return()
  }
  group_2 <- trimws(strsplit(visr.param.group_2,split = ',')[[1]])
  for (cluster in c(group_1,group_2)){
    if (!(cluster %in% levels(gbmData@ident))){
      visr.message(sprintf("Cluster %s doesn't exist.If ignore, DE analysis will be skipped. Current cluster ids are: %s", cluster,curr_ids),type = 'warning')
      return()
    }
  }
  if (length(group_2) == 0){group_2 <- NULL}
  return(list(group_1,group_2))
}

find_markers <- function(gbmData, target,comparison){
  group_1.markers <- FindMarkers(object = gbmData, ident.1 = target, ident.2 = comparison, 
                                 test.use = visr.param.DE_test_method,min.pct = visr.param.min_pct,
                                 logfc.threshold = visr.param.min_logfc, only.pos = T)
  group_1.markers <- group_1.markers[order(group_1.markers$p_val_adj,decreasing = F),]
  table <- data.frame(group_1.markers)
  gene <- rownames(table)
  cluster <- paste(target, collapse = ",")
  pct_diff <- table$pct.1 - table$pct.2
  table <- cbind(gene,cluster,table,pct_diff)
  table <- table[,c("gene","cluster","avg_logFC","pct.1","pct.2","pct_diff","p_val","p_val_adj")]
  
  return(table)
}

# DE anlysis of one dataset
diff_exp <- function(gbmData){
  if (length(levels(gbmData@ident)) == 1){
    visr.message("Cluster cells first. If ignore, DE analysis will be skipped.", type = 'warning')
    return()
  }
  print(paste("Running DE analysis"))
  if (visr.param.Choose_Clusters == "selected"){
    groups <- get_groups(gbmData)
    if (is.null(groups)){return()}
    table <- find_markers(gbmData, target = groups[[1]], comparison = groups[[2]])
    
  }else{
    DE.markers <- FindAllMarkers(object = gbmData, test.use = visr.param.DE_test_method,only.pos = T,
                                 min.pct = visr.param.min_pct, logfc.threshold = visr.param.min_logfc)
    # table <- data.frame(DE.markers %>% group_by(cluster) %>% top_n(top_genes_n, avg_logFC))
    table <- data.frame(DE.markers)
    pct_diff <- table$pct.1 - table$pct.2
    table <- cbind(table,pct_diff)
    table <- table[,c("gene","cluster","avg_logFC","pct.1","pct.2","pct_diff","p_val","p_val_adj")]
  }
  # output table
  write.table(x = table, file = paste(output_folder, DE_output,sep = "/"), row.names = F, quote = F, sep = "\t")
}

# integrated analysis find conserved genes
diff_exp_conserved <- function(gbmData){
  if (length(levels(gbmData@ident)) == 1){
    visr.message("Cluster cells first. If ignore, DE analysis will be skipped.", type = 'warning')
    return()
  }
  dataset_labels <- levels(as.factor(gbmData@meta.data$group))
  print(paste("Running DE analysis"))
  if (visr.param.Choose_Clusters == "selected"){
    groups <- get_groups(gbmData)
    if (is.null(groups)){return()}
    group_1.markers <- FindConservedMarkers(object = gbmData, ident.1 = groups[[1]], ident.2 = groups[[2]],
                                            test.use = visr.param.DE_test_method,min.pct = visr.param.min_pct,
                                            logfc.threshold = visr.param.min_logfc, grouping.var = "group", only.pos = T)
    
    
    group_1.markers$max_pval_adj <- pmax(group_1.markers[,paste0(dataset_labels[1],"_p_val_adj")],
                                         group_1.markers[,paste0(dataset_labels[2],"_p_val_adj")])
    
    table <- group_1.markers[order(group_1.markers$max_pval_adj,decreasing = F),]
    
    gene <- rownames(table)
    gene <- data.frame(gene)
    cluster <- paste(groups[[1]], collapse = ",")
    table <- cbind(gene,cluster,table)
    
  }else{
    table <- data.frame()
    for (id in levels(gbmData@ident)){
      markers <- FindConservedMarkers(object = gbmData, ident.1 = id, test.use = visr.param.DE_test_method,
                                      min.pct = visr.param.min_pct, logfc.threshold = visr.param.min_logfc,
                                      grouping.var = "group")
      markers$max_pval_adj <- pmax(markers[,paste0(dataset_labels[1],"_p_val_adj")],markers[,paste0(dataset_labels[2],"_p_val_adj")])
      gene <- rownames(markers)
      gene <- data.frame(gene)
      markers <- cbind(gene,cluster=id,markers)
      table <- rbind(table,markers)
    }
  }
  
  groups <- colnames(table)[c(3,8)]
  groups <- gsub(pattern = '.{6}$', replacement = '', x = groups) #remove last 6 chars "_p_val"
  for (group in groups){
    pct_diff <- data.frame(table[,paste0(group,"_pct.1")] - table[,paste0(group,"_pct.2")])
    colnames(pct_diff) <- paste0(group,"_pct_diff")
    table <- cbind(table, pct_diff)
  }
  table <- table[,c(1,2,3,4,5,6,16,7,8,9,10,11,17,12,13,14,15)]
  
  # output table
  write.table(x = table, file = paste(output_folder, DE_output,sep = "/"), row.names = F, quote = F, sep = "\t")
}

# integrated analysis find DE genes across conditions
diff_exp_across <- function(gbmData){
  if (length(levels(gbmData@ident)) == 1){
    visr.message("Cluster cells first. If ignore, DE analysis will be skipped.", type = 'warning')
    return()
  }
  print(paste("Running DE analysis"))
  gbmData@meta.data$temp_id <- as.factor(paste0(gbmData@ident, "_", gbmData@meta.data$group))
  gbmData <- SetAllIdent(gbmData, id = "temp_id")
  
  groups <- levels(as.factor(gbmData@meta.data$group))
  sub_cluster_1 <- paste0(visr.param.de_cluster, "_", groups[1])
  sub_cluster_2 <- paste0(visr.param.de_cluster, "_", groups[2])
  
  table1 <- find_markers(gbmData, sub_cluster_1,sub_cluster_2)
  table2 <- find_markers(gbmData, sub_cluster_2,sub_cluster_1)
  table <- rbind(table1,table2)
  
  # output table
  write.table(x = table, file = paste(output_folder, DE_output,sep = "/"), row.names = F, quote = F, sep = "\t")
}



