source("visrutils.R")
visr.library("ggplot2")

# start parameter definition
visr.app.start("K-means Clustering", debugdata = iris)
visr.category("clustering parameters")
visr.param("columns", type = "multi-column-numerical")
visr.param("k", label = "Number of Clusters (K)", default = 3L, min = 1L)
visr.param("algorithm", label = "Algorithm", items = c("Hartigan-Wong", "Lloyd", "MacQueen"))
visr.param("sortids", label="Sort cluster IDs by size",
           info = "Sort clusters based on cluster size before assigning IDs. So largest cluster will always get the cluster ID = 1",
           default = TRUE)
visr.category("output")
visr.param("plottype", label = "Summary Plot", items = c("barplot", "splom"), item.labels = c("Barplot of Cluster Sizes", "Scatterplot Matrix"))
visr.param("output.clusterid", label = "Column name to output cluster IDs", type = "output-column", default = "clusterid")
visr.param("output.clustermeans", label = "Table name to output cluster means", type = "output-table", default = "")
visr.app.end(writefile = TRUE, printjson=TRUE)

visr.applyParameters()

#run kmeans
cluster_data<-subset(visr.input, select = visr.param.columns)
kmeans_result <- kmeans(cluster_data, visr.param.k, algorithm = visr.param.algorithm)

#prepare results
clusterid <- kmeans_result$cluster
clustercenters <- data.frame(kmeans_result$centers, row.names = NULL)
clustercenters[,"clusterID"] <- as.factor(rownames(clustercenters))
rownames(clustercenters) <- NULL

# sort results
if (visr.param.sortids) {
  ids_sorted <- dimnames(sort(table(clusterid), decreasing = T))[[1]]
  sorted_clusterid <- clusterid
  sorted_clustercenters <- clustercenters
  for (c in 1:length(ids_sorted)) {
    sorted_clusterid[which(clusterid == ids_sorted[c])] <- c
    sorted_clustercenters[ids_sorted[c], "clusterID"] <- c
  }
  clusterid <- sorted_clusterid
  clustercenters <- sorted_clustercenters
}

visr.param.output.clusterid <- as.factor(clusterid)
visr.param.output.clustermeans <- clustercenters

if (visr.param.plottype == "barplot") {
  #barplot(table(factor(clusterid)))
  c <- ggplot(visr.input, aes(factor(visr.param.output.clusterid))) +
       theme_bw() +
       geom_bar(fill = "grey40") +
       labs(title = paste("k-means clusters for k =", visr.param.k), x = "cluster id", y = "count")
  print(c)
} else {
  # warning: scatterplot matrices for large data will take forever to render
  plot(cluster_data,
       main = paste("Scatter plot matrix for k =", visr.param.k),
       col = as.integer(visr.param.output.clusterid))
}
