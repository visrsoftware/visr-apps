rm(list=ls())
source("visrutils.R")

# start parameter definition
visr.app.start("kmeans", debugdata = iris)
visr.category("clustering parameters")
visr.param("columns", type = "multi-column-numerical")
visr.param("k", default = 3)
visr.param("algorithm", items = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"))
visr.category("output")
visr.param("plot.type", items = c("scatter plot", "histogram of cluster sizes"))
visr.param("plot.title", default = "kmeans results")
visr.param("output.clusterid", type = "output-column")
visr.app.end(printjson=TRUE)

visr.applyParameters()

cluster_data<-subset(visr.input, select = visr.param.columns)
visr.param.output.clusterid <- kmeans(cluster_data, visr.param.k, algorithm = visr.param.algorithm)$cluster

# plotting options
if (visr.param.plot.type == "scatter plot") {
    plot(cluster_data, main = visr.param.plot.title, col = as.integer(visr.param.output.clusterid))
} else {
  clustersTable <- table(visr.param.output.clusterid)
  lbls<-as.character(clustersTable)
  bplt<-barplot(clustersTable , xlab="cluster ID", ylab = "cluster size", main=visr.param.plot.title)
  text(y = 0, x = bplt, labels=lbls, xpd=TRUE, adj=c(0.5, -1))
}
