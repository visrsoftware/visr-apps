source("visrutils.R")

# start parameter definition
visr.app.start("kmeans", debugdata = iris)
visr.category("clustering parameters")
visr.param("columns", type = "multi-column-numerical")
visr.param("k", default = 3L)
visr.param("algorithm", items = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"))
visr.category("output")
visr.param("plot.title", default = "kmeans results")
visr.param("output.clusterid", type = "output-column")
visr.app.end(writefile = TRUE, printjson=TRUE)

visr.applyParameters()

cluster_data<-subset(visr.input, select = visr.param.columns)
visr.param.output.clusterid <- kmeans(cluster_data, visr.param.k, algorithm = visr.param.algorithm)$cluster

plot(cluster_data, main = visr.param.plot.title, col = as.integer(visr.param.output.clusterid))
