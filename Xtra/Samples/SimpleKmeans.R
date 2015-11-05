source("visrutils.R")

visr.applyParameters()

cluster_data<-subset(visr.input, select = param.columns)

output.clusterid <- kmeans(cluster_data, 
                           param.k, 
                           algorithm = param.algorithm)$cluster

#plot(cluster_data, main = param.plot.title, col = as.integer(output.clusterid))

clustersTable <- table(output.clusterid)
lbls<-as.character(clustersTable)

bplt<-barplot(clustersTable , xlab="cluster ID", ylab = "cluster size", main=param.plot.title)
text(y = 0, x = bplt, labels=lbls, xpd=TRUE, adj=c(0.5, -1))  
