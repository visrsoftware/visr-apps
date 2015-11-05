# write your R code below

#print(colnames(visr.input))
#print(visr.input.colnames)

cluster_input <- c("index")
cluster_k <- 4
cluster_algorithm<-"Hartigan-Wong"

visr.applyParameters()

cluster_data<-subset(visr.input, select = cluster_input)
cluster_id <- kmeans(cluster_data, cluster_k, algorithm=cluster_algorithm)$cluster
plot(subset(visr.input, select = cluster_input), main="Kmeans", col=as.integer(cluster_id))
