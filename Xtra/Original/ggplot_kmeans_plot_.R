# write your R code below
library("ggplot2")
source("multiplot.R")
print(colnames(visr.input))
print(visr.input.colnames)


cluster_input <- c("index")
cluster_k <- 4
cluster_algorithm<-"Hartigan-Wong"
scatter_x <- "H3K27me3"
scatter_y <- "H3K4me3"
histogram_x <- "RNA_seq"
histogram_xlim <- c(0,5)
text_size <- 20


visr.applyParameters()

cluster_data<-subset(visr.input, select = cluster_input)
visr.input$cluster <- as.character(kmeans(cluster_data, cluster_k, algorithm=cluster_algorithm)$cluster)
gg1<-ggplot(visr.input,aes_string(x=scatter_x,y=scatter_y, colour="cluster")) + geom_point(alpha = 0.5) + theme(text = element_text(size=text_size))
gg2<-ggplot(visr.input, aes_string(x=histogram_x, fill="cluster")) + geom_density(position="fill") + xlim(histogram_xlim) + theme(text = element_text(size=text_size))
multiplot(gg1, gg2, cols=2)
