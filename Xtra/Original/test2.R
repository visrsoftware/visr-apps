# write your R code below
library("ggplot2")
source("multiplot.R")
print(colnames(visr.input))
print(visr.input.colnames)

cluster_x <- data.frame(visr.input$H3K27me3, visr.input$H3K4me3/4) 
cluster_k <- 4
scatter_x <- 'H3K27me3'
scatter_y <- 'H3K4me3'
histogram_x <- 'RNA_seq'
histogram_xlim <- xlim(0,5)
text_size <- 20

# visr.Params 
#{
# "cluster_x": {"label": "input", "type":"Columns", "min":1},
# "cluster_k": {"label": "k", "type":"int", "min":1, "default":3},
# "scatter_x": {"label": "scatter plot x", "type":"ColName"},
# "scatter_y": {"label": "scatter plot y", "type":"ColName"},
# "histogram_x": {"label": "histogram x", "type":"ColName"},
# "histogram_xlim": {"label": "range", "type":"DoubleRange"}
# "text_size" : {"label": "text size", "type":"int", "default":20}
#}


visr.input$cluster <- as.character(kmeans(cluster_x, cluster_k)$cluster)
gg1<-ggplot(visr.input,aes_string(x=scatter_x,y=scatter_y, colour="cluster")) + geom_point(alpha = 0.5) + theme(text = element_text(size=text_size))
gg2<-ggplot(visr.input, aes_string(x=histogram_x, fill="cluster")) + geom_density(position="fill") + histogram_xlim + theme(text = element_text(size=text_size))
multiplot(gg1, gg2, cols=2)
