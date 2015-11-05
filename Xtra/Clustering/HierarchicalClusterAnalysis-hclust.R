
input_table <- USArrests
input_distmethod <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")[1] # the distance measure to be used
input_clustmethod <- c("ward", "single", "complete", "average", "mcquitty", "median", "centroid")[1] # the agglomeration method to be used
input_labels <- NULL # A character vector of labels for the leaves of the tree
input_hang <- 0.1 # The fraction of the plot height by which labels should hang below the rest of the plot
input_axes <- TRUE
input_main <- ""
input_sub <- NULL
input_xlab <- NULL
input_ylab <- "Height"
input_k <- 5 # an integer scalar or vector with the desired number of groups
input_h <- 0.65 # numeric scalar or vector with heights where the tree should be cut
input_column <- "Murder"

visr.applyParameters()

output_dist <- dist(x = input_table, method = input_distmethod)
output_hclust<-hclust(d = output_dist, method = input_clustmethod)

if (input_main =="") input_main = "Cluster Dendrogram"
{{
  plot(x = output_hclust, 
       labels = input_labels,
       hang = input_hang,
       axes = input_axes,
       main = input_main,
       sub = input_sub,
       xlab = input_xlab,
       ylab = input_ylab)
}}

# The cluster analysis can be "sliced" horizontally to produce unique clusters either by specifying a similarity or the number of clusters desired.
# To cut the tree at a specific similarity, specify the explicit "h" argument second with the specified similarity (or "height").

{{
  if (is.numeric(input_k)){
    output_cutree <- cutree(output_hclust, k = input_k)
  } else{
    output_cutree<-cutree(output_hclust, h = input_h)
  }
}}    

# To label the dendrogram with the group IDs
plot(output_hclust, labels = as.character(output_cutree))

# Given the clusters, you can use the cluster IDs
table(input_table[,input_column], output_cutree)

# To look at the distribution of plot elevations within clusters
boxplot(input_table[,input_column] ~ output_cutree)

