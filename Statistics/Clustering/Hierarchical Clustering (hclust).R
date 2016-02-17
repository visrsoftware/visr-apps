visr.applyParameters()

output_dist <- dist(x = input_table[,input_columns], method = input_distmethod)
output_hclust<-hclust(d = output_dist, method = input_clustmethod)


# The cluster analysis can be "sliced" horizontally to produce unique clusters either by specifying a similarity or the number of clusters desired.
# To cut the tree at a specific similarity, specify the explicit "h" argument second with the specified similarity (or "height").

{{
  if (is.numeric(input_k)){
    output_cutree <- cutree(output_hclust, k = input_k)
  } else{
    output_cutree<-cutree(output_hclust, h = input_h)
  }
}}    


# Given the clusters, you can use the cluster IDs
output_table <- table(input_table[,input_column], output_cutree)
print(capture.output(print(output_table)),collapse='\n')


{{
  if (plotstage == "Cluster dendrogram"){
    if (input_labels == TRUE){
      plot(x = output_hclust,
           labels = NULL,
           hang = input_hang,
           axes = input_axes,
           main = input_main,
           sub = input_sub,
           xlab = input_xlab,
           ylab = input_ylab)
    }else{
      plot(x = output_hclust,
           labels = FALSE,
           hang = input_hang,
           axes = input_axes,
           main = input_main,
           sub = input_sub,
           xlab = input_xlab,
           ylab = input_ylab)
    }
  }else{
    if (plotstage == "Cluster dendrogram labeled with cluster ID"){
      plot(output_hclust, labels = as.character(output_cutree)) # To label the dendrogram with the group IDs
    }else{
      if (plotstage == "boxplot"){
        boxplot(input_table[,input_column] ~ output_cutree, main = "boxplot") # To look at the distribution of plot elevations within clusters
      }
    }
  }
}}