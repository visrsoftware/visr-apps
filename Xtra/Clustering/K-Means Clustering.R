
usePackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {install.packages(pkg, repos = "http://cran.us.r-project.org");require(pkg, character.only = TRUE)}}

usePackage("NbClust")


visr.applyParameters()

#1 standardize data

data <- input_table[,input_columns]
if (input_scale){ data = scale(data)}

#2 determine number of clusters

#3 K-means cluster analysis

{{
  wssplot <- function(data, nc, seed=1234){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)
  }
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
}
}}

oldpar <- par()  # make a copy of current graphical settings
set.seed(1234) # determining the best number of clusters
nc2 <- NbClust(data, min.nc = 2, max.nc = input_nc, method="kmeans")
table(nc2$Best.n[1,]) # restore original settings
par(oldpar)



{{
  output_kmeans <- kmeans(data, 
                        centers = input_centers,
                        iter.max = input_itermax,
                        nstart = input_nstart,
                        algorithm = input_algorithm)
}}

output_clusterid <- output_kmeans$cluster
print(capture.output(print(output_kmeans)),collapse='\n')
print(paste("total sum of squares", eval(output_kmeans$totss),"        between sum of squares", eval(output_kmeans$betweenss),"        total within sum of squares", eval(output_kmeans$tot.withinss)))


{{
  if (input_choiceofplot == "Plot the within groups sums of squares vs. the number of clusters extracted"){
    wssplot(data, nc = input_nc) # A bend in the graph can suggest the appropriate number of clusters.
  }else{
    if (input_choiceofplot == "Recommended number of clusters using 26 criteria provided by the NbClust package"){
      barplot(table(nc2$Best.n[1,]), 
              xlab="Numer of Clusters", ylab="Number of Criteria",
              main="Number of Clusters Chosen by 26 Criteria")
    }else{
      if (input_choiceofplot == "K-means clustering"){
        plot(data, main="K-means clustering", col=as.integer(output_clusterid))
      }
    }
  }
}}