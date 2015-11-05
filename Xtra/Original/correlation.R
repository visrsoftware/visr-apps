usePackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {install.packages(pkg, repos = "http://cran.us.r-project.org");require(pkg, character.only = TRUE)}}
usePackage("gplots")

visr.applyParameters()

input_colRotate<-25
input_rowRotate<-45
input_hasKey<-TRUE
input_keySize<-1

output_corr <- cor(subset(visr.input, select = input_columns), method=input_method)

{{
  if (input_plot_type == "heatmap") {
    heatmap.2(output_corr, srtCol=input_colRotate, srtRow=input_rowRotate, 
              keysize = input_keySize, labRow = input_columns, col=heat.colors(256), 
              key=input_hasKey, symkey=FALSE, density.info="none", trace="none",
              cexRow=1.5, cexCol=1.5, , margins =c(10,13)
              )
    #heatmap.2(output_corr, srtCol=25, srtRow=45, keysize = 1, labRow = corr_input, col=heat.colors(256), key=TRUE, symkey=FALSE, density.info="none", trace="none", dendrogram ="row", margins =c(10,13))
  } else if (input_plot_type == "MDS") {
    d <- dist(output_corr) # euclidean distances between the rows
    fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
    x<-fit$points[,1]
    y<-fit$points[,2]
    plot(x, y)#, xlim=c(-1,4))
    text(x, y, row.names(output_corr), cex=0.8, pos=4, col="blue")
  }
}}

