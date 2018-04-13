# use package dimRed? https://cran.r-project.org/web/packages/dimRed/index.html

source("visrutils.R")

################################################
# Parameters
################################################

visr.app.start("Dimensionality Reduction",
               info = "t-Distributed Stochastic Neighbor Embedding (t-SNE)",
               debugdata = iris)

visr.category("Input")
visr.param("input_columns", type = "multi-column-numerical",
           debugvalue = c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width" ))


METHOD_PCA = "PCA"
METHOD_MDS = "MDS"
METHOD_TSNE = "TSNE"
visr.param("method", items = c(METHOD_PCA, METHOD_MDS, METHOD_TSNE))
visr.param("tsne_perplexity", default = 30,
           active.condition = paste0("visr.param.method=='", METHOD_TSNE, "'"))
visr.param("tsne_pca_scale", default = TRUE,
           info="Should data be scaled before pca is applied?",
           active.condition = paste0("visr.param.method=='", METHOD_TSNE, "'"))



# Output parameters
visr.category("Dimensionality Reduction Output")

visr.param("output_column_prefix", label = "Column name prefix for output", type = "output-multi-column", default = "Coord")

visr.app.end(printjson = TRUE, writefile = TRUE)
visr.applyParameters()

if (!visr.isGUI()) { # RStudio debug
  visr.input <- visr.readDataTable("~/SFU/visrseq-prototypes/Data/BobRoss/BobRoss.tsv")
  visr.param.input_columns <- c(colnames(visr.input)[6:ncol(visr.input)])
  visr.param.method <- METHOD_TSNE
}

################################################
# Compute clusters
################################################

visr.library("Rtsne")

x <- visr.input[, visr.param.input_columns]

x_matrix <- as.matrix(sapply(x, as.numeric))
# remove columns where all values are the same
unique_cols <- apply(x_matrix, 2, function(x) length(unique(x)) == 1)
x_matrix <- x_matrix[, which(!unique_cols)]

if (ncol(x_matrix) <= 2) {
  visr.param.output_column_prefix <- as.data.frame(x_matrix)
} else {
  # remove rows with NA and Inf
  x_sum <- apply(x_matrix, 1, sum)
  ok_index <- which(!is.na(x_sum) & !is.infinite(x_sum))
  x_matrix <- x_matrix[ok_index,]
  if (visr.param.method == METHOD_PCA) {
    x2d <- prcomp(x_matrix)$x[, 1:2] # first two principal components
  } else if (visr.param.method == METHOD_MDS) {
    d <- dist(x_matrix) #, method = "euclidean")
    fit <- cmdscale(d, eig = TRUE) # k = 2
    x2d <- fit$points
  } else if (visr.param.method == METHOD_TSNE) {

    visr.library("Rtsne")
    set.seed(31415) # Set a fixed seed to get reproducible plot
    tsne_out <- Rtsne(x_matrix, check_duplicates = FALSE,
                      perplexity = visr.param.tsne_perplexity,
                      pca_scale = visr.param.tsne_pca_scale)
    rm(.Random.seed, envir=globalenv()) # back to random seed
    x2d <- tsne_out$Y
  }
  colnames(x2d) <- paste0(visr.param.output_column_prefix, c("1", "2"))
  plot(x2d)

  x_ret <- cbind(x_sum, x_sum)
  x_ret[ok_index, 1] <- x2d[,1]
  x_ret[ok_index, 2] <- x2d[,2]
  colnames(x_ret) <- c("1", "2")
  visr.param.output_column_prefix <- x_ret
}
