source("visrutils.R")

##############################################################################
# Parameters
##############################################################################

visr.app.start("Clustering",
               info = "Multiple clustering algorithms",
               debugdata = iris)

##############################################################################
# Input
visr.category("Input")
visr.param("input_columns", type = "multi-column-numerical",
           debugvalue = c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width" ))

CLUSTERING_KMEANS = "k-means"
CLUSTERING_SPECTRAL = "Spectral Clustering"
CLUSTERING_DBSCAN = "DBSCAN"

visr.category("Clustering methods")
visr.param("clustering_method", items=c(CLUSTERING_KMEANS, CLUSTERING_SPECTRAL, CLUSTERING_DBSCAN))

##############################################################################
# k-means parameters
visr.category("k-means Parameters", active.condition = paste0("visr.param.clustering_method =='", CLUSTERING_KMEANS, "'"))
visr.param("kmeans_k", default = 3L, min = 1L, max = 100L)
visr.param("kmeans_algorithm", items = c("Hartigan-Wong", "Lloyd", "MacQueen"))
visr.param.kmeans_itermax <- 50 # not exposed to UI

# spectral clustering parameters
visr.category("Spectral Clustering Parameters", active.condition = paste0("visr.param.clustering_method =='", CLUSTERING_SPECTRAL, "'"))
visr.param("spectral_k", default = 3L, min = 1L, max = 100L)
visr.param("spectral_kernel",
           items = c("rbfdot", "polydot", "vanilladot", "tanhdot", "laplacedot", "besseldot", "anovadot", "splinedot", "stringdot"),
           item.labels = c("Radial Basis (Gaussian)", "Polynomial", "Linear", "Hyperbolic", "Laplacian", "Bessel", "ANOVA RBF", "Spline", "String"))

##############################################################################
# DBScan parameters
visr.category("DBScan Parameters", active.condition = paste0("visr.param.clustering_method =='", CLUSTERING_DBSCAN, "'"))

visr.param("dbscan_epsilon", default = 1.0, min = 0.0, max = 1000.0,
           info = "Size of the epsilon neighborhood", debugvalue = 0.4)
visr.param("dbscan_min_pts", default = 5L, min = 2L, max = 1000L,
           info = "number of minimum points in the eps region (for core points)")
visr.param("dbscan_border_pts", default = TRUE)

visr.category("Other Parameters")
visr.param("sortids", label="Sort Clusters",
           info = "Sort clusters based on cluster size before assigning IDs. So largest cluster will always get the cluster ID = 1",
           default = TRUE)


##############################################################################
# Output parameters
visr.category("Clustering Output")
visr.param("cluster_ids", type = "output-column", default = "clusterID", debugvalue = "cid")
visr.param("quality_criteria", type = "output-table", options = "importRowNames=false", default = "quality_criteria", debugvalue = "qc")

##############################################################################
visr.category("Plot Options")
PLOT_TYPE_TSNE = "2D scatter plot: t-SNE"
PLOT_TYPE_PCA = "2D scatter plot: PCA"
PLOT_TYPE_MDS = "2D scatter plot: MDS"
visr.param("plot_type", items = c(PLOT_TYPE_PCA, PLOT_TYPE_MDS, PLOT_TYPE_TSNE), default=PLOT_TYPE_TSNE)
visr.param("draw_hull", default = TRUE)
visr.param("plotxy", info = "Column name prefix to output XY coordinates of the plot", type = "output-multi-column", debugvalue = "plotxy")

visr.app.end(printjson = TRUE, writefile = TRUE)
visr.applyParameters()

if (!visr.isGUI()) { # RStudio debug
  visr.input <- visr.readDataTable("~/SFU/visrseq-prototypes/Data/BobRoss/BobRoss.tsv")
  visr.param.input_columns <- c(colnames(visr.input)[6:ncol(visr.input)])
  visr.param.clustering_method = CLUSTERING_KMEANS
  visr.param.kmeans_k = 10
  #visr.param.dbscan_min_pts <- 10
  #visr.param.dbscan_epsilon <- 1.65
  #visr.param.dbscan_border_pts <- TRUE
  visr.param.plot_type <- PLOT_TYPE_TSNE
}

##############################################################################
# Compute clusters
##############################################################################

visr.library("dbscan") # used for dbscan and also hullplot and adjustcolor

x <- visr.input[, visr.param.input_columns]
x_matrix <- as.matrix(sapply(x, as.numeric))

if (visr.param.clustering_method == CLUSTERING_KMEANS) {
  ####### k-means
  kmeans_result <- kmeans(x, centers = visr.param.kmeans_k, iter.max = visr.param.kmeans_itermax, algorithm = visr.param.kmeans_algorithm)
  cluster_ids <- kmeans_result$cluster
} else if (visr.param.clustering_method == CLUSTERING_SPECTRAL) {
  ####### spectral
  visr.library("kernlab")
  specc_result <- specc(x_matrix, centers = visr.param.spectral_k, kernel = visr.param.spectral_kernel)
  cluster_ids <- specc_result@.Data
} else if (visr.param.clustering_method == CLUSTERING_DBSCAN) {
  ####### DBSCAN
  cluster_result <- dbscan(x, eps = visr.param.dbscan_epsilon, minPts = visr.param.dbscan_min_pts, borderPoints = visr.param.dbscan_border_pts)
  cluster_ids <- cluster_result$cluster
}

if (visr.param.sortids) {
  ids_sorted <- dimnames(sort(table(cluster_ids), decreasing = T))[[1]]
  sorted_clusterid <- cluster_ids
  for (c in 1:length(ids_sorted)) {
    sorted_clusterid[which(cluster_ids == ids_sorted[c])] <- c
  }
  cluster_ids <- sorted_clusterid
}


visr.param.cluster_ids <- as.factor(cluster_ids)

################################################
# Output plots
################################################
col = adjustcolor(visr.var.customPalette25, alpha.f = 0.7)

if (ncol(x) <= 2) {
  x2d = x
} else if (visr.param.plot_type == PLOT_TYPE_PCA) {
  x2d <- prcomp(x)$x[, 1:2] # first two principal components
} else if (visr.param.plot_type == PLOT_TYPE_MDS) {
  d <- dist(x) #, method = "euclidean")
  fit <- cmdscale(d, eig = TRUE) # k = 2
  x2d <- fit$points
  colnames(x2d) <- c("MDS-1", "MDS-2")
} else if (visr.param.plot_type == PLOT_TYPE_TSNE) {
  visr.library("Rtsne")
  # ecb = function(x,y){plot(x, col = col[cluster_ids%%length(col) + 1L], pch=19)} #plot(x,t='n'); text(x,labels=iris$Species, col=colors[iris$Species]) }
  set.seed(31415) # Set a fixed seed to get reproducible plot
  tsne_out <- Rtsne(x_matrix, check_duplicates = FALSE) #tsne(x_matrix, initial_dims = ncol(x_matrix), max_iter = 1000, epoch = 25, epoch_callback = ecb, perplexity=0.75)
  rm(.Random.seed, envir=globalenv()) # back to random seed
  x2d <- tsne_out$Y
  colnames(x2d) <- c("t-SNE-1", "t-SNE-2")
}

visr.library("ggplot2")
x2d_colnames <- colnames(x2d)
colnames(x2d) <- c("x", "y")
plot(x2d, pch=19, col = col[cluster_ids%%length(col) + 1L])

p <-
  ggplot(data.frame(x2d), aes(x, y)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray")) +
  geom_point(colour = col[cluster_ids%%length(col) + 1L]) +
  labs(x = x2d_colnames[1]) +
  labs(y = x2d_colnames[2])

if (visr.param.draw_hull) {
  # convex hull:
  # hullplot(x2d, cluster_ids, cex=1, pch=19, col=col)
  DRAW_SMOOTH = TRUE
  if (DRAW_SMOOTH) {
  visr.library("ggalt")
  p <- p + geom_encircle(aes(group=cluster_ids),
                         s_shape=0.5,
                         expand=0.01,
                         fill = col[cluster_ids%%length(col) + 1L],
                         colour = "white",
                         alpha=0.2) + guides(fill=FALSE)
  } else {
    visr.library("alphahull")
    fortify.ashape <- function(ashape_res) {
      xdf <- data.frame(ashape_res$edges)
      xdf <- do.call(
        rbind,
        lapply(1:nrow(xdf), function(i) {
          rbind(
            data.frame(x=xdf$x1[i], y=xdf$y1[i]),
            data.frame(x=xdf$x2[i], y=xdf$y2[i])
          )
        })
      )
      xdf <- xdf[order(-1 * atan2(
        xdf$y - mean(range(xdf$y)),
        xdf$x - mean(range(xdf$x)))), c("x", "y")]
      xdf <- rbind.data.frame(xdf[nrow(xdf),], xdf[1:(nrow(xdf)-1),])
      xdf
    }

    for (cl in 1:max(cluster_ids)) {
      cl_rows <- which(cluster_ids == cl)
      if (length(cl_rows) > 1) {
        x2d_nodup <- x2d[cl_rows, ]
        # remove duplicate points so the ahull function doesn't error out
        nodup_rows <- which(!duplicated(as.matrix(as.data.frame(x2d_nodup))))
        if (length(nodup_rows) > 2) {
          x2d_nodup <- x2d_nodup[nodup_rows, ]
          # colnames(x2d_nodup) <- c("x", "y")
          alphashape1 <- ashape(x2d_nodup,alpha=10)
          p <- p + geom_polygon(data=alphashape1, aes(x, y),
                                fill=col[cl%%length(col) + 1L], alpha=0.1)
        }
      }
    }
  }
}
print(p)
#plot(x2d, pch=19, col = col[cluster_ids%%length(col) + 1L])
colnames(x2d) <- x2d_colnames
visr.param.plotxy <- x2d

################################################
# Compute internal clustering criteria
################################################
if (visr.param.quality_criteria != "") {
  visr.library("clusterCrit")
  # visr.library("cluster")
  # cat(paste(getCriteriaNames(TRUE), collapse = '",\n"'))
  critNames = c(
    "Ball_Hall",
    "Banfeld_Raftery",
    "C_index",
    "Calinski_Harabasz",
    "Davies_Bouldin",
    # "Det_Ratio", # all NAN
    "Dunn",
    "Gamma",
    "G_plus",
    "GDI11",
    "GDI12",
    "GDI13",
    "GDI21",
    "GDI22",
    "GDI23",
    "GDI31",
    "GDI32",
    "GDI33",
    "GDI41",
    "GDI42",
    "GDI43",
    "GDI51",
    "GDI52",
    "GDI53",
    # "Ksq_DetW", # all 0
    # "Log_Det_Ratio", # all 0
    "Log_SS_Ratio",
    "McClain_Rao",
    "PBM",
    "Point_Biserial",
    "Ray_Turi",
    # "Ratkowsky_Lance", # all NAN
    # "Scott_Symons", # all -Inf
    "SD_Scat",
    "SD_Dis",
    # "S_Dbw", #NAN
    "Silhouette",
    "Tau",
    "Trace_W",
    #"Trace_WiB",
    "Wemmert_Gancarski",
    "Xie_Beni"
    )
  #critNames = "all"
  critValues <- intCriteria(x_matrix, cluster_ids, critNames)
  critValues[["NumClusters"]] <- length(levels(as.factor(cluster_ids)))
  visr.param.quality_criteria <- as.data.frame(critValues)
  colnames(visr.param.quality_criteria) <- c(critNames, "NumClusters")
}
