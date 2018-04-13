source("visrutils.R")

################################################
# Parameters
################################################

visr.app.start("DBScan",
               info = paste0("DBSCAN (Density-based spatial clustering of applications with noise) clustering algorithm. ",
                             "DBSCAN estimates the density around each data point by counting the number of points in a user-specified eps-neighborhood ",
                             "and applies a used-specified minPts thresholds to identify core, border and noise points. ",
                             "In a second step, core points are joined into a cluster if they are density-reachable ",
                             "(i.e., there is a chain of core points where one falls inside the eps-neighborhood of the next). ",
                             "Finally, border points are assigned to clusters."),
               debugdata = iris)
visr.category("Input")

visr.param("input_columns", type = "multi-column-numerical",
           debugvalue = c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width" ))

visr.category("Algorithm Parameters")

visr.param("epsilon_distance", default = 1.0, min = 0.0, max = 1000.0,
           info = "Size of the epsilon neighborhood", debugvalue = 0.4)
visr.param("min_points", default = 5, min = 2, max = 1000,
           info = "number of minimum points in the eps region (for core points)")
visr.param("assign_border_points", default = TRUE)

visr.category("Output")

PLOT_TYPE_TSNE = "2D scatter plot: t-SNE"
PLOT_TYPE_PCA = "2D scatter plot: PCA"
PLOT_TYPE_MDS = "2D scatter plot: MDS"
visr.param("plot_type", items = c(PLOT_TYPE_PCA, PLOT_TYPE_MDS, PLOT_TYPE_TSNE))
visr.param("plot_covex_hull", default = FALSE)
visr.param("cluster_ids", type = "output-column", debugvalue = "cid")
visr.param("quality_criteria", type = "output-table", options = "importRowNames=false", debugvalue = "qc")

visr.app.end(printjson = TRUE, writefile = TRUE)

visr.applyParameters()

################################################
# Compute clusters
################################################

visr.library("dbscan") # used for dbscan and also hullplot and adjustcolor

x <- visr.input[, visr.param.input_columns]
cluster_result <- dbscan(x, eps = visr.param.epsilon_distance, minPts = visr.param.min_points, borderPoints = visr.param.assign_border_points)
cluster_ids <- cluster_result$cluster
visr.param.cluster_ids <- cluster_ids

################################################
# Output plots
################################################

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
  visr.library("tsne")
  x2d <- tsne(x,  initial_dims = ncol(x))
  colnames(x2d) <- c("t-SNE-1", "t-SNE-2")
}


col = adjustcolor(visr.var.customPalette25, alpha.f = 0.5)
if (visr.param.plot_covex_hull) {
  hullplot(x2d, cluster_result, cex=1, pch=19, col=col)
} else {
  plot(x2d, pch=19, col = col[cluster_ids%%length(col) + 1L])
}


################################################
# Compute internal clustering criteria
################################################
# critNames = c("ball_hall","banfeld_raftery","c_index","calinski_harabasz","davies_bouldin","det_ratio","dunn","gamma","g_plus","gdi11","gdi12","gdi13","gdi21","gdi22","gdi23","gdi31","gdi32","gdi33","gdi41","gdi42","gdi43","gdi51","gdi52","gdi53","ksq_detw","log_det_ratio","log_ss_ratio","mcclain_rao","pbm","point_biserial","ray_turi","ratkowsky_lance","scott_symons","sd_scat","sd_dis","s_dbw","silhouette","tau","trace_w","trace_wib","wemmert_gancarski","xie_beni")
if (visr.param.quality_criteria != "") {
  visr.library("clusterCrit")
  critNames = "all"
  critValues <- intCriteria(as.matrix(sapply(x, as.numeric)), cluster_ids, critNames)
  critValues[["NumClusters"]] <- length(levels(as.factor(cluster_ids)))
  visr.param.quality_criteria <- as.data.frame(critValues)
}
