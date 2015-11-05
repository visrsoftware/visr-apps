source("getExpression.R", chdir=TRUE)
source("removeNArows.R", chdir=TRUE)
source("prunePaths.R", chdir=TRUE)
source("clipperBackbone.R", chdir=TRUE)
source("timeClip.R", chdir=TRUE)

myTimeClipSpaced <- function (expr, times, graph, npc = 1, robust = FALSE, root = NULL, 
          trZero = 0.001, signThr = 0.05, maxGap = 1, eqids = c()) 
{
  if (NROW(expr) == 0) {
    warning("Expression matrix has 0 rows.")
    return(NULL)
  }
  expr <- t(getExpression(expr, times, TRUE))
  expGenes <- row.names(expr)
  genes <- nodes(graph)
  genes <- intersect(genes, expGenes)
  if (length(genes) == 0) 
    stop("There is no intersection between expression feature names and the node names on the graph.")
  graph <- subGraph(genes, graph)
  ct <- cliqueSpacedTimeCourseTest(expr, times, graph, npc, 
                                   robust, root, eqids)
  if (is.null(ct)) {
    return(NULL)
  }
  jtp <- getJunctionTreePaths(graph, root)
  if (is.null(jtp)) 
    return(NULL)
  clipped <- runCoreClipper(ct, jtp, trZero, signThr, maxGap)

  clpprNames <- c("startIdx", "endIdx", "maxIdx", "lenght", 
                  "maxScore", "aScore", "activation", "impact", "involvedCliques", 
                  "cliqueOnPath", "involvedGenes", "pathGenes")
  if (!is.matrix(clipped)) 
    clipped <- as.matrix(clipped)
  clipped <- t(clipped)
  colnames(clipped) <- clpprNames
  clipped <- addNames(clipped)
  clipped <- removeNArows(clipped)
  if (NROW(clipped) == 0) 
    return(NULL)
  clipped <- as.data.frame(clipped, stringsAsFactors = FALSE)
  list(clipped = clipped, pathMat = extractTimeDependencies(clipped, 
                                                            ct),ct=ct)
}