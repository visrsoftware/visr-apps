

# Copyright 2012 Paolo Martini <paolo.martini@unipd.it>
#
#
# This file is part of clipper.
#
# clipper is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License
# version 3 as published by the Free Software Foundation.
#
# clipper is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public
# License along with clipper. If not, see <http://www.gnu.org/licenses/>.

runCoreClipper <- function(cliqueTest, pathList, trZero, thr, maxGap){  
  cliques <- cliqueTest$cliques
  clNames <- nameCliques(cliques)
  alphas <- cliqueTest$alpha
  names(alphas) <- clNames

  sapply(pathList, function(idx) {
    formatBestSubPath(alphas, idx, trZero, thr, maxGap)
  })
}

extractTimeDependencies <- function(clipped, ct) {
    res <- vector(length=NROW(clipped), "list")
    for (i in 1:NROW(clipped)){
        cliques <- as.numeric(unlist(strsplit(clipped[i,"involvedCliques"], ",")))
        res[[i]] <- do.call(rbind, ct$pcs[cliques])
    }
    return(res)
}

timeClipSpaced <- function(expr, times, graph, npc=1, robust=FALSE, root=NULL, trZero=0.001, signThr=0.05, maxGap=1, eqids=c()){
    
  if (NROW(expr)==0){
    warning("Expression matrix has 0 rows.")
    return(NULL)
  }
  
  expr <- t(getExpression(expr, times, TRUE))
  
  expGenes <- row.names(expr)
  genes <- nodes(graph)
  genes <- intersect(genes, expGenes)
  if (length(genes)== 0)
    stop("There is no intersection between expression feature names and the node names on the graph.")
  graph <- subGraph(genes, graph)
  
  ct <- cliqueSpacedTimeCourseTest(expr, times, graph, npc, robust, root, eqids)
  
  if (is.null(ct)){
    return(NULL)
  }
  
  jtp <- getJunctionTreePaths(graph, root)
  
  if (is.null(jtp))
    return(NULL)
  
  clipped <- runCoreClipper(ct, jtp, trZero, signThr, maxGap)
  clpprNames <- c("startIdx", "endIdx", "maxIdx", "lenght", "maxScore", "aScore", "activation", "impact", "involvedCliques", "cliqueOnPath", "involvedGenes", "pathGenes")
  
  if (!is.matrix(clipped))
    clipped <- as.matrix(clipped)
  
  clipped <- t(clipped)
  colnames(clipped) <- clpprNames
  clipped <- addNames(clipped)
  clipped <- removeNArows(clipped)
  if (NROW(clipped)==0)
    return(NULL)
  clipped <- as.data.frame(clipped, stringsAsFactors=FALSE)
  list(clipped=clipped, pathMat=extractTimeDependencies(clipped, ct))
}

timeClipEqui <- function(expr, times, graph, npc=1, robust=FALSE, root=NULL, trZero=0.001, signThr=0.05, maxGap=1){
  if (NROW(expr)==0){
    warning("Expression matrix has 0 rows.")
    return(NULL)
  }
  
  expr <- t(getExpression(expr, times, TRUE))
  expGenes <- row.names(expr)
  genes <- nodes(graph)
  genes <- intersect(genes, expGenes)
  if (length(genes)== 0)
    stop("There is no intersection between expression feature names and the node names on the graph.")
  graph <- subGraph(genes, graph)
  
  ct <- cliqueEquiTimeCourseTest(expr, times, graph, npc, robust, root)
  
  if (is.null(ct)){
    return(NULL)
  }
  
  jtp <- getJunctionTreePaths(graph, root)
  
  if (is.null(jtp))
    return(NULL)
  
  clipped <- runCoreClipper(ct, jtp, trZero, signThr, maxGap)
  clpprNames <- c("startIdx", "endIdx", "maxIdx", "lenght", "maxScore", "aScore", "activation", "impact", "involvedCliques", "cliqueOnPath", "involvedGenes", "pathGenes")
  
  if (!is.matrix(clipped))
    clipped <- as.matrix(clipped)
  
  clipped <- t(clipped)
  colnames(clipped) <- clpprNames
  clipped <- addNames(clipped)
  clipped <- removeNArows(clipped)
  if (NROW(clipped)==0)
    return(NULL)
  clipped <- as.data.frame(clipped, stringsAsFactors=FALSE)
  list(clipped=clipped, pathMat=extractTimeDependencies(clipped, ct))
}
