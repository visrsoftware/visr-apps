# Copyright 2014 Paolo Martini <paolo.martini@unipd.it>
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

pathwayTimeCourse <- function(expr, times, graph, npc, robust=FALSE, eqids=c(2,3,4,6,8,10,12,14,16,18,19,20,21,22)){
  expr <- getExpression(expr, times, TRUE)
  
  genes <- nodes(graph)
  genes <- intersect(genes, colnames(expr))
  
  if (length(genes)== 0)
    stop("There is no intersection between expression feature names and the node names on the graph.")
  
  graph <- subGraph(genes, graph)
  expr <- expr[, genes, drop=FALSE]
  
  components <- computePCAs(expr, npc, robust) #

  alphas <- apply(components, 2 , function(pc) {
    
    gfit <- lm(pc ~ times)
    
    res <- gfit$res[eqids]
    corr<-cor(res[-1],res[-length(eqids)])
    
    if (corr < 0) {
      pc <- -pc
      corr <- abs(corr)
    }
    
    ggfit <- gls(pc ~ times+I(times^2), correlation=corCAR1(corr,form=~times))
    
    summary(ggfit)$tTable[2,4]
  })
  
  idx <- which(alphas == min(alphas))
  
  list(alpha=min(alphas), pc=components[,idx, drop=FALSE], wpc=idx)
  
}


pathwayTimeCourseEqui <- function(expr, time, graph, npc, robust=FALSE){
    expr <- getExpression(expr, time, TRUE)
  
  genes <- nodes(graph)
  genes <- intersect(genes, colnames(expr))
  
  if (length(genes)== 0)
    stop("There is no intersection between expression feature names and the node names on the graph.")
  
  graph <- subGraph(genes, graph)
  expr <- expr[, genes, drop=FALSE]
  
  components <- computePCAs(expr, npc, robust) #
  
  alphas <- apply(components, 2 , function(pc) {
    
    ggfit <- gls(pc ~ time+I(time^2), correlation=corAR1(form=~time))
    
    summary(ggfit)$tTable[2,4]
  })
  
  idx <- which(alphas == min(alphas))
  
  list(alpha=min(alphas), pc=components[,idx, drop=FALSE], wpc=idx)
  
}


pathwayTimeCourseTopologycal <- function(expr, times, graph, npc, robust=FALSE, equidistributed=TRUE, eqids=c(2,3,4,6,8,10,12,14,16,18,19,20,21,22)){  
  expr <- getExpression(expr, times, TRUE)
  
  genes <- nodes(graph)
  genes <- intersect(genes, colnames(expr))
  
  if (length(genes)== 0)
    stop("There is no intersection between expression feature names and the node names on the graph.")
  
  graph <- subGraph(genes, graph)
  expr <- expr[, genes, drop=FALSE]

  cliques <- maxClique(moralize(graph))$maxCliques
  cliques <- lapply(cliques, function(x) match(x, nodes(graph)))
  maxcliques <- max(sapply(cliques, length))

  shrink <- NCOL(expr) < maxcliques 
  shrink <- TRUE                               # Always shrink
  
  components <- computePCAsTopo(expr, npc, shrink, cliques) #Covariance are shinked and estimated with qpIPF

  
  
  alphas <- apply(components, 2 , function(pc) {
    
    gfit <- lm(pc ~ times)
    
    res <- gfit$res[eqids]
    corr<-cor(res[-1],res[-length(eqids)])
    
    if (corr < 0) {
      pc <- -pc
      corr <- abs(corr)
    }
    
    ggfit <- gls(pc ~ times+I(times^2), correlation=corCAR1(corr,form=~times))
    
    summary(ggfit)$tTable[2,4]
  })
  
  idx <- which(alphas == min(alphas))
  
  list(alpha=min(alphas), pc=components[,idx, drop=FALSE], wpc=idx)
  
}

computePCAsTopo <- function(exp, npc, shrink, cliques){
  if (NCOL(exp) < npc)
    npc <- NCOL(exp)
  
  covmat <- estimateExprCov(exp, shrink)
  scovmat <- qpIPF(covmat, cliques)
  
  pcaCov<-eigen(scovmat)
  scalee <- scale(exp)
  eigenvector <- pcaCov$vectors
  scores <- scalee%*%eigenvector
  scores[, npc, drop=FALSE]
}
