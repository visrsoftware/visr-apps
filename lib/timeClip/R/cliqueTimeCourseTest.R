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

computePCAs <- function(exp, npc, robust){
  if (NCOL(exp) < npc)
    npc <- NCOL(exp)
  
  if(robust){
    PcaHubert(exp,k=npc)@scores
  } else {
    pc1 <- prcomp(exp)$x[,1:npc]
    as.matrix(pc1,ncol=npc)
  }
}

runSpacedTimeCourseTest <- function(expr, time, graph,  npc, robust, root, eqids) {
  cliques <- extractCliquesFromDag(graph, root=root)
  
  pcs <- vector(length=length(cliques), "list")
  wpcs <- NULL
  iterator <- 1
  p <- NULL
  
  for (cli in cliques) {

    cliExp <- expr[, cli, drop=FALSE]
    
    components <- computePCAs(cliExp, npc, robust)

    alphas <- apply(components, 2 , function(pc) {
      gfit <- lm(pc ~ time)
      
      res <- gfit$res[eqids]
      corr<-cor(res[-1],res[-length(eqids)])
      
      if (corr<0) {
        pc <- -pc
        corr <- abs(corr)
      }
      
      ggfit <- gls(pc ~ time+I(time^2), correlation=corCAR1(corr,form=~time))
      
      summary(ggfit)$tTable[2,4]
    })
    
    idx <- which(alphas == min(alphas))
    pcs[[iterator]] <- as.numeric(components[, idx, drop=FALSE])
    wpcs <- c(wpcs, idx)
    iterator = iterator + 1
    
    p <- c(p, min(alphas))
  }
  list(alpha=p, cliques=cliques, pcs=pcs, wpcs=wpcs)
}

cliqueSpacedTimeCourseTest <- function(expr, time, graph, npc=1, robust=FALSE, root=NULL, eqids=c(2,3,4,6,8,10,12,14,16,18,19,20,21,22)) {
  expr <- getExpression(expr, classes=time, TRUE)

  genes <- nodes(graph)
  genes <- intersect(genes, colnames(expr))
  
  if (length(genes)== 0)
    stop("There is no intersection between expression feature names and the node names on the graph.")
  
  graph <- subGraph(genes, graph)
  expr <- expr[, genes, drop=FALSE]

  tryCatch(runSpacedTimeCourseTest(expr, time, graph, npc, robust, root, eqids),
           error=function(e) list(alpha=NA,
             cliques=NA,
             pcs=e$call,
             wpcs=e$message))
  
}

runEquiTimeCourseTest <- function(expr, time, graph,  npc, robust, root) {
  cliques <- extractCliquesFromDag(graph, root=root)
  
  pcs <- vector(length=length(cliques), "list")
  wpcs <- NULL
  iterator <- 1
  p <- NULL
  
  for (cli in cliques) {
    
    cliExp <- expr[, cli, drop=FALSE]
    
    components <- computePCAs(cliExp, npc, robust)
    
    alphas <- apply(components, 2 , function(pc) {
      
      ggfit <- gls(pc ~ time+I(time^2), correlation=corAR1(form=~time))
      
      summary(ggfit)$tTable[2,4]
    })
        
    idx <- which(alphas == min(alphas))
    pcs[[iterator]] <- as.numeric(components[, idx, drop=FALSE])
    wpcs <- c(wpcs, idx)
    iterator = iterator + 1
    
    p <- c(p, min(alphas))
  }
  list(alpha=p, cliques=cliques, pcs=pcs, wpcs=wpcs)
}

cliqueEquiTimeCourseTest <- function(expr, time, graph, npc=1, robust=FALSE, root=NULL) {
  expr <- getExpression(expr, classes=time, TRUE)

  genes <- nodes(graph)
  genes <- intersect(genes, colnames(expr))
  
  if (length(genes)== 0)
    stop("There is no intersection between expression feature names and the node names on the graph.")
  
  graph <- subGraph(genes, graph)
  expr <- expr[, genes, drop=FALSE]
  
  runEquiTimeCourseTest(expr, time, graph, npc, robust, root)
}
