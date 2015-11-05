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

markMultiple <- function(g) {
  d <- edgeData(g)
  if (length(d) == 0)
    return(g)

  ns <- names(d)
  for (i in 1:length(d)) {
    tp <- d[[i]]$edgeType
    if (length(grep(";", tp, fixed=T)) > 0) {
      nodes <- unlist(strsplit(ns[[i]], "|", fixed=T))
      edgeData(g, nodes[1], nodes[2], "edgeType") <- "multiple"
    }
  }

  return(g)
}


plotTimeInCytoscape <- function(graph, timePath, path=1, color="#6699FF", main="graph", layout="jgraph-spring", root=NULL, addSymbol=FALSE){
  if (!requireNamespace("RCytoscape"))
    stop("the RCytoscape package is missing")

  if (length(timePath)!=2) {
    warning(timePath)
    return(NULL)
  }

  if (length(timePath$pathMat) < path)
    stop("Path not present.")

  if (path <= 0)
    stop("Wrong path specification.")
  
  jt <- getJunctionTree(graph, root)

  if (is.null(jt)) {
    warning("Unable to compute junctionTree.")
    return(NULL)
  }
  
  cliques <- extractCliquesFromDag(graph, root)
  cliquesStr <- sapply(cliques, function(x) paste(x, collapse=","))

  g <- markMultiple(jt)
  
  g <- initEdgeAttribute(g, "edgeType", "char", "undefined")
  g <- initEdgeAttribute(g, "weight", "numeric", 1)
  cy <- CytoscapeConnection()

  if (main %in% as.character(getWindowList(cy)))
    deleteWindow(cy, main)

  ## setEdgeAttributesDirect (w, 'weight', 'int', list('1'='2'), 1)
  
  w <- new.CytoscapeWindow(main, g)
  displayGraph(w)
  setDefaultNodeColor(w, "#FFFFFF")
  setDefaultNodeBorderColor(w, "#8888FF")
  setDefaultBackgroundColor(w, "#FFFFFF")
  
  theDir <- getwd()

  ids <- createPies(timePath, path, cliquesStr)
  
  cliquesOfThePath <- unlist(strsplit(timePath$clipped[path,"cliqueOnPath"],","))
  setNodeImageDirect(w, cliquesOfThePath, paste("file://",theDir, "/tmp/grey.png", sep=""))

  setNodeAttributesDirect(w, "genes", "string", as.character(1:length(cliques)), cliquesStr)

  if (addSymbol==TRUE) {
    if (!requireNamespace("org.Mm.eg.db"))
      stop("the org.Mm.eg.db package is missing")
    sym <- org.Mm.egSYMBOL
    
    cliquesSYM <- sapply(cliques, function(x) {
      y <- mget(x, sym, ifnotfound=NA)
      paste(y, collapse=",")
    })
    
    setNodeAttributesDirect(w, "symbols", "string", as.character(1:length(cliques)), cliquesSYM)
  }
  
  for ( i in 1:length(ids)) {
    setNodeImageDirect(w, ids[i], paste("file://",theDir, "/tmp/", ids[i], ".png", sep=""))
  }

  setNodeShapeDirect (w, cliquesOfThePath, 'ellipse')
  
  layoutNetwork(w, layout.name=layout)
  redraw(w)
}

createPies <- function(timePath, path, cliquesStr) {
  dir.create("tmp", showWarnings=F)
  myPalette <- colorRampPalette(c("green", "red"))
  
  cliquesID <- unlist(strsplit(timePath$clipped[path,"involvedCliques"],","))
  trends <- timePath$pathMat[[path]]
  
  png("tmp/grey.png",  bg="transparent")
  par(mar=c(0,0,0,0), bg="transparent")
  pie(rep(1,length(trends[1,])), col="grey", clockwise=T, labels=NA)
  dev.off()
  
  for (i in 1:length(cliquesID)){
    png(paste("tmp/", cliquesID[i], ".png", sep=""),  bg="transparent")
    palette(myPalette(8))
    par(mar=c(0,0,0,0), bg="transparent")
    myCols<-cut(trends[i,], 8, labels=myPalette(8))
    pie(rep(1,length(myCols)), col=myCols, clockwise=T, labels=c(rep(NA,length(myCols))))
    dev.off()
  }
  return(cliquesID)
}

getJunctionTree <- function(graph, root=NULL) {
  if (sum(diag(as(graph, "matrix"))) != 0) {
    graph <- removeSelfLoops(graph)
  }
  ripped <- rip(triangulate(moralize(graph)), root = root)
  if (length(ripped) == 0) {
    warning("The DAG provided can not be ripped. Please check if your input graph is a DAG.")
    return(NULL)
  }
  cliques <- ripped$cliques
  if (length(cliques) == 1) {
    warning("The DAG presents only one clique.")
    return(NULL)
  }
  parents <- ripped$parents
  edges <- cbind(parents, 1:length(parents))
  edges <- edges[edges[, 1] != 0, , drop = FALSE]
  if (nrow(edges) == 0) {
    warning("The DAG presents cliques that are not connected.")
    return(NULL)
  }
  ftM2graphNEL(edges, edgemode="undirected")
}
