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

extractCliquesFromDag <- function(dag, root=NULL) {
  if (sum(diag(as(dag,"matrix")))!=0){
    dag <- removeSelfLoops(dag)
  }
  moral <- moralize(dag)
  tg    <- triangulate(moral)
  ripped <- rip(tg, root=root)
  if (length(ripped)==0)
    return(NULL)
  ripped$cliques
}

symbolCount <- function(x,complete=NULL) {
  uname <- unique(x)
  if (!is.null(complete)) {
    if (length(setdiff(uname,complete)) != 0)
      stop("Invalid complete specification")
    sc <- rep(0,length(complete))
    names(sc) <- as.character(complete)
  } else {
    sc <- rep(0,length(uname))
    names(sc) <- as.character(uname)
  }
  
  for (i in x) {
    sc[i] <- sc[i]+1
  }
  return(sc)
}

extractStarts <- function(edgeMatrix, genes) {
  if(!is.matrix(edgeMatrix))
    return(edgeMatrix[1])
  
  genes <- as.character(genes)
  sons <- symbolCount(edgeMatrix[,2], genes)
  fathers <- symbolCount(edgeMatrix[,1], genes)
  names(fathers[sons == 0 & fathers != 0])
}

extractEnds <- function(edgeMatrix, genes) {
  if(!is.matrix(edgeMatrix))
    return(edgeMatrix[1])
  genes <- as.character(genes)
  sons <- symbolCount(edgeMatrix[,2], genes)
  fathers <- symbolCount(edgeMatrix[,1], genes)
  names(fathers[fathers == 0 & sons != 0])
}

edgeList <- function(edgeMat){
  if (!(is.matrix(edgeMat))){
    return(edgeMat)
  }
  edgeVector <- NULL
  for (i in 1:NROW(edgeMat)){
    edgeVector<-c(edgeVector,edgeMat[i,])
  }
  return(edgeVector)
}

nameCliques <- function(cliques) {
  sapply(cliques, function(x) {
    paste(sort(x), collapse=";")
  })
}
