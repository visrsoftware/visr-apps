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

addNames <- function(m, sep=";"){
  matrixNames <- apply(m[,1:2, drop=FALSE], 1, function(x) paste(x, collapse=sep))
  if (length(matrixNames) != length(unique(matrixNames)))
    stop("Duplicated Names")
  row.names(m) <- matrixNames
  return(m)
}

computeSimilarity <- function(x,y){
  n1 <- length(setdiff(x,y))
  n2 <- length(setdiff(y,x))
  if (n1 < n2)
    return(n1/length(x))
  return(n2/length(y))
}

getBest <- function(m, col, firsts=1){
  if (NROW(m) == 1)
    return(m)
  sorter <- as.numeric(m[,col])
  names(sorter) <- row.names(m)
  m[names(sort(sorter, decreasing=T))[firsts],, drop=FALSE]
}

prunePaths <- function(pathSummary, thr=NULL, clust=NULL, sep=";"){

  if (is.null(pathSummary)){
    warning("pathSummary is NULL")
    return(NULL)
  }

  if(NROW(pathSummary) == 0)
    return(NULL)
  
  if(NROW(pathSummary) == 1)
    return(pathSummary)

  pathNames <- row.names(pathSummary)

  similarityMatrix            <- matrix(NA, length(pathNames), length(pathNames))
  row.names(similarityMatrix) <- pathNames
  colnames(similarityMatrix)  <- pathNames
  
  ways <- pathSummary[,11]
  names(ways) <- pathNames
  
  for (i in 1:length(pathNames)-1) {
    for (j in (i+1):length(pathNames)) {
      x <- unlist(strsplit(ways[i],sep))
      y <- unlist(strsplit(ways[j],sep))
      v <- computeSimilarity(x,y)
      similarityMatrix[i,j] <- v
      similarityMatrix[j,i] <- v
    }
  }
  diag(similarityMatrix) <- 0
  hData <- hclust(as.dist(similarityMatrix))
  
  if (!is.null(clust)){
    pdf(clust)
    plot(hData, main="", xlab="paths", ylab="threshold")
    dev.off()
  }
  
  if (is.null(thr)) {
    return(pathSummary)
  }
  
  hc <- cutree(hData, h = thr)
  
  clusters <- tapply(1:length(hc), hc, function(x) names(hc[x]))
  slimMatrix <- NULL
  for (grp in clusters) {
    slimMatrix <- rbind(slimMatrix, getBest(pathSummary[grp,, drop=FALSE], 5))
    pathSummary <- pathSummary[!(row.names(pathSummary) %in% grp),, drop=FALSE]
  }
  rbind(slimMatrix, pathSummary)
}
