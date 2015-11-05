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

getJunctionTreePaths <- function(graph, root=NULL) {
  if (sum(diag(as(graph,"matrix")))!=0){
    graph <- removeSelfLoops(graph)
  }
  ripped <- rip(triangulate.default(moralize.default(graph)), root=root)
  if (length(ripped)==0){
    warning("The DAG provided can not be ripped. Please check if your input graph is a DAG.")
    return(NULL)
  }
  cliques <- ripped$cliques

  if (length(cliques) == 1){
    warning("The DAG presents only one clique.")
    return(list(1))
  }
  
  parents <- ripped$parents
  edges <- cbind(parents,1:length(parents))
  edges <- edges[edges[,1] != 0,, drop=FALSE]
  
  if (nrow(edges) == 0){
    warning("The DAG presents cliques that are not connected.")
    return(NULL)
  }
  
  startingCliques <- extractStarts(edges, 1:length(parents))
  endingCliques   <- extractEnds(edges, 1:length(parents))
  junctionTree    <- graph(edgeList(edges), directed=FALSE)

  paths <- NULL
  for (s in startingCliques) {
    paths <- c(paths, get.all.shortest.paths(junctionTree, s, endingCliques)$res)
  }
  return(paths)
}
