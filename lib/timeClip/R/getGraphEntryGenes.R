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

getGraphEntryGenes <- function(graph, byCliques=FALSE, root=NULL){
  graphi <- igraph.from.graphNEL(graph)
  edgeM  <- get.edgelist(graphi)
  genes  <- nodes(graph)
  starts <- extractStarts(edgeM, genes)

  cliques <- extractCliquesFromDag(graph, root=root)

  if (is.null(cliques)){
    warning("No cliques available or the DAG provided can not be ripped. Please check if your input graph is a DAG.")
    return(NULL)
  }
  
  startOnClique <- lapply(cliques, function(x) intersect(starts,x))

  if (byCliques)
    return(startOnClique)
  return(starts)
}

