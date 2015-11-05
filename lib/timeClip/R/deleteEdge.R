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

deleteEdge <- function(graph, from, to){
  edgeL <- graph@edgeL
  vertex <- nodes(graph)
  from <- as.character(from)
  to <- as.character(to)

  if (!(from %in% vertex))
    stop(paste(from, " not found in graph.", sep=""))

  pos <- match(to, vertex)
  if (is.na(pos))
    stop(paste(to, " not found in graph.", sep=""))
  
  vEdges <- edgeL[[from]]$edges
  nestedPos <- match(pos, vEdges)
  if (is.na(nestedPos))
    stop(paste(from," - ",to, ": edge not found.",sep=""))

  surrogate <- vEdges[-nestedPos]
  
  if (length(surrogate) == 0)
    surrogate <- integer(0)
  edgeL[[from]]$edges <- surrogate
  graph@edgeL <- edgeL
  return(graph)
}
