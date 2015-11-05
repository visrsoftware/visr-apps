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

removeSelfLoops <- function(graph){
  edgeL <- graph@edgeL
  for (i in 1:length(edgeL)) {
    pos <- match(i,edgeL[[i]]$edges)
    if (!(is.na(pos)))
      edgeL[[i]]$edges <- edgeL[[i]]$edges[-pos]
  }  
  graph@edgeL <- edgeL
  return(graph)
}
