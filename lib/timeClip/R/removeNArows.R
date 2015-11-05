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

removeNArows <- function(m, col=5){
  if (is.null(m)){
    warning("m is NULL")
    return(NULL)
  }
  
  if (NROW(m) == 1) {
    if (is.na(m[col]))
      return()
    return(m)
  } else {
    toRemove <- NULL
    for (i in 1:NROW(m)) {
      if (is.na(m[i,col]))
        toRemove <- c(toRemove, i)
    }
    if (!is.null(toRemove))
      return(m[-toRemove,,drop=FALSE])
    return(m)
  }
}
