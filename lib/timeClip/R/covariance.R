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

estimateExprCov <- function(expr, shrink) {
  if (shrink)
    unclass(cov.shrink(expr, verbose=FALSE))
  else
    cov(expr)
}

estimateCov <- function(exp1, exp2, shrink) {
  ncl1 <- NROW(exp1)
  ncl2 <- NROW(exp2)
  
  cov1 <- estimateExprCov(exp1, shrink)  
  cov2 <- estimateExprCov(exp2, shrink)

  s <- (cov1*(ncl1-1) + cov2*(ncl2-1)) / (ncl1 + ncl2 - 2)
  list(s1=cov1, s2=cov2, s=s)  
}
