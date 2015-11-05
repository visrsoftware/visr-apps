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

setGeneric("getExpressionImpl", function(expr) standardGeneric("getExpressionImpl"))

setMethod("getExpressionImpl",
          signature("matrix"),
          function(expr) return(expr))

setMethod("getExpressionImpl",
          signature("ExpressionSet"),
          function(expr) exprs(expr))


getExpression <- function(expr, classes, timecourse=FALSE) {
    expr <- getExpressionImpl(expr)
    
    if (is.null(rownames(expr)))
        stop("Gene names not specified.")
    if (!timecourse) {        
        if (NCOL(expr) != length(classes))
            stop("Class vector length and sample number differs.")
        
        if (!all((classes == 1 | classes == 2)))
            stop("Class vector should be made by either '1' or '2'.")
        
        if (sum(classes==1) < 3)
            stop("Too few sample on class 1.")
        
        if (sum(classes==2) < 3)
            stop("Too few sample on class 1.")
    } else {
        if (any(duplicated(classes)))
            stop("Please time points must be unique.")
    }

    t(expr)
   
}
