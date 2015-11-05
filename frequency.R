freq = function(data) {
	n = ncol(data)
	m = nrow(data)
	o = data.frame(row.names = 1:n)
 	l = levels(as.factor(as.matrix(data)))
 	lc = data.frame(row.names = 1:length(l))
 	
 	for (i in 1:(n)) {
 		lc$count = 0
 		lc$sym = l
 		for (j in 1:length(l)) {
 			lc$count[j] = length(which(as.factor(data[,i]) == lc$sym[j] ))
 		}
 		index = which.max(lc$count)
 		
 		o$sym[i] = lc$sym[index]
 		o$count[i] = lc$count[index]
 		o$freq[i] = lc$count[index] / m
 	}
 	
 	return(o)
}
