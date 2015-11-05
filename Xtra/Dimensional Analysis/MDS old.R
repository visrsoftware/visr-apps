#http://www.statmethods.net/advstats/mds.html
#http://sape.inf.usi.ch/quick-reference/ggplot2/colour

input_columns <- "start"

visr.applyParameters()

d <- dist(input_table[, input_columns]) # euclidean distances between the rows
fit <- cmdscale(d, eig=TRUE, k=2) # k is the number of dim
output_x<-fit$points[,1]
output_y<-fit$points[,2]
plot(output_x, output_y, xlab="Coordinate 1", ylab="Coordinate 2")
text(output_x, output_y, colnames(input_table[, input_columns]), cex=0.8, pos=4, col="blue")

