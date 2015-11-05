# write your R code below.
library("ggplot2")
gp<-qplot(chr, H3K4me3, data=visr.input,xlab="chr",main="Histogram", geom=c("boxplot"), fill=strand)
print(gp)


