# write your R code below.
library("ggplot2")
df<-data.frames(c0, c2,c3)
gp<-qplot(data=df,
	x=c2,
	xlab="chr start",maifn="Histogram")
print(gp)

