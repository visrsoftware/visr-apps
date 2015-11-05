# write your R code below
library("ggplot2")

carrots <- data.frame(length = rnorm(100, 6, 2))
cukes <- data.frame(length = rnorm(50, 7, 2.5))

#Now, combine your two dataframes into one.  First make a new column in each.
carrots$veg <- 1
cukes$veg <- 2

#and combine into your new data frame vegLengths
vegLengths <- rbind(carrots, cukes)
vegLengths$vegc <- as.character(vegLengths$veg)

#now make your lovely plot
gg1<-ggplot(vegLengths, aes(length, fill = vegc)) + geom_density(alpha = 0.2)
print(gg1)