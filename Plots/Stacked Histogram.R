#dev.off() # reset the device

##inputPath = "/Users/hyounesy/Downloads"
##indexPath = paste(inputPath, "/Deseq_RunTime.txt", sep = "")

#tt<-read.table(indexPath, header=TRUE,sep="\t", check.names = F)
#notNA<-!apply(is.na(tt), 2, any)
#tt<-tt[,c(names(which(notNA)))]


#input_table<-iris
#visr.param.columns<-c("Petal.Width","Species")
#visr.param.sort1<-"Petal.Width"
#visr.param.sort2<-"Sepal.Length"
#visr.param.sort3<-"Petal.Length"
#visr.param.enableBinning<-TRUE
#visr.param.numBins<-10
#visr.param.leftMargin <- 15
#visr.param.labelCex <- 10


source("visrutils.R")
visr.applyParameters()


tt <- subset(input_table, select = visr.param.columns)

ylabs<-colnames(tt)#c("input_pvalue","input_foldchanges")

# Binning numeric values for better (human) pattern recognition
if (visr.param.enableBinning) {
  for (yy in 1:length(tt)) {
  	if(is.numeric(tt[[ylabs[yy]]])) {
  		tt[[ylabs[yy]]] <- tapply(tt[[ylabs[yy]]], cut(tt[[ylabs[yy]]], visr.param.numBins))
  	}
  }
}

### Sorting
# method: 9
# p: 14
# adj: 15
# sum_u: 28
# sum_d: 29
# mds1: 30
# mds2: 31
#

if (visr.param.sort1 != "")
{
  sortI <- c(visr.param.sort1) #Vector of indices used to sort the table.

  if (visr.param.sort2 != "")
    sortI <- c(sortI, visr.param.sort2)

  if (visr.param.sort3 != "")
    sortI <- c(sortI, visr.param.sort3)

  sortV <- cbind(input_table[sortI])
  o <- do.call(order, sortV)
} else {
  sortI <- c("index")
  o <- c(1:nrows(tt))
}

###


y1<-tt[[""]]
par(las=2)
par(mfrow=c(length(ylabs),1),mar=c(0, visr.param.leftMargin, 0, 1))#,las=1)

x<-c(1:nrow(tt))

for (yy in 1:length(ylabs)) {
  color <- "cornflowerblue"
	pos <- match(colnames(tt)[yy], sortI)
	if (!is.na(pos))
    color = paste("chocolate", pos, sep = "") # Plots are colored according to their position in the sorting order

	ylabel = ylabs[[yy]]#paste(ylabs[[yy]], yy, sep="\n")
  y = tt[[ylabs[yy]]][o]  #  ٫range(x$VALUE, na.rm=TRUE)
  if (!is.numeric(y))
    y <- as.integer(as.factor(y))
  plot(x, y, type="h", xaxt="n", yaxt="n", xlab=NA, ylab=NA, col=color, axes = T) #;text(0, 1.5,paste(yy, ylabs[[yy]]))
  mtext(ylabel, side=2, line=1, cex = visr.param.labelCex*0.1)
}
