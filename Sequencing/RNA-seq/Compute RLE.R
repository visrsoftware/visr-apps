#rm(list=ls())
#input_table <- read.table("/Users/hyounesy/Research/Data/VisRseq/error/041616-Ref-seq.xls",sep="\t",header=TRUE)
#input_columns<-c("WT1count","WT2count","WT3count")
#input_columns<-c("WT1RPKM","WT2RPKM","WT3RPKM")

source("visrutils.R")

visr.applyParameters()

{{
  generateRLE <- function (object) {
    nARR = dim(object)[2]
    nGEN = dim(object)[1]
    y = apply(object, 1, median)
    y_median<<-y
    rle = matrix(nrow = nGEN, ncol = nARR)
    for (i in 1:nARR) {
      x = object[, i]
      rle[, i] = log(x/y)
    }
    colnames(rle) = colnames(object)
    return(rle)
  }
}}

output_RLE<-generateRLE(subset(input_table, select = input_columns))
par(mfrow=c(1,1))
boxplot(x = output_RLE, xlab = "Output", ylab = "RLE", main = "RLE output", outline = FALSE)
