source("visrutils.R")

visr.app.start("Relative Log Expression (RLE)", info = "The log-ratio of the read count of each gene to the median read count across samples")
visr.param("columns", type = "multi-column-numerical", label = "Columns to calculate RLE for")
visr.param("addvalue", type = "double", label = "Increment all counts by", default = 0.0, info="The number to add to all counts (necessary if there is 0 in data)")
visr.param("scale", type = "boolean", label = "Scale by library sizes", default = TRUE, info="Calculate normalization factors and scale the raw library sizes.")
visr.param("output_Normalized", type = "output-multi-column", default="SizeNorm_", label = "Size Normalized Output Columns Prefix", info ="values of columns normalized by the size factor")
visr.param("output_RLE", type = "output-multi-column", default="RLE_", label = "RLE Output Columns Prefix", info ="RLE values of columns")
visr.app.end(printjson=TRUE, writefile=T)
visr.applyParameters()
#exit()

generateRLE <- function (object) {
  nARR = dim(object)[2]
  nGEN = dim(object)[1]
  y = apply(object, 1, median)
  y_median<<-y
  rle = matrix(nrow = nGEN, ncol = nARR)
  for (i in 1:nARR) {
    x = object[, i]
    rle[, i] = log(x/y, 10)
  }
  colnames(rle) = colnames(object)
  return(rle)
}

X <- subset(visr.input, select = visr.param.columns)
X <- X + visr.param.addvalue

if (visr.param.scale) {
  visr.biocLite("edgeR")
  f <- calcNormFactors(X) # Calculate normalization factors to scale the raw library sizes.
  X <- X*f
}

visr.param.output_Normalized<-X
visr.param.output_RLE<-generateRLE(X)
par(mfrow=c(1,1))
boxplot(x = visr.param.output_RLE, xlab = "Data Columns", ylab = "RLE Value", main = "Relative Log Expression (RLE)", outline = FALSE)
