#setwd('~/Research/svn/Projects/GeneCalc/visr/srcR')
source("visrutils.R")

visr.library("gplots")
visr.library("RColorBrewer")
visr.library("functional")

#input_table <- read.table("~/Downloads/R42_Timecourse_Zscores7filtered\ for\ D1\ and\ D3mir-snordFiltered.txt", header=T,sep="\t",check.names=F)
#input_table <- read.table("~/Research/svn/Projects/GeneCalc/Distribution/data/Muscle.txt", header=T,sep="\t",check.names=F)

#visr.param.columns<-c("Tom-dmg-d1.FPKM","Tom-dmg-d2.FPKM","Tom-dmg-d3.FPKM","Tom-dmg-d4.FPKM","Tom-dmg-d5.FPKM","Tom-dmg-d7.FPKM","Tom-dmg-d10.FPKM","Tom-dmg-d14.FPKM")
#visr.param.columns<-c("WholeMuscle","MusclePellet","Linminus","Linplus","LinSca1minus","LinSca1plus")

### Exporting parameters to R environment...
visr.applyParameters()

isRowDendogram = (visr.param.dendrogram == "both" || visr.param.dendrogram == "row")
isColDendogram = (visr.param.dendrogram == "both" || visr.param.dendrogram == "column")

columnLabels <- ""
if (visr.param.labCol) {
  columnLabels <- visr.param.columns
}

if (visr.param.main=="")
  visr.param.main <- NULL
if (visr.param.xlab=="")
  visr.param.xlab <- NULL
if (visr.param.ylab=="")
  visr.param.ylab <- NULL


visr.param.color_map<-gsub("\\s","", strsplit(visr.param.color_map,",")[[1]])
visr.param.color_map<-gsub("0x","#", visr.param.color_map)


##################################
#  Prepare data
##################################

# workaround to allow heatmap to be drawn when there is a single column
if (length(visr.param.columns) == 1)
  visr.param.columns = c(visr.param.columns, visr.param.columns)

heatmapData <- subset(input_table, select = visr.param.columns)
heatmapMatrix <- as.matrix(heatmapData)

if (isRowDendogram && nrow(input_table) > 10000) {
  #TODO: use fastcluster for large tables. https://cran.r-project.org/web/packages/fastcluster/fastcluster.pdf
  visr.message(paste("You have requested row dendogram for a large number of rows (", nrow(input_table), "). This may take a while. Continue (ignore)?"), type="error")
}

#debug: heatmapMatrix[1,] <- NA

heatmapRows <- c(1 : nrow(heatmapMatrix))

# remove NA and NAN values
if (visr.param.narm) {
  heatmapRowMeans <- rowMeans(heatmapMatrix, na.rm = visr.param.narm)
  heatmapRows <- heatmapRows[(!is.na(heatmapRowMeans) & !is.nan(heatmapRowMeans))]
}

heatmapMatrix <- heatmapMatrix[heatmapRows,]

#### hierarchical clustering on columns
hcCol <- NULL
if (isColDendogram)
  hcCol<-hclust(dist(t(heatmapMatrix)))

#### hierarchical clustering on rows
hcRow <- NULL
sortedHeatmapRows <- heatmapRows
if (isRowDendogram) {
  hcRow<-hclust(dist(heatmapMatrix))
  sortedHeatmapRows <- heatmapRows[rev(hcRow$order)]
} else if (visr.param.sort_column != "") {
  # sort the rows by the specified column
  heatmapRows <- heatmapRows[order(input_table[heatmapRows, visr.param.sort_column], decreasing=visr.param.sort_decreasing)]
  heatmapMatrix <- heatmapMatrix[heatmapRows,]
  sortedHeatmapRows <- heatmapRows
}


# subset of rows between startIndex and endIndex
if (visr.param.startIndex < 0)
  visr.param.startIndex = length(heatmapRows) + visr.param.startIndex + 1

if (visr.param.endIndex < 0)
  visr.param.endIndex = length(heatmapRows) + visr.param.endIndex + 1

if (visr.param.startIndex < 1 || visr.param.endIndex > length(heatmapRows) || visr.param.startIndex > visr.param.endIndex)
  visr.message(paste("Invalid range Indices start index = ", visr.param.startIndex , ", end index = ", visr.param.endIndex), type="error")

# check if the range has changed
if (visr.param.startIndex > 1 || visr.param.endIndex < length(heatmapRows)) {
  # recalculate the range
  heatmapRows <- sortedHeatmapRows[visr.param.startIndex:visr.param.endIndex]
  heatmapMatrix <- as.matrix(heatmapData)
  heatmapMatrix <- heatmapMatrix[heatmapRows,]

  if (isRowDendogram) {
    hcRow<-hclust(dist(heatmapMatrix))
  }
}


##################################
#  Color Normalization
##################################

# remove old breaks
if (exists("breaks"))
  rm(breaks)

minData <- min(heatmapMatrix)
maxData <- max(heatmapMatrix)

if (visr.param.manualscaling == TRUE) {  #if clamp(), update minData/ maxData with input clamp min/max
  if (visr.param.manualmin >= visr.param.manualmax)
    visr.message("'clamp max' should be greater than 'clamp min'", type="error")

  minData <- visr.param.manualmin
  maxData <- visr.param.manualmax
}

#create breaks for color map
if (visr.param.normalize == "log10(x)" || visr.param.normalize == "log10(x+1)") { #if log ()
  if (visr.param.normalize == "log10(x+1)") {
    minData <- minData + 1
    maxData <- maxData + 1
  }

  # prevent log(0)
  if (minData <= 0) {
    visr.message("There are values <= 0. Log scale will not work.", type="error")
  }
  breaks = 10^seq(log10(minData), log10(maxData), by = (log10(maxData) - log10(minData)) / visr.param.colormap_count)
} else {
  breaks=seq(minData, maxData, by = ((maxData - minData) / visr.param.colormap_count))
}

# if not clamping, append breaks with the actual min and max;
# otherwise, clamp with the input clamp min and max
if (!visr.param.manualscaling) {
  minData <- min(heatmapMatrix)
  maxData <- max(heatmapMatrix)
  if (head(breaks,1) > minData)
    breaks=append(breaks, minData, 0)
  if (tail(breaks,1) < maxData)
    breaks = append(breaks, maxData)
}
heatmap_color <- colorRampPalette(visr.param.color_map)(n = length(breaks)-1)



#fnorm<-function(x) {
#  return (pmin((x-min(x))/(2*(median(x)-min(x))), 1.0))
#}
#
#if (visr.param.normalize == "linear [0 .. 1]") {
#  heatmapMatrix<-apply(heatmapMatrix, 2, fnorm)
#}
#
#if(visr.param.normalize == "log10(x)") {
#  for(i in 1:ncol(heatmapMatrix)) {
#    if (is.numeric(heatmapMatrix[,i]))
#      heatmapMatrix[,i] <- log10(heatmapMatrix[,i])
#  }
#}
#else if(visr.param.normalize == "log10(x+1)") {
#  for(i in 1:ncol(heatmapMatrix)) {
#    if (is.numeric(heatmapMatrix[,i]))
#      heatmapMatrix[,i] <- log10(heatmapMatrix[,i]+1)
#  }
#}
#heatmapMatrix <- data.matrix(heatmapMatrix[apply(heatmapMatrix, 1, Compose(is.finite, all)),])

#if (visr.param.manualscaling == TRUE) {
#  heatmapMatrix[heatmapMatrix <= visr.param.manualmin] = visr.param.manualmin
#  heatmapMatrix[heatmapMatrix >= visr.param.manualmax] = visr.param.manualmax
#}

#if (visr.param.label == "")
#{
#  visr.param.labRow <- ""
#   m_1 <- heatmapMatrix
#} else {
#
#}

##################################
#  Color labels
##################################

visr.param.labRow <- visr.param.label
rowLabels <- ""
if (visr.param.labRow != "") {
	# get the row labels
  #rowLabels <- subset(input_table, select = c(visr.param.columns, visr.param.label))
  rowLabels <- input_table[heatmapRows, visr.param.label]
}

if(visr.param.key == FALSE)
  visr.param.keysize = 0.5


layout(mat = rbind(4:3,2:1), widths = c(1.5,4), heights = c(1.5,4))


# Set visr.param.bordercolor with alpha zero, if no border
if (!visr.param.border) {
  visr.param.bordercolor<-"#00000000"
}

dendogramRow <- FALSE
dendogramCol <- FALSE

if (!is.null(hcRow))
  dendogramRow <- as.dendrogram(hcRow)

if (!is.null(hcCol))
  dendogramCol <- as.dendrogram(hcCol)

heatmapResult <-
  heatmap.2 (
    x = heatmapMatrix,
    # dendrogram control
    Rowv = dendogramRow, #dd,#isRowDendogram,
    Colv = dendogramCol, #isColDendogram,
    #distfun = dist,
    #hclustfun = hclust,
    dendrogram = visr.param.dendrogram,
    #symm = visr.param.symm,

    # data scaling
    scale = visr.param.scale,
    na.rm= visr.param.narm,

    # image plot
    #revC = identical(Colv, "Rowv"),
    #add.expr,

    # mapping data to colors
    #breaks = visr.param.breaks,
	  breaks=breaks,
    #symbreaks= visr.param.symbreaks,

    # colors
    col = heatmap_color,

    # block sepration
    colsep = 1:ncol(heatmapMatrix),
    rowsep = 1:nrow(heatmapMatrix),
    sepcolor= visr.param.bordercolor,
    sepwidth= c(visr.param.borderwidth,visr.param.borderwidth),

    # cell labeling
    #cellnote,
    #notecex=1.0,
    #notecol="cyan",
    #na.color=par("bg"),

    # level trace
    trace = visr.param.trace,
    tracecol = visr.param.tracecol,
    #hline=median(breaks),
    #vline=median(breaks),
    #linecol=tracecol,

    # Row/Column Labeling
    margins = c(visr.param.margin1, visr.param.margin2),
    #ColSideColors,
    #RowSideColors,
    cexRow = visr.param.rowLabelSize,
    cexCol = visr.param.colLabelSize,
    #cexRow = 0.2 + 1/log10(nr),
    #cexCol = 0.2 + 1/log10(nc),
    labRow = rowLabels,
    labCol = columnLabels,
    srtRow = visr.param.srtRow,
    srtCol = visr.param.srtCol,
    #adjRow = c(0,NA),
    #adjCol = c(NA,0),
    #offsetRow = 0.5,
    #offsetCol = 0.5,

    # color key + density info
    key = visr.param.key,
    keysize = visr.param.keysize,
    density.info = visr.param.densityinfo,
    #denscol=visr.param.tracecol,
    #symkey = min(x < 0, na.rm=TRUE) || symbreaks,
    #densadj = 0.25,

    # plot labels
    main = visr.param.main,
    xlab = visr.param.xlab,
    ylab = visr.param.ylab,
	  key.xlab="value", key.ylab="count", key.title = "",

    #key.xtickfun = function() {
    #return (list (at    = parent.frame()$scale01(c(key.breaks[1],key.breaks[length(key.breaks)])),
    #              labels= c(as.character(key.breaks[1]), as.character(tail(key.breaks,1)))))}

    # plot layout
    #lmat = NULL,
    #lhei = NULL,
    #lwid = NULL

    #lmat=rbind(c(2),c(3),c(1),c(4)),lhei=c(1,1,10,1), lwid=c(1)
  )


