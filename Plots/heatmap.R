#setwd('~/Research/svn/Projects/GeneCalc/visr/srcR')
source("visrutils.R")

visr.library("gplots")
visr.library("RColorBrewer")
visr.library("functional")

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

if (visr.param.scale == "rows-> mean(0) sd(1)") {
  heatmapMatrix <- (heatmapMatrix - rowMeans(heatmapMatrix)) / apply(heatmapMatrix, 1, sd)
} else if (visr.param.scale == "columns-> mean(0) sd(1)") {
  heatmapMatrixT <- t(heatmapMatrix)
  heatmapMatrix <- t((heatmapMatrixT - rowMeans(heatmapMatrixT)) / apply(heatmapMatrixT, 1, sd))
}
  

if (isRowDendogram && nrow(input_table) > 10000) {
  #TODO: use fastcluster for large tables. https://cran.r-project.org/web/packages/fastcluster/fastcluster.pdf
  visr.message(paste("You have requested row dendogram for a large number of rows (", nrow(input_table), "). This may take a while. Continue (ignore)?"), type="error")
}

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


if (visr.param.showRowTicks) {
  rowIndexLabels<-rep("", nrow(input_table))
  rowIndexLabels[sortedHeatmapRows] <- as.character(seq(1, length(sortedHeatmapRows)))
  rowIndexLabels[-sortedHeatmapRows[c(visr.param.startIndex, visr.param.endIndex, seq(0, length(sortedHeatmapRows),visr.param.rowTickInterval))]] <- "    "
}

if (nchar(visr.output.orderIndex) > 0) { # an output column name specified
  visr.output.orderIndex <- rep(-1, nrow(input_table))
  visr.output.orderIndex[sortedHeatmapRows] <- seq(1, length(sortedHeatmapRows))
}

# check if the range has changed
if (visr.param.startIndex > 1 || visr.param.endIndex < length(heatmapRows)) {
  # recalculate the range
  heatmapRows <- sortedHeatmapRows[visr.param.startIndex:visr.param.endIndex]
  heatmapMatrix <- as.matrix(heatmapData)
  heatmapMatrix <- heatmapMatrix[heatmapRows,]
  sortedHeatmapRows <- heatmapRows

  if (isRowDendogram) {
    if (visr.param.reorderRows) {
      hcRow<-hclust(dist(heatmapMatrix))
    } else {
      # don't apply a dendrogram as it will reorder
      hcRow<-NULL
      isRowDendogram <- FALSE
    }
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

if (visr.param.clampvalues == TRUE) {  #if clamp(), update minData/ maxData with input clamp min/max
  if (visr.param.clampmin >= visr.param.clampmax)
    visr.message("'clamp max' should be greater than 'clamp min'", type="error")

  minData <- visr.param.clampmin
  maxData <- visr.param.clampmax
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
if (!visr.param.clampvalues) {
  minData <- min(heatmapMatrix)
  maxData <- max(heatmapMatrix)
  minBreak <- min(lapply(list(breaks), diff)[[1]]) # get the minimum diff between the breaks
  if (head(breaks,1) > minData) {
    if (head(breaks,1) - minData < minBreak) {
      breaks = breaks[-1] # remove the first break element
    }
    breaks=append(breaks, minData, 0) # append to the start
  }

  if (tail(breaks,1) < maxData) {
    if (maxData - tail(breaks,1) < minBreak) {
      # to avoid very small breaks due to the rounding errors. which will cause the error: "Error in seq.default(min.raw, max.raw, by = min(diff(breaks)/4)) : 'by' argument is much too small"
      breaks = breaks[-length(breaks)] # remove the last break element
    }
    breaks = append(breaks, maxData) # append to end
  }
}
heatmap_color <- colorRampPalette(visr.param.color_map)(n = length(breaks)-1)

##################################
#  row labels
##################################

visr.param.labRow <- visr.param.label
rowLabels <- ""
if (visr.param.labRow != "") {
  # get the row labels
  rowLabels <- input_table[heatmapRows, visr.param.label]
}

if (visr.param.showRowTicks) {
  rowTickLabels <- rowIndexLabels[heatmapRows]
  rowLabels <- paste(rowTickLabels, rowLabels, sep = "   ")
}

if(visr.param.key == FALSE)
  visr.param.keysize = 0.5


layout(mat = rbind(4:3,2:1), widths = c(1.5,4), heights = c(1.5,4))


# Set visr.param.bordercolor with alpha zero, if no border
if (!visr.param.border) {
  visr.param.bordercolor<-"#00000000"
}


##################################
#  Draw heatmap
##################################

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
    dendrogram = ifelse(isRowDendogram, ifelse(isColDendogram, "both","row"), ifelse(isColDendogram, "column", "none")),
    #symm = visr.param.symm,

    # data scaling
    scale = "none", #visr.param.scale,
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


