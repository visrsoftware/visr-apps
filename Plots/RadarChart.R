# Library
source("visrutils.R")

################################################
# Parameters
################################################

visr.app.start("RadarChart",
               info = "Drawing radar chart (a.k.a. spider plot)",
               debugdata = iris)

visr.category("Input")
visr.param("input_columns", type = "multi-column-numerical",
           debugvalue = c("Sepal.Length","Sepal.Width","Petal.Length","Petal.Width" ))


visr.app.end(printjson = TRUE, writefile = TRUE)
visr.applyParameters()


visr.library("fmsb")

if (!visr.isGUI()) {
  visr.input <- visr.readDataTable("~/SFU/visrseq-prototypes/Data/joy_clustering_all_300/runsInfo.txt")
  visr.param.input_columns <- c(colnames(visr.input)[18:25]) #18:59

  visr.setSelectedRows(c(10:13))
}

safeMax <- function(x) {
  okIndex = which(!is.infinite(x) & !is.na(x) & x < 3e+300)
  if (length(okIndex) == 0)
    return(1)
  resultMax <- max(x[okIndex])
  resultMin <- min(x[okIndex])
  if (resultMax == resultMin) {
    return(resultMin + 1)
  }
  return(resultMax)
}

safeMin <- function(x) {
  okIndex = which(!is.infinite(x) & !is.na(x) & x > -3e+300)
  if (length(okIndex) == 0)
    return(0)
  return(min(x[okIndex]))
}

if (length(visr.param.input_columns) > 2 && !is.null(visr.getSelectedRows())) {
  # The number of columns (variables) must be more than 2.
  data <- visr.input[, visr.param.input_columns]
  #data <- apply(data, 2, normalize)
  dataMax <- apply(data, 2, safeMax)
  dataMin <- apply(data, 2, safeMin)
  data <- data[visr.getSelectedRows(), ]

  numSelected = length(visr.getSelectedRows())

  # To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each topic to show on the plot!
  data=rbind(dataMax, dataMin , data)

  # The default radar chart proposed by the library:
  # print(data)
  # radarchart(data)


  #0: no axis label. (default)
  #1: center axis label only.
  #2:  around-the-chart label only.
  #3: means both center and around-the-chart (peripheral) labels.
  #4: *.** format of 1.
  #5: *.** format of 3. Default is 0
  visr.param.axis_type = 0
  visr.param.num_segments = 5
  visr.param.font_size = 10

  # Custom the radarChart !
  old_mar <- par("mar")
  par(mar=c(0,0,0,0))

  pcols = c(rgb(191, 55, 6, 255, max=255),
            rgb(2, 166, 191, 255, max=255),
            rgb(108, 2, 194, 255, max=255),
            rgb(189, 146, 4, 255, max=255))

  pfcols = c(rgb(191, 55, 6, 64, max=255),
             rgb(2, 166, 191, 64, max=255),
             rgb(108, 2, 194, 64, max=255),
             rgb(189, 146, 4, 64, max=255))

  radarchart(data,
             axistype = visr.param.axis_type,
             seg = visr.param.num_segments,
             pty = 16, # point shape. use 32 for no points
             #custom polygon
             pcol=pcols, # vector of color codes for plot data
             plty=1, # vector of types for plot data:
             pfcol=pfcols , #vector of color codes for filling polygons
             plwd=3 ,
              #custom the grid
             cglcol="grey",
             cglty=2,
             axislabcol="grey",
             cglwd=0.5,
             #custom labels
             vlcex=visr.param.font_size * 0.1
  )
  par(mar=old_mar)
}
