# help with xlim: http://stackoverflow.com/questions/13842560/get-xlim-from-a-plot-in-r

input_table <- faithful
input_column <- "eruptions"
input_breaks <- "Sturges" 
# a vector giving the breakpoints between histogram cells, or a function to compute the vector of breakpoints,or a single number giving the number of cells for the histogram,or a character string naming an algorithm to compute the number of cells ,or a function to compute the number of cells
input_freq <- NULL # logical; if TRUE, the histogram graphic is a representation of frequencies, the counts component of the result; if FALSE, probability densities, component density, are plotted (so that the histogram has a total area of one). Defaults to TRUE if and only if breaks are equidistant (and probability is not specified).
#input_probability <- !input_freq 
input_includelowest <- TRUE # logical; if TRUE, an x[i] equal to the breaks value will be included in the first (or last, for right = FALSE) bar. This will be ignored (with a warning) unless breaks is a vector.
input_right <- TRUE # logical; if TRUE, the histogram cells are right-closed (left open) intervals.
input_density <- NULL # the density of shading lines, in lines per inch. 
input_angle <- 45 # the slope of shading lines, given as an angle in degrees (counter-clockwise)
input_col <- NULL # a colour to be used to fill the bars. The default of NULL yields unfilled bars.
input_border <- NULL # the color of the border around the bars. The default is to use the standard foreground color.
input_main <- ""
input_xlab <- input_column
input_ylab <- NULL
input_xlim <- c(min(input_table[,input_column]),max(input_table[,input_column]))
input_ylim <- NULL 
input_axes <- TRUE # logical. If TRUE (default), axes are draw if the plot is drawn.
input_plot <- TRUE # logical. If TRUE (default), a histogram is plotted, otherwise a list of breaks and counts is returned. In the latter case, a warning is used if (typically graphical) arguments are specified that only apply to the plot = TRUE case.
input_labels <- FALSE # logical or character. Additionally draw labels on top of bars, if not FALSE
#input_nclass <- NULL # numeric (integer). For S(-PLUS) compatibility only, nclass is equivalent to breaks for a scalar or character argument
input_warnunused <- TRUE # logical. If plot = FALSE and warn.unused = TRUE, a warning will be issued when graphical parameters are passed to hist.default()


visr.applyParameters()

if (input_main == "") {input_main = paste("Histogram of" ,input_column)}


{{
  hist(x = input_table[,input_column],
       breaks = input_breaks,
       freq = input_freq,
       #probability = input_probability,
       include.lowest = input_includelowest,
       right = input_right,
       density = input_density,
       angle = input_angle,
       col = input_col,
       border = input_border,
       main = input_main,
       xlab = input_xlab,
       ylab = input_ylab,
       xlim = input_xlim,
       ylim = input_ylim,
       axes = input_axes,
       plot = input_plot,
       labels = input_labels,
       #nclass = input_nclass,
       warn.unused = input_warnunused)
}}

