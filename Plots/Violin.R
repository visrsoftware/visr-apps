source("visrutils.R")
visr.library("ggplot2")
visr.library("reshape2")

# start parameter definition
visr.app.start("Violin Plot", debugdata = ToothGrowth)
# visr.param("y", label = "numeric variable (y)",
#           type = "column-numerical", debugvalue = "len")

visr.param("y", label = "Y value column(s)",
           type = "multi-column-numerical", debugvalue = c('len','dose'))

visr.param("x", label = "X grouping column",
           type = "column", debugvalue = "supp")


visr.param("factorx", label = "Treat X group as fators",
           default=TRUE, debugvalue = TRUE)

visr.param("add_dot_plot", label = "Add dots", default=FALSE, debugvalue = TRUE)

visr.param("add_jitter", label = "Add jitter points",
           info = "Adds a small amount of random variation to the location of each point, useful for handling overplotting.",
           default=FALSE, debugvalue = TRUE,
           active.condition = "visr.param.add_dot_plot == true")

visr.param("dot_binwidth", label = "Dot bins size", default=0, min = 0,
           info = "Size of the bins. Set to 0 to use 1/30th of data range as the bin size",
           active.condition = "visr.param.add_jitter == false")

visr.param("jitter_amount", label = "Jitter amount", default=0.05, min = 0,
           info = "Amount of jitter: added in both positive and negative directions, so the total spread is twice the value specified here.",
           active.condition = "visr.param.add_jitter == true")

visr.param("dot_size", label = "Dot size", default=0.4, min = 0,
           active.condition = "visr.param.add_dot_plot == true")


visr.param("trim", label = "Trim the tails to data range", default=FALSE, debugvalue = TRUE)

visr.param("boxplot", label = "Add boxplot", default=FALSE, debugvalue = TRUE)

visr.param("bpwidth", label = "Boxplot width", default=0.1, min = 0, debugvalue = 0.2,
           active.condition = "visr.param.boxplot == true")

visr.param("mean_sdl", label = "Add mean +/- sd", default = FALSE, debugvalue = TRUE)

visr.param("vertical", label = "Vertical stacking", default=FALSE, debugvalue = FALSE)

visr.param("border", label = "Colored borders", default = TRUE)

visr.param("fill", label = "Collored fill", default = FALSE)

visr.param("legendpos", label = "Legend position",
           items=c("none", "left", "right", "top", "bottom"), default="right", debugvalue="left")

visr.param("xlabel", label = "X label", default = "")
visr.param("ylabel", label = "Y label",default = "")
visr.param("title", label = "Title", default = "")

#visr.param("k", default = 3L)
#visr.param("algorithm", items = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"))
#visr.category("output")
#visr.param("plot.title", default = "kmeans results")
#visr.param("output.clusterid", type = "output-column")
visr.app.end(writefile = TRUE, printjson=TRUE)

visr.applyParameters()

# Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

if (visr.param.x != "" && visr.param.factorx)
  visr.input[,visr.param.x] <- as.factor(visr.input[,visr.param.x])

if (length(visr.param.y) > 1) {
  # reshape data first:

  if (visr.param.x == "") {
    input_melt <- melt(visr.input[, visr.param.y])
  } else {
    input_melt <- melt(visr.input[, c(visr.param.y, visr.param.x)], id = visr.param.x)
  }
  # Basic violin plot
  p <- ggplot(input_melt, aes(x = variable, y = value))
  if (visr.param.x != "") {
    p <- p + facet_wrap(as.formula(paste("~", visr.param.x)))
  }
  visr.param.x <- "variable"
  visr.param.y <- "value"
} else {
  # Basic violin plot
  if (visr.param.x != "") {
    p <- ggplot(visr.input, aes_string(x = visr.param.x, y = visr.param.y))
  } else {
    p <- ggplot(visr.input) + aes(x = 'x') + aes_string(y = visr.param.y)
  }
}

p <- p + geom_violin(trim = visr.param.trim)

if (visr.param.add_dot_plot) {
  if (visr.param.add_jitter) {
    p <- p + geom_jitter(position=position_jitter(visr.param.jitter_amount), size = visr.param.dot_size) #shape=16,
  }
  else {
    if (visr.param.dot_binwidth == 0)
      p <- p + geom_dotplot(binaxis='y', stackdir='center', dotsize=visr.param.dot_size)
    else
      p <- p + geom_dotplot(binaxis='y',binwidth = visr.param.dot_binwidth, stackdir='center', dotsize=visr.param.dot_size)
  }
}

if (visr.param.border && visr.param.x != "")
  p <- p + aes_string(color = visr.param.x) +
  labs(color=visr.param.xlabel)

if (visr.param.fill && visr.param.x != "")
  p <- p + aes_string(fill = visr.param.x) +
  labs(fill=visr.param.xlabel)

#p <- p + scale_color_manual(values=c("#999999", "#E69F00"))

if (nchar(visr.param.xlabel) == 0)
  visr.param.xlabel = visr.param.x

if (nchar(visr.param.ylabel) == 0)
  visr.param.ylabel = visr.param.y

if (nchar(visr.param.title) == 0)
  visr.param.title = paste("Plot of", visr.param.ylabel, "by", visr.param.xlabel)


p <- p + labs(title = visr.param.title,
              x = visr.param.xlabel,
              y = visr.param.ylabel)
              #color = visr.param.xlabel)

if (visr.param.vertical)# Rotate the violin plot
  p <- p + coord_flip()

if (visr.param.boxplot)# add boxplot
  p <- p + geom_boxplot(width=visr.param.bpwidth)

if (visr.param.mean_sdl)
  p <- p + stat_summary(fun.data=data_summary,
                   geom="pointrange")


p <- p + theme_bw() + theme(
        axis.line = element_line(colour = "black"),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank()
        )
p <- p + theme(legend.position=visr.param.legendpos)


# violin plot with dot plot
# p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)
# violin plot with jittered points: 0.2 degree of jitter in x direction
# p + geom_jitter(shape=16, position=position_jitter(0.2))

print(p)


