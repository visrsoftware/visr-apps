source("visrutils.R")
visr.library("ggplot2")


# start parameter definition
visr.app.start("Violin Plot", debugdata = ToothGrowth)
visr.param("y", label = "numeric variable (y)",
           type = "column-numerical", debugvalue = "len")
visr.param("x", label = "group (x)",
           type = "column", debugvalue = "dose")
visr.param("factorx", label = "treat group (x) as fator?",
           default=TRUE, debugvalue = TRUE)
visr.param("trim", label = "trim?",
           default=FALSE, debugvalue = TRUE)
visr.param("vertical", label = "vertical stacking?",
           default=FALSE, debugvalue = FALSE)
visr.param("boxplot",
           label = "add boxplot?", default=FALSE, debugvalue = TRUE)
visr.param("bpwidth", label = "boxplot width",
           default=0.1, debugvalue = 0.2)
visr.param("mean_sdl", label = "add mean +/- sd?",
           default = FALSE, debugvalue = TRUE)
visr.param("border", label = "color border",
           default = TRUE)
visr.param("fill", label = "filled",
           default = FALSE)
visr.param("legendpos", label = "legend position",
           items=c("none", "left", "right", "top", "bottom"), default="right", debugvalue="left")
visr.param("xlabel", label = "x label", default = "")
visr.param("ylabel", label = "y label",default = "")
visr.param("title", label = "title", default = "")

#visr.param("k", default = 3L)
#visr.param("algorithm", items = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"))
#visr.category("output")
#visr.param("plot.title", default = "kmeans results")
#visr.param("output.clusterid", type = "output-column")
visr.app.end(writefile = TRUE, printjson=TRUE)

visr.applyParameters()

if (nchar(visr.param.xlabel) == 0)
  visr.param.xlabel = visr.param.x

if (nchar(visr.param.ylabel) == 0)
  visr.param.ylabel = visr.param.y

if (nchar(visr.param.title) == 0)
  visr.param.title = paste("Plot of", visr.param.ylabel, "by", visr.param.xlabel)


# Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

if (visr.param.factorx)
  visr.input[,visr.param.x] <- as.factor(visr.input[,visr.param.x])

# Basic violin plot
p <- ggplot(visr.input, aes(x = visr.input[,visr.param.x],
                            y = visr.input[,visr.param.y]
                            )) +
  geom_violin(trim = visr.param.trim)

if (visr.param.border)
  p <- p + aes(color = visr.input[,visr.param.x]) +
  labs(color=visr.param.xlabel)
if (visr.param.fill)
  p <- p + aes(fill = visr.input[,visr.param.x]) +
  labs(fill=visr.param.xlabel)

#p <- p + scale_color_manual(values=c("#999999", "#E69F00"))

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


