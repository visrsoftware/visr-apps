source("visrutils.R")

visr.app.start("Simple Boxplot", debugdata=mtcars)
visr.param("y", type="column-numerical", debugvalue="mpg")
visr.param("group", type="column", debugvalue="cyl")
visr.app.end(printjson=TRUE, writefile=TRUE)
visr.applyParameters()

boxplot(visr.input[[visr.param.y]]~visr.input[[visr.param.group]], xlab = visr.param.group, ylab = visr.param.y)