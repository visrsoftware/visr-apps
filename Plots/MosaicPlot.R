source("visrutils.R")

# http://www.statmethods.net/advgraphs/mosaic.html

usePackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {install.packages(pkg, repos = "http://cran.us.r-project.org");require(pkg, character.only = TRUE)}}

usePackage("vcd")

visr.applyParameters()


if (input_main=="") input_main <- "Mosaic plot"

output_freqtable <- table(input_table[,input_columns])

print(capture.output(print(output_freqtable)),collapse='\n')

mosaic(output_freqtable,
       convars = input_convars,
       main = input_main,
       sub = input_sub,
       shade = input_shade,
       legend = input_legend)

output_freq <- as.data.frame(output_freqtable)