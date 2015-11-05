# write your R code below
paste(capture.output(print(summary(visr.input))),collapse='\n')
plot_input <- c("index")
visr.applyParameters()
plot(subset(visr.input, select = plot_input), main="Summary")