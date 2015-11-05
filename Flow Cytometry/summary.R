# write your R code below
#paste(capture.output(print(summary(input_table))),collapse='\n')
library(scales)
visr.applyParameters()
inputtable <- data.frame(visr.input) 
dataTable <- subset(visr.input, select = input_columns)


## define a colour pallette

visr.param.color_map<-gsub("\\s","", strsplit(visr.param.color_map,",")[[1]])
visr.param.color_map<-gsub("0x","#", visr.param.color_map)  # color_palette is like "#D73027" "#FC8D59" "#FEE090"...

if(visr.param.point_shape == "filled circle"){
	pch = 20
} else if(visr.param.point_shape == "circle"){
	pch = 1
}

final_color <- "black"

# scatter plot matrix of continuous variables, coloured by input_color
if(input_color != "") {
	if(input_color_factor) {
		palette(visr.param.color_map)
		input_color_modified <- as.factor(inputtable[,input_color])
		final_color <- input_color_modified
	} else {
		input_color_modified <- inputtable[,input_color]
		heatmap_color <- colorRampPalette(visr.param.color_map)(n = length(input_color_modified)-1)
		final_color <- heatmap_color
	}
}

plot(dataTable, main=input_title, col=alpha(final_color, input_point_alpha), pch=pch, cex=input_point_size/3)
