# write your R code below
library("ggplot2")
print(colnames(visr.input))
print(visr.input.colnames)
visr.input$gccol <- kmeans(data.frame(visr.input$H3K27me3, visr.input$H3K4me3/3), 4)$cluster
#plot(visr.input$H3K27me3, visr.input$H3K4me3, main="Regions per Chromosome", col=gccol)
gp<-qplot(data=visr.input,x=RColumnName)#, aes(fill=chr)) 
print(gp)