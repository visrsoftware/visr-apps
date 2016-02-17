#input_table <- iris
#input_show_percentage<-FALSE
#input_column<-"Species"
#input_title_size<-20
#input_label_size<-15
#input_horizontal<-TRUE

visr.applyParameters()


#counts<-prop.table(table((output_cluster)))*length(output_cluster)


if (input_title == "") input_title <- paste("Bar Plot of ", input_column)

mytable <- table(input_table[,input_column])
lbls<-as.character(mytable)


if (input_show_percentage) {
  mytable<-mytable/sum(mytable)*100
  lbls <- as.character(round(mytable/sum(mytable)*100))
  lbls <- paste(lbls, "%",sep="") # ad % to labels 
}


par(cex=input_label_size/20)
par(cex.main=input_title_size/input_label_size)

bplt<-barplot(mytable, main=input_title, horiz=input_horizontal)


if (input_horizontal) {
  text(x= mytable, y= bplt, labels=lbls, xpd=TRUE, adj=c(0,0.5))
} else {
  text(y= mytable, x= bplt, labels=lbls, xpd=TRUE, adj=c(0.5,0))  
}
