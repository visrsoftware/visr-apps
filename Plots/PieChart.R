
visr.applyParameters()

mytable <- table(input_table[,input_column])

lbls <- names(mytable)
{{
if (input_show_percentage) {
  pct <- round(mytable/sum(mytable)*100)
  lbls <- paste(lbls, "", pct) # add percents to labels 
  lbls <- paste(lbls, "%",sep="") # ad % to labels 
} else {
  lbls <- paste(lbls, " ", mytable, sep="")
}
}}

par(cex=input_label_size/20)
par(cex.main=input_title_size/input_label_size)

if (input_title == "") input_title<-paste("Pie Chart of ", input_column)


{{
if (input_rainbow_color) {
  pie(mytable, labels = lbls, main=input_title, col=rainbow(length(lbls)))
  legend("left", legend = names(mytable), fill=rainbow(length(lbls)), cex=1)
} else {
  pie(mytable, labels = lbls, main=input_title)
}
}}
