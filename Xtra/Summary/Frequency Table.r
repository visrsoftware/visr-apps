

###applyParameters


if (input_main=="") input_main <- "Mosaic plot"

output_freqtable <- table(input_table[,input_columns])
#ftable(output_freqtable)
print(capture.output(print(output_freqtable)),collapse='\n')
output_proptable <- prop.table(output_freqtable)
print(capture.output(print(output_proptable)),collapse='\n')

{{
  mosaicplot(output_freqtable,
             main = input_main,
             color = input_color,
             shade = input_shade)
}}

output_allmargin<-c()
{{
  if (input_margin == TRUE){
    for (i in 1:length(input_columns)) {
      assign(eval(paste("output_margin", i, sep="")),margin.table(output_freqtable, i))
      output_allmargin <- c(output_allmargin, capture.output(eval(parse(text = paste("output_margin", i, sep="")))))
    }
  }
}}
print(capture.output(print(output_allmargin)),collapse='\n')

output_allprop<-c()
{{
  if (input_prop == TRUE){
    for (i in 1:length(input_columns)) {
      assign(eval(paste("output_prop", i, sep="")),prop.table(output_freqtable, i))
      output_allprop <- c(output_allprop, capture.output(eval(parse(text = paste("output_prop", i, sep="")))))
    }
  }
}}
print(capture.output(print(output_allprop)),collapse='\n')
