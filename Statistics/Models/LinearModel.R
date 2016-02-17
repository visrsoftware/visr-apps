

visr.applyParameters()


input_formula<-paste(c(input_y, paste(input_columns,collapse = '+')), collapse="~")

{{
  if (input_weights != ""){
    output_lm <- lm(eval(input_formula),
                   data = input_table,
                   weights = input_table[,input_weights])
  }else{
    output_lm <- lm(eval(input_formula),
                       data = input_table)
  }
}}

print(capture.output(print(summary(output_lm))),collapse='\n')
par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
plot(output_lm, ask = FALSE)

output_residuals <- output_lm$residuals     
output_fittedvalues <- output_lm$fitted.values




