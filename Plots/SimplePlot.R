visr.applyParameters()
{{
  param_color="black"
  if (input_color!="") param_color = as.factor(input_table[,input_color])
  plot(input_table[,input_x], 
       input_table[,input_y],
       col=param_color, 
       main=input_title, 
       xlab=input_x, 
       ylab=input_y,
       log=input_log)
}}