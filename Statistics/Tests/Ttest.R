
visr.applyParameters()

{{
  if (input_y == ""){
      output_t<- t.test(input_table[,input_x], 
             alternative = input_alternative,
             mu = input_mu,
             paired = input_paired,
             var_equal = input_varequal,
             conf.level = input_conflevel) 
  }else{
    output_t <- t.test(input_table[,input_x], 
           input_table[,input_y], 
           alternative = input_alternative,
           mu = input_mu,
           paired = input_paired,
           var_equal = input_varequal,
           conf.level = input_conflevel)
  }  
}}
print(capture.output(print(output_t)),collapse='\n')
