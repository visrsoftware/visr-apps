
visr.applyParameters()

{{
  if (input_y != ""){
    output_wilcox <- wilcox.test(input_table[,input_x],
                                 input_table[,input_y],
                                 alternative = input_alternative,
                                 mu = input_mu,
                                 paired = input_paired,
                                 exact = input_exact,
                                 correct = input_correct,
                                 conf.int = input_confint,
                                 conf.level = input_conflevel) 
    
  }else{
    output_wilcox <- wilcox.test(input_table[,input_x],
                                 alternative = input_alternative,
                                 mu = input_mu,
                                 paired = input_paired,
                                 exact = input_exact,
                                 correct = input_correct,
                                 conf.int = input_confint,
                                 conf.level = input_conflevel) 
  }
}}
print(capture.output(print(output_wilcox)),collapse='\n')