

visr.applyParameters()

{{
  output_friedman <- friedman.test(y = input_table[,input_variable],
                                   groups = input_table[,input_groups],
                                   blocks = input_table[,input_blocks])
}}
print(capture.output(print(output_friedman)),collapse='\n')