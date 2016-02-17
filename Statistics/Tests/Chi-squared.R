

visr.applyParameters()

output_table = table(input_table[,input_x], input_table[,input_y]) 
print(capture.output(print(output_table)),collapse='\n') # the contingency table 

{{
  output_chisq <- chisq.test(output_table, 
                             correct = input_correct, 
                             simulate.p.value = input_simulatepvalue,
                             B = input_B) 
}}
print(capture.output(print(output_chisq)),collapse='\n') 