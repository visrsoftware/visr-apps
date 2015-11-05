


visr.applyParameters()

output_kruskal <- kruskal.test(x = input_table[,input_x], g = input_table[,input_g])
print(capture.output(print(output_kruskal)),collapse='\n')