
###applyParameters

{{
  if (input_choice == "Correlation"){
    output_co <- cor(x = input_table[,input_x],
                     use = input_use,
                     method = input_method)
  }
}}

{{
  if (input_choice == "Covariance"){
    output_co <- cov(x = input_table[,input_x],
                     use = input_use,
                     method = input_method)
  }
}}

{{
  if (input_choice == "Variance"){
    output_co <- var(x = input_table[,input_x],
        na.rm = input_narm, 
        use = input_use)
  }
}}

{{
  if(input_choice == "cov2cor"){
    output_co <- cov2cor(x = input_table[,input_x])
  }
}}

print(capture.output(print(output_co)),collapse='\n')