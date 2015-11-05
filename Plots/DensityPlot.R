

visr.applyParameters()


{{ 
  if(input_from != 0 | input_to != 0){
  plot(density(x = input_table[,input_column],
               bw = input_bw,
               adjust = input_adjust,
               kernel = input_kernel,
               give.Rkern = input_giverkern,
               n = input_n,
               from = input_from,
               to = input_to,
               cut = input_cut),
       main = input_main,
       sub = input_sub)
}else{
  plot(density(x = input_table[,input_column],
               bw = input_bw,
               adjust = input_adjust,
               kernel = input_kernel,
               give.Rkern = input_giverkern,
               n = input_n,
               cut = input_cut),
       main = input_main,
       sub = input_sub)
}
}}
