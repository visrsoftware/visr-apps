usePackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {install.packages(pkg, repos = "http://cran.us.r-project.org");require(pkg, character.only = TRUE)}}
usePackage("corrgram")

visr.applyParameters()

if(input_order == "none") {input_order <- FALSE}
if(input_order == "") {input_order <- NULL}
{{
  corrgram(input_table[,input_columns], 
         type = input_type,
         order=input_order, 
         lower.panel=input_lowerpanel,
         upper.panel=input_upperpanel, 
         text.panel=panel.txt,
         diag.panel=panel.minmax, 
         main= input_main,
         sub = input_sub)
}}