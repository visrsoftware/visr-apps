
usePackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {install.packages(pkg, repos = "http://cran.us.r-project.org");require(pkg, character.only = TRUE)}}

usePackage("nFactors")


visr.applyParameters()

output_ev <- eigen(cor(input_table[,input_columns])) # get eigenvalues

input_subject <- nrow(input_table[,input_columns]) # number of subjects
input_var <- ncol(input_table[,input_columns]) #  number of variables

{{
  output_pa <- parallel(subject = input_subject,
                        var = input_var,
                        rep = input_rep,
                        cent = input_cent,
                        model = input_model)
}}

output_nS <- nScree(x = output_ev$values, aparallel = output_pa$eigen$qevpea)
plotnScree(output_nS)