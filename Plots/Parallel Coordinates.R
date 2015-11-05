
usePackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {install.packages(pkg, repos = "http://cran.us.r-project.org");require(pkg, character.only = TRUE)}}

usePackage("MASS")

visr.applyParameters()


parcoord(x= input_table[,input_columns], col = input_table[,input_col], lty = input_lty, var.label = input_varlabel)