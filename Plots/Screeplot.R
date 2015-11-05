
usePackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {install.packages(pkg, repos = "http://cran.us.r-project.org");require(pkg, character.only = TRUE)}}

usePackage("psych")

visr.applyParameters()

{{
  scree(input_table[,input_columns],
        factors = input_factors,
        pc = input_pc,
        hline = input_hline,
        main = input_main)
}}

