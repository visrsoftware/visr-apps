
usePackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {install.packages(pkg, repos = "http://cran.us.r-project.org");require(pkg, character.only = TRUE)}}

usePackage("psych")

visr.applyParameters()

{{
  output_fa <- fa(r = input_table[,input_columns],
     nfactors = input_nfactors,
     n.iter = input_niter, 
     rotate = input_rotate,
     scores = input_scores, 
     SMC = input_smc,
     covar = input_covar,
     residuals = input_residuals, 
     missing = input_missing,
     impute = input_impute,
     min.err = input_minerr,
     max.iter = input_maxiter,
     warnings = input_warnings, 
     fm = input_fm,
     alpha = input_alpha,
     p = input_p,
     oblique.scores = input_obliquescores)
     #use = input_use
}}

print(capture.output(print(output_fa)),collapse='\n')

# plot factor 1 by factor 2 
load <- output_fa$loadings[,1:2] 
plot(load,type="n") # set up plot 
text(load,labels=names(input_table[,input_columns]),cex=.7) # add variable names
