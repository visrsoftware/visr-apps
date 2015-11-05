# https://www.stat.washington.edu/research/reports/2012/tr597.pdf

usePackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {install.packages(pkg, repos = "http://cran.us.r-project.org");require(pkg, character.only = TRUE)}}

usePackage("mclust")


input_data <- faithful
input_G <- NULL # An integer vector specifying the numbers of mixture components (clusters) for which the BIC is to be calculated.
input_modelnames <- NULL # A vector of character strings indicating the models to be fitted in the EM phase of clustering.
input_prior <- NULL
input_x <- NULL # An object of class "mclustBIC". If supplied, mclustBIC will use the settings in x to produce another object of class "mclustBIC", but with G and modelNames as specified in the arguments.

input_g <- NULL
input_modelNames <- NULL

input_G2 <- NULL
input_modelnames2 <- NULL
input_symbols <- NULL
input_colors <- NULL
input_xlab <- NULL
input_ylab <- "BIC"
input_ylim <- NULL
input_legendArgs <- list(x="bottomright", ncol=2, cex=1)

visr.applyParameters()

{{
  output_BIC <- mclustBIC(data = input_data, 
                          G = input_G, 
                          modelNames = input_modelnames, 
                          prior = input_prior, 
                          x = input_x)
}}

{{
output_summary <- summary(object = output_BIC, 
                          data = input_data,
                          G = input_g,
                          modelNames = input_modelNames)
}}
output_summary 




{{
plot(x = output_BIC, 
     G = input_G2, 
     modelNames = input_modelnames,
     symbols = input_symbols,
     colors = input_colors,
     xlab = input_xlab,
     ylab = input_ylab,
     ylim = input_ylim, 
     legendArgs = list(x = "bottomright", ncol = 5))
}}
