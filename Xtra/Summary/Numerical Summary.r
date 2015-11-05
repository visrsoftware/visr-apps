usePackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {install.packages(pkg, repos = "http://cran.us.r-project.org");require(pkg, character.only = TRUE)}}

usePackage("RcmdrMisc")

###applyParameters
                    
input_statistics <- NULL
if (input_mean == TRUE) {input_statistics = c(input_statistics,"mean")}
if (input_sd == TRUE) {input_statistics = c(input_statistics,"sd")}
if (input_se == TRUE) {input_statistics = c(input_statistics,"se(mean)")}
if(input_IQR == TRUE) { input_statistics = c(input_statistics, "IQR")}
if(input_quantiles == TRUE) { input_statistics = c(input_statistics, "quantiles")}
if(input_cv == TRUE){input_statistics = c(input_statistics, "cv")}
if(input_skewness == TRUE) {input_statistics = c(input_statistics, "skewness")}
if(input_kurtosis == TRUE){input_statistics = c(input_statistics, "kurtosis")}

if (input_quantiles == TRUE) input_quantiles2 = eval(parse(text = paste("c(",input_quantiles2,")"))) 
{{
  if (input_groups == ""){
    if (input_skewness == TRUE | input_kurtosis == TRUE) {
      output_numsummary <- numSummary(data = input_table[,input_columns], 
                                      statistics = input_statistics, 
                                      quantiles = input_quantiles2, 
                                      type = input_type)
    }else{
      output_numsummary <- numSummary(data = input_table[,input_columns], 
                                      statistics = input_statistics, 
                                      quantiles = input_quantiles2)
    }
  }else{
    if (input_skewness == TRUE | input_kurtosis == TRUE) {
      output_numsummary <- numSummary(data = input_table[,input_columns], 
                                      statistics = input_statistics, 
                                      quantiles = input_quantiles2, 
                                      type = input_type,
                                      groups = input_table[,input_groups])
    }else{
      output_numsummary <- numSummary(data = input_table[,input_columns], 
                                      statistics = input_statistics, 
                                      quantiles = input_quantiles2, 
                                      groups = input_table[,input_groups])
    }
  }
  
}}


print(capture.output(print(output_numsummary)),collapse='\n')