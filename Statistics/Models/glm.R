
usePackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {install.packages(pkg, repos = "http://cran.us.r-project.org");require(pkg, character.only = TRUE)}}

usePackage("boot")

visr.applyParameters()

input_formula<-paste(c(input_y, paste(input_columns,collapse = '+')), collapse="~")

{{
  if (input_family == "quasi"){
    if (input_weights != ""){
      output_glm <- glm(formula = eval(input_formula),
                        data = input_table,
                        weights = input_table[,input_weights],
                        family = eval(parse(text = paste(input_family,"(",input_link,",",input_variance,")"))))
    }else{
      output_glm <- glm(formula = eval(input_formula),
                        data = input_table,
                        family = eval(parse(text = paste(input_family,"(",input_link,",",input_variance,")"))))
      
    }
  }else{
    if (input_weights != ""){
      output_glm <- glm(formula = eval(input_formula),
                        data = input_table,
                        weights = input_table[,input_weights],
                        family = eval(parse(text = paste(input_family,"(",input_link,")"))))
    }else{
      output_glm <- glm(formula = eval(input_formula),
                        data = input_table,
                        family = eval(parse(text = paste(input_family,"(",input_link,")"))))
      
    }
  }
}}

print(capture.output(print(summary(output_glm))),collapse='\n') # display results
print(capture.output(print(confint(output_glm))),collapse='\n') # 95% CI for the coefficients

output_fittedvalues <- predict(output_glm, type="response") # predicted values
output_residuals <- residuals(output_glm, type="deviance") # residuals

output_glmdiag <- glm.diag(output_glm)
glm.diag.plots(output_glm,output_glmdiag) #Diagnostics plots for generalized linear models
