
visr.applyParameters()

input_col <- 1

if (input_row == 0) {input_row <- dim(input_table)[1]}

if (input_choice == "Normal distribution") {output_sample <- as.data.frame(matrix(rnorm(input_row*input_col, mean=0, sd=1), ncol=input_col))}

if(input_choice == "t distribution"){output_sample <- as.data.frame(matrix(rt(input_row*input_col, df=input_df), ncol=input_col))}

if(input_choice == "Chi-squared distribution"){output_sample <- as.data.frame(matrix(rchisq(input_row*input_col, df=input_df), ncol=input_col))}

if(input_choice == "F distribution"){output_sample <- as.data.frame(matrix(rf(input_row*input_col, df1=input_df1, df2=input_df2), ncol=input_col))}

if(input_choice == "Exponential distribution"){output_sample <- as.data.frame(matrix(rexp(input_row*input_col, rate=input_rate), ncol=input_col))}

if(input_choice == "Uniform distribution"){output_sample <- as.data.frame(matrix(runif(input_row*input_col, min=input_min, max=input_max), ncol=input_col))}

if(input_choice == "Beta distribution"){output_sample <- as.data.frame(matrix(rbeta(input_row*input_col, shape1=input_shape1, shape2=input_shape2),ncol=input_col))}

if(input_choice == "Cauchy distribution"){output_sample <- as.data.frame(matrix(rcauchy(input_row*input_col, location=input_location1, scale=input_scale1), ncol=input_col))}

if(input_choice == "Logistic distribution"){output_sample <- as.data.frame(matrix(rlogis(input_row*input_col, location=input_location2, scale=input_scale2), ncol=input_col))}

if(input_choice == "Lognormal distribution"){output_sample <- as.data.frame(matrix(rlnorm(input_row*input_col, meanlog=input_meanlog, sdlog=input_sdlog),  ncol=input_col))}

if(input_choice == "Gamma distribution"){output_sample <- as.data.frame(matrix(rgamma(input_row*input_col, shape=input_shape3, scale=input_scale3),  ncol=input_col))}

if(input_choice == "Weibull distribution"){output_sample <- as.data.frame(matrix(rweibull(input_row*input_col, shape=input_shape4, scale=input_scale4),  ncol=input_col))}

if(input_choice == "Gumbel distribution"){output_sample <- as.data.frame(matrix(log(rweibull(input_row*input_col, shape=input_shape5, scale=input_scale5)),ncol=input_col))}

output_col <- output_sample[,1]

output_table <- output_sample
hist(output_col,xlab = "sample", main = "Histogram of the sample generated")
