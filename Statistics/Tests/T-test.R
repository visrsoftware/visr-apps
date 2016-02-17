
input_table <- mtcars
input_x <- "mpg"
input_y <- ""
input_group <- "vs"
input_alternative <- c("two.sided", "less", "greater")[1] # the alternative hypothesis
input_mu <- 0 # a number indicating the true value of the mean (or difference in means if you are performing a two sample test)
input_paired <- FALSE # a logical indicating whether you want a paired t-test
input_conflevel <- 0.95 # confidence level of the interval

visr.applyParameters()

{{
  if (input_y == ""){
    if (input_group == ""){
      t.test(input_table[,input_x], 
             alternative = input_alternative,
             mu = input_mu,
             paired = input_paired,
             conf.level = input_conflevel) 
    }else{
      t.test(input_table[,input_x] ~ input_table[,input_group],
             data = input_table,
             alternative = input_alternative,
             mu = input_mu,
             paired = input_paired,
             conf.level = input_conflevel)
    }
  }else{
    t.test(input_table[,input_x], 
           input_table[,input_y], 
           alternative = input_alternative,
           mu = input_mu,
           paired = input_paired,
           conf.level = input_conflevel)
  }  
}}

