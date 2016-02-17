
visr.applyParameters()

if (input_row == 0) {input_row <- dim(input_table)[1]}

set.seed(input_seed)

if (input_choice == "Binomial distribution") {output_sample <- as.data.frame(matrix(rbinom(input_row*1, size=input_size1, prob=input_prob1), ncol=1))}

if (input_choice == "Poisson distribution") {output_sample <- as.data.frame(matrix(rpois(input_row*1, lambda=input_lambda), ncol=1))}

if (input_choice == "Geometric distribution") {output_sample <- as.data.frame(matrix(rgeom(input_row*1, prob=input_prob2), ncol=1))}

if (input_choice == "Hypergeometric distribution") {output_sample <- as.data.frame(matrix(rhyper(input_row*1, m=input_m, n=input_n, k=input_K),ncol=1))}

if (input_choice == "Negative binomial distribution") {output_sample <- as.data.frame(matrix(rnbinom(input_row*1, size=input_size2, prob=input_prob3), ncol=1))}

output_col <- output_sample[,1]

output_table <- output_sample

hist(output_col,xlab = "sample", main = "Histogram of the sample generated")
