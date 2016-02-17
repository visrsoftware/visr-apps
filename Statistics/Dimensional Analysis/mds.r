source("visrutils.R")
visr.applyParameters()

if (input_main == "")  input_main = "Multidimensional Scaling"

mdsMatrix <- input_table[, input_columns]

if (input_method != "hamming") {
  d <- dist(mdsMatrix, method = input_method) # distances between the rows
} else {
  n <- nrow(mdsMatrix)
  d <- matrix(nrow=n, ncol=n)
  for(i in seq_len(n))
    for(j in seq(i, n))
      d[j, i] <- d[i, j] <- sum(mdsMatrix[i,] != mdsMatrix[j,])
}

fit <- cmdscale(d, eig = TRUE, k = input_k) # k is the number of dim
output_x <- fit$points[,1]
output_y <- fit$points[,2]

if (input_labels != "") {
  plot(output_x, output_y, xlab = input_xlab, ylab = input_ylab, main = input_main, type="n")
  text(output_x, output_y, labels = input_table[,input_labels], cex = input_cex, pos = input_pos, col = input_col)
} else {
  plot(output_x, output_y, xlab = input_xlab, ylab = input_ylab, main = input_main)
}
