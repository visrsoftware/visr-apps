usePackage<-function (pkg) {
if (!require(pkg, character.only = TRUE)) 
  {install.packages(pkg, repos = "http://cran.us.r-project.org");
   require(pkg, character.only = TRUE)
  }
}

usePackage("gplots")

input_table <- mtcars
### 2D plot
#input_x <- "mpg"
#input_y <- "hp"
#input_color <- "cyl"
#plot(input_table[,input_x], input_table[,input_y], col=input_table[,input_color])

### kmeans
input_columns <- c("mpg", "hp")
input_k <- 3
input_algorithm <- c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen")[1]
input_itermax <- 10 #the maximum number of iterations allowed.
input_nstart <- 1 # how many random sets should be chosen?

visr.applyParameters()
{{
output_cluster <- kmeans(x = input_table[,input_columns], 
                        centers=input_k,
                        nstart=input_nstart,
                        algorithm=input_algorithm,
                        iter.max = input_itermax
                        )$cluster
}}

