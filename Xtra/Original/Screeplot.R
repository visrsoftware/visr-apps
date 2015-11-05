# http://courses.ttu.edu/isqs6348-westfall/images/6348/Exploratory_Factor_Analysis_in_R.pdf\

usePackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {install.packages(pkg, repos = "http://cran.us.r-project.org");require(pkg, character.only = TRUE)}}

usePackage("psych")

input_table <- mtcars # a correlation matrix or a data matrix. If data, then correlations are found using pairwise deletions.
input_columns <- c("mpg","cyl","disp", "hp","drat")
input_factors <- TRUE#  draw the scree for factors
input_pc <- TRUE #  draw the scree for components
input_hline <- NULL # if null, draw a horizontal line at 1, otherwise draw it at hline (make negative to not draw it)
input_main <- "Scree plot"
input_add <- FALSE # Should multiple plots be drawn?

visr.applyParameters()

{{
  scree(input_table[,input_columns],
        factors = input_factors,
        pc = input_pc,
        hline = input_hline,
        main = input_main,
        add = input_add)
}}

