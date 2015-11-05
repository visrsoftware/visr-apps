# http://homepages.inf.ed.ac.uk/bwebb/statistics/ANOVA_in_R.pdf
# http://www.cookbook-r.com/Statistical_analysis/ANOVA/#problem
# http://www.r-tutor.com/category/statistical-concept/anova
# http://personality-project.org/r/r.anova.html
# http://www.statmethods.net/stats/anova.html

Gender <- c("Female", "Female", "Female", "Male", "Male", "Male") 
Income <- c("High", "High", "Medium", "Medium", "Low", "Low")
Believer <- c("Yes", "No", "Yes", "No", "Yes", "No")
Count <- c(435, 147, 375, 134, 257, 234)
afterlife <- data.frame(Gender, Income, Believer, Count)
rm(Gender, Income, Believer, Count)


input_table <- mtcars
input_y <- "am"
input_x1 <- "cyl"
input_x2 <- "disp"
input_x3 <- "hp"

input_weights <- "drat"

