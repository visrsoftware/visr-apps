usePackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {install.packages(pkg, repos = "http://cran.us.r-project.org");require(pkg, character.only = TRUE)}}

usePackage("scatterplot3d")


visr.applyParameters()

# if(input_x = "") error_message <- "variable for horizontal axis is missing."
# if(input_y = "") error_message <- "variable for vertical axis is missing."
# if(input_z = "") error_message <- "variable for out-of-screen axis is missing."

if (input_xlab == "") {input_xlab <- input_x}
if (input_ylab == "") {input_ylab <- input_y}
if (input_zlab == "") {input_zlab <- input_z}
if (input_main == "") {input_main <- "3D Scatterplot"}

{{
  output_3ds <-scatterplot3d(input_table[,input_x],input_table[,input_y],input_table[,input_z], 
                           xlab = input_xlab,ylab = input_ylab,zlab = input_zlab,
                           type = input_type,
                           #col = input_col,
                           pch=input_pch, 
                           highlight.3d=input_highlight3d,
                           scale.y = input_scaley,
                           angle = input_angle,
                           grid = input_grid,
                           box = input_box,
                           main=input_main, sub = input_sub)
}}
{{
  if (input_fit == TRUE){
  fit <- lm(input_table[,input_z]~ input_table[,input_x] + input_table[,input_y]) 
  output_3ds$plane3d(fit)
}
}}
