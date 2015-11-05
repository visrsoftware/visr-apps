

usePackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {install.packages(pkg, repos = "http://cran.us.r-project.org");require(pkg, character.only = TRUE)}}

usePackage("car")


visr.applyParameters()


if(input_xlab == ""){input_xlab <- input_column}

{{
  if(input_group == ""){
  densityPlot(input_table[,input_column], 
              bw="SJ", 
              adjust=1, 
              kernel= input_kernel,
              grid= input_grid,
              show.bw = input_showbw, 
              rug = input_rug,
              xlab = input_xlab,
              ylab = input_ylab)
title(main = input_main)  
}else{

group <- as.factor(input_table[,input_group])
if(input_legendtitle == ""){input_legendtitle <- input_group}

densityPlot(input_table[,input_column],group, 
            bw="SJ", 
            adjust=1, 
            kernel= input_kernel,
            grid= input_grid,
            legend.location = input_legendlocation,
            legend.title = input_legendtitle,
            show.bw = input_showbw, 
            rug = input_rug,
            main = input_main,
            xlab = input_xlab,
            ylab = input_ylab)
}
}}

