usePackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {install.packages(pkg, repos = "http://cran.us.r-project.org");require(pkg, character.only = TRUE)}}

usePackage("ggplot2")
usePackage("reshape")


visr.applyParameters()


if(input_xlab == ""){input_xlab <- NULL}
if(input_ylab == ""){input_ylab <- NULL}

{{
  if(length(input_column) != 1) {
    # "Melt" data:
  data <- melt(input_table[,input_column])
  p <- ggplot(data) + 
    geom_density(aes(x = value, colour = variable,fill= variable), alpha= input_transparency) + 
    labs(title = input_main,x = input_xlab, y = input_ylab) +
    theme(legend.position = input_legendlocation)
  }else{
    data <- melt(input_table[,input_column])
    p <- ggplot(data) + 
      geom_density(aes(x = value), alpha= input_transparency) + 
      labs(title = input_main,x = input_xlab, y = input_ylab) +
      theme(legend.position = input_legendlocation)
  }
}}
print(p)
