usePackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {install.packages(pkg, repos = "http://cran.us.r-project.org");require(pkg, character.only = TRUE)}}
usePackage("ggplot2")

# input_table <- read.csv("C:/Users/hexu/Desktop/VisRSeq.0.72.4/VisRSeq/data/result.txt",sep = "\t")

# input_logfold <- "log2FoldChange"
# input_p <- "pvalue"
# #input_padj <- "padj"
# input_label <- "Gene"
# input_cutlog <- 2
# input_cutp <- 0.05
# input_highlight <- c("Bonferroni","FDR")[1]
# 
# input_text = FALSE
# input_alpha <- 0.4
# input_size <- 1.75
# input_textsize <- 1.2

visr.applyParameters()

input_n <- length(input_table[,input_logfold])
{{
  if (input_highlight == "Bonferroni") {
    ##Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
    input_table$threshold = as.factor(abs(input_table[,input_logfold]) > input_cutlog & input_table[,input_p] < 0.05/input_n)
  }else{
    # use the false discovery rate as a cut-off (FDR)
    input_table$threshold = as.factor(abs(input_table[,input_logfold]) > input_cutlog & input_table[,input_p] < 0.05)
  }
}}
##Construct the plot object
{{
  g = ggplot(data=input_table, 
             aes(x=input_table[,input_logfold], 
                 y=-log10(input_table[,input_p]), 
                 colour=threshold)) +
    geom_point(alpha=input_alpha, size=input_size) +
    xlim(c(-10, 10)) + 
    ylim(c(0, 15)) +
    xlab("log2 fold change") + 
    ylab("-log10 p-value")
}}

##Add text to the plot
{{
  if (input_text == TRUE) {
    g = g + geom_text(aes(x=input_table[,input_logfold], 
                          y=-log10(input_table[,input_p]),
                          label=input_table[,input_label], 
                          size=input_textsize), 
                      colour="black")
  }
}}

print(g)
