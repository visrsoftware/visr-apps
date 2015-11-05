biocPackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {source("http://bioconductor.org/biocLite.R"); biocLite(pkg);require(pkg, character.only = TRUE)}}
biocPackage("baySeq")
usePackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {install.packages(pkg, repos = "http://cran.us.r-project.org");require(pkg, character.only = TRUE)}}
usePackage("parallel")

# input_table <- read.csv("C:/Users/hexu/Desktop/VisRSeq/VisRSeq/data/counts_tab.txt",sep = "\t")
# input_columns <- c("CT.PA.1","CT.PA.2","KD.PA.3","KD.PA.4") #,"CT.SI.5","KD.SI.6","CT.SI.7")
# #input_sample1 <- c("CT.PA.1","CT.PA.2")
# #input_sample2 <- c("KD.PA.3","KD.PA.4")
# input_groups <- "1,1,2,2"
# input_groups0 <- "1,1,1,1"
# 
# data(simData)
# input_table <- simData
# colnames(input_table)<-c("a1","a2","a3","a4","a5","b1","b2","b3","b4","b5")
# input_columns <- c("a1","a2","a3","a4","b1","b2","b3","b4")
# input_groups <- "1,1,1,1,2,2,2,2"
# input_groups0 <- "1,1,1,1,1,1,1,1"
# #or
# data(pairData)
# input_table <- pairData
# # The first four columns in these data are paired with the second four columns
# colnames(input_table)<-c("a1","a2","a3","a4","b1","b2","b3","b4")
# input_columns1 <- c("a1","a2","a3","a4")
# input_columns2 <- c("b1","b2","b3","b4")
# input_groups <- "1,1,2,2"
# input_groups0 <- "1,1,1,1"

input_normaliseData <- TRUE
input_samplesize <- 1000
input_equalDispersions <- TRUE
input_estimation <- c("QL","ML","edgeR")[1]
#input_zeroML <- FALSE
input_consensus <- FALSE

input_pET <- c("BIC","none","iteratively")[1]

input_number <- 10

visr.applyParameters()

{{
  if(input_paired == FALSE) {
    if(length(input_columns) != length(eval(parse(text = paste("c(",input_groups,")"))))) error_message<-"the number of samples is different from the number of treatment conditions specified in DE"
    if(length(input_columns) != length(eval(parse(text = paste("c(",input_groups0,")"))))) error_message<-"the number of samples is different from the number of treatment conditions specified in NDE"
  }else{
    if(length(input_columns1) != length(eval(parse(text = paste("c(",input_groups,")"))))) error_message<-"the number of samples in samples1 is different from the number of treatment conditions specified in DE"
    if(length(input_columns1) != length(eval(parse(text = paste("c(",input_groups0,")"))))) error_message<-"the number of samples in samples1 is different from the number of treatment conditions specified in NDE"
    if(length(input_columns2) != length(eval(parse(text = paste("c(",input_groups,")"))))) error_message<-"the number of samples in samples2 is different from the number of treatment conditions specified in DE"
    if(length(input_columns2) != length(eval(parse(text = paste("c(",input_groups0,")"))))) error_message<-"the number of samples in samples2 is different from the number of treatment conditions specified in NDE"    
  }
}}


top1 <- NULL
cl <- makeCluster(8) 

{{
  if (input_paired == FALSE){
    #par(mfrow = c(1,2))
    replicates <- as.character(eval(parse(text = paste("c(",input_groups,")"))))
    groups <- list(NDE = eval(parse(text = paste("c(",input_groups0,")"))),DE = eval(parse(text = paste("c(",input_groups,")"))))
    CD <- new("countData", 
              data = as.matrix(input_table[,input_columns]), 
              replicates = replicates, 
              groups = groups)
    #     CD <- new("countData", data = data)
    #     groups(CD) <- groups
    #     replicates(CD) <- replicates
    libsizes(CD) <- getLibsizes(CD) # library size
    plotMA.CD(CD, samplesA = "1", samplesB = "2", normaliseData = input_normaliseData)
    CD <- getPriors.NB(CD,  # estimate an empirical distribution by bootstrapping
                       samplesize = input_samplesize, 
                       equalDispersions = input_equalDispersions,
                       estimation = input_estimation, 
                       #zeroML = input_zeroML,
                       consensus = input_consensus,
                       cl = cl)
    CD <- getLikelihoods(CD, pET = input_pET, cl = cl) # posterior likelihoods
    # estimating the proportions of differentially expressed counts
    CD@estProps
    # the (log) posterior likelihoods that each count belongs to a group defined by the @group slot of the input object
    output_posterior1 <- CD@posteriors[,1]
    output_posterior2 <- CD@posteriors[,2]
    
    #  the proportion of differential expressed counts
    CD@estProps[2]
    
    #  the top candidates for differential expression
    top <- topCounts(CD, group = "DE", number = input_number)
    # Takes posterior likelihoods and returns the counts with highest (or lowest) likelihood of association
    output_table <- topCounts(CD, group = "DE", number = (dim(input_table[,input_columns])[1]))
    
    # plot the posterior likelihoods against the log-ratios of the two sets of samples
    # (where this would be non-infinite) or log
    # values (where all data in the other sample group consists of zeros).
    #plotPosteriors(CD, group = "DE")
    
  }else{
    replicates <- eval(parse(text = paste("c(",input_groups,")")))
    groups <- list(NDE = eval(parse(text = paste("c(",input_groups0,")"))),DE = eval(parse(text = paste("c(",input_groups,")"))))
    CD <- new("countData", 
              data = list(as.matrix(input_table[,input_columns1]), as.matrix(input_table[,input_columns2])),
              replicates = replicates,
              groups = groups,
              densityFunction = bbDensity)
    libsizes(CD) <- getLibsizes(CD) # library size
    #     plotMA.CD(CD, samplesA = "1", samplesB = "2", normaliseData = input_normaliseData)
    CD <- getPriors(CD, samplesize = input_samplesize, cl = cl)
    CD <- getLikelihoods(CD, pET = input_pET, nullData = TRUE, cl = cl)
    # The use of 'nullData = TRUE' in this context allows us to identify
    # pairs which show no differential expression between replicate groups, but does show deviation from a one-to-one ratio
    # of data between pairs.
    
    # the (log) posterior likelihoods that each count belongs to a group defined by the @group slot of the input object
    output_posterior1 <- CD@posteriors[,1]
    output_posterior2 <- CD@posteriors[,2]
    top <- topCounts(CD, group = 2, number = input_number)
    top1 <- topCounts(CD, group = 1, number = input_number)
    # top candidates for differential expression between replicate groups
    output_table1 <- topCounts(CD, group = 2,number = (dim(input_table[,input_columns1])[1]))
    # consistent differential expression between the pairs
    output_table2 <- topCounts(CD, group = 1,number = (dim(input_table[,input_columns1])[1]))
    
  }
}}

print(capture.output(cat("\nsize factors\n")),collapse='\n')

# print(capture.output(cat("\nestimating the proportions of differentially expressed counts\n")),collapse='\n')
# print(capture.output(print(CD@estProps)),collapse='\n')

print(capture.output(cat("\nthe top candidates for differential expression\n")),collapse='\n')
print(capture.output(print(top)),collapse='\n')

print(capture.output(cat("\nfor paired data, consistent differential expression between the pairs\n")),collapse='\n')
print(capture.output(print(top1)),collapse='\n')

