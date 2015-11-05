biocPackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {source("http://bioconductor.org/biocLite.R"); biocLite(pkg);require(pkg, character.only = TRUE)}}
biocPackage("NOISeq")

# data(Marioni)
# input_table <- read.csv("C:/Users/hexu/Desktop/old versions/VisRSeq.0.72.4/VisRSeq/data/counts_tab.txt",sep = "\t")
# # input_columns <- c("CT.PA.1","CT.PA.2","KD.PA.3","KD.PA.4") #,"CT.SI.5","KD.SI.6","CT.SI.7")
# input_sample1 <- c("CT.PA.1","CT.PA.2")
# input_sample2 <- c("KD.PA.3","KD.PA.4")
# #input_groups <- "1,1,2,2"
# #input_groups0 <- "1,1,1,1"
# #input_pregroup1 <- "1,2,1,2,2,1,2,1,2,1"
# #input_pregroup2 <- "1,2,1,2,2,1,2,3,4,3"

# input_replicates <- c("technical","biological","no")[1]
# input_k <- 0.5
# input_norm <- c("rpkm","uqua","tmm","n")[1]
# input_lc <- 0
# input_pnr <- 0.2
# input_nss <- 5
# input_v <- 0.02
# 
# input_r <- 50
# input_adj <- 1.5
# input_nclust <- 15
# input_filter <- 1
# input_cv <- 500
# input_cpm <- 1
# 
# input_q <- 0.8

visr.applyParameters()

par(mfrow = c(1,2))

data1 = input_table[,input_sample1]
data2 = input_table[,input_sample2]
mycounts = cbind(data1,data2)
input_group = c(rep("A",length(input_sample1)),rep("B",length(input_sample2)))

myfactors = data.frame(group = input_group)

mydata <- readData(data=mycounts,factors=myfactors)
mydata

# # Global saturation plot to compare two samples, respectively.
# mysaturation = dat(mydata, k = 0, ndepth = 7, type = "saturation")
# explo.plot(mysaturation, toplot = 1, samples = 1:2, yleftlim = NULL, yrightlim = NULL)
# 
# # Number of features with low counts for each sample.
# explo.plot(mycountsbio, toplot = 1, samples = NULL, plottype = "barplot")
# 
# # RNA composition plot 
# mycd = dat(mydata, type = "cd", norm = FALSE, refColumn = 1)
# explo.plot(mycd)


# # normalization (choose one of three methods)
# myRPKM = rpkm(assayData(mydata)$exprs, long = mylength, k = 0, lc = 1)
# myUQUA = uqua(assayData(mydata)$exprs, long = mylength, lc = 0.5, k = 0)
# myTMM = tmm(assayData(mydata)$exprs, long = 1000, lc = 0)
# head(myRPKM[,1:4])


# # filtering
# myfilt = filtered.data(mycounts, factor = myfactors$Tissue, norm = FALSE, depth = NULL, method = 1, cv.cutoff = 100, cpm = 1)

# results

{{
  if (input_replicates != "biological") {
    # technical replicates or no replicates at all
    mynoiseq = noiseq(mydata, 
                      factor="group",
                      replicates = input_replicates,
                      k = input_k, 
                      norm = input_norm,  
                      lc = input_lc,
                      pnr = input_pnr, 
                      nss = input_nss, 
                      v = input_v)
  }else{  
    # data with biological replicates
    mynoiseq = noiseqbio(mydata, 
                         k = input_k, 
                         norm = input_norm,  
                         factor="group",
                         lc = input_lc, 
                         r = input_r, 
                         adj = input_adj, 
                         nclust = input_nclust,
                         plot = FALSE,
                         a0per = 0.9, 
                         random.seed = 12345, 
                         filter = input_filter,
                         cv.cutoff = input_cv, 
                         cpm = input_cpm)
  }
}}

# view results
# head(mynoiseq@results[[1]])
print(capture.output(cat("\nresult\n")),collapse='\n')
print(capture.output(print(head(mynoiseq@results[[1]]))),collapse='\n')

# all the differentially expressed features
deg = degenes(mynoiseq, q = input_q, M = NULL)
# only the differentially expressed features that are more expressed in condition 1 than in condition 2 (M = "up") 
deg1 = degenes(mynoiseq, q = input_q, M = "up")
# only the differentially expressed features that are under-expressed in condition 1 with regard to condition 2 (M = "down"):
deg2 = degenes(mynoiseq, q = input_q, M = "down")

output_cluster<-rep(0, dim(mycounts)[1])
output_cluster[as.numeric(rownames(deg1))]<-1
output_cluster[as.numeric(rownames(deg2))]<- -1
counts<-prop.table(table((output_cluster)))*length(output_cluster)

print(capture.output(cat("\nsummary\n")),collapse='\n') 
print(capture.output(print(counts)),collapse='\n')

# Summary plot of the expression values for both conditions (black), where differentially expressed genes are highlighted (red)
DE.plot(mynoiseq, q = input_q, graphic = "expr", log.scale = TRUE, main = "Expression plot")

# M is the log-fold-change and D is the absolute value of the difference between conditions
DE.plot(mynoiseq, q = input_q, graphic = "MD", main = "MD plot")

# Manhattan plot
# DE.plot(mynoiseq, chromosomes = c(1,2), log.scale = TRUE,join = FALSE, q = input_q, graphic = "chrom", main = "Manhattan plot")

# Distribution of differentially expressed features per chromosomes or biotypes
# DE.plot(mynoiseq, chromosomes = NULL, q = input_q, graphic = "distr")

