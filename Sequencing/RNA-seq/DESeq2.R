# rm(list=ls())
# assign("last.warning", NULL, envir = baseenv())
# setwd('/Users/hyounesy/Research/svn/Projects/GeneCalc/visr/srcR')
# options(repos='http://cran.rstudio.com/')
# .libPaths(c("/Users/hyounesy/VisRseq/RLibs", .libPaths()))

# input_table <- read.table("/Users/hyounesy/Research/Data/VisRseq/error/041616-Ref-seq.xls",sep="\t",header=TRUE)
# input_group1<-c("WT1count","WT2count","WT3count")
# input_group2<-c("KO1count","KO2count","KO3count")
# input_advanced<-TRUE
# input_columns<-c("WT1count","WT2count","WT3count","KO1count","KO2count","KO3count")
# input_prefactor1<-"1,1,1,2,2,2"
# input_prefactor2<-""
# input_interaction<-FALSE
# input_ruvseq<-TRUE
# input_ruvn<-5000
# input_contrast<-""
# input_fdr<-0.1
# input_n<-10
# input_trans<-"Variance stabilizing transformation"
# input_plot<-""
# input_cores<-1
# input_test<-"Wald"
# input_fitType<-"local"

source("visrutils.R")
visr.biocLite("DESeq2")
visr.biocLite("BiocParallel")
visr.library("RColorBrewer")
visr.library("gplots")

visr.applyParameters()

if (!input_advanced && length(input_group1) == 1) visr.message("There are no replicates in group 1.", type="warning")
if (!input_advanced && length(input_group2) == 1) visr.message("There are no replicates in group 2.", type="warning")

if (input_advanced) local_pre1 <- eval(parse(text = paste("c(",input_prefactor1,")")))
if (input_advanced && length(input_columns) != length(local_pre1)) visr.message("The number of samples is different from the number of treatment conditions specified in factor 1", type="error")
if (input_advanced && identical(local_pre1[duplicated(local_pre1)], numeric(0))) visr.message("There are no replicates in the data for factor 1. To proceed: choose method -> blind; sharing mode -> fit-only.", type="warning")

if (input_advanced && input_prefactor2 != "") local_pre2 <- eval(parse(text = paste("c(",input_prefactor2,")")))
if (input_advanced && input_prefactor2 != "" && length(input_columns) != length(local_pre2)) visr.message("The number of samples is different from the number of treatment conditions specified in factor 2", type="error")
if (input_advanced && input_prefactor2 != "" && identical(local_pre2[duplicated(local_pre2)], numeric(0))) visr.message("There are no replicates in the data for factor 2. To proceed: choose method -> blind; sharing mode -> fit-only.", type="warning")
####1 input ####

par(mfrow = c(2,2))

{{
  if (!input_advanced) {
    input_columns <- c(input_group1, input_group2)
    input_factor1 <- factor(c(rep(1,length(input_group1)), rep(2,length(input_group2))))
    colData = data.frame( row.names = colnames(input_table[,input_columns]), factorOne = input_factor1)
    dds <- DESeqDataSetFromMatrix(countData = ceiling(input_table[,input_columns]),
                                  colData = colData,
                                  design = ~ factorOne)
  }
}}

{{
  if (input_advanced) {
    if (input_prefactor2 == ""){
      input_factor1 = as.character(eval(parse(text = paste("c(",input_prefactor1,")"))))
      colData = data.frame( row.names = colnames(input_table[,input_columns]), factorOne = input_factor1)
      dds <- DESeqDataSetFromMatrix(countData = ceiling(input_table[,input_columns]),
                                    colData = colData,
                                    design = ~ factorOne)
    }else{
      input_factor1 = as.character(eval(parse(text = paste("c(",input_prefactor1,")"))))
      input_factor2 = as.character(eval(parse(text = paste("c(",input_prefactor2,")"))))
      colData = data.frame( row.names = colnames(input_table[,input_columns]), factorOne = input_factor1, factorTwo = input_factor2)
      if (input_interaction == FALSE){
        dds <- DESeqDataSetFromMatrix(countData = ceiling(input_table[,input_columns]),
                                      colData = colData,
                                      design = ~ factorOne + factorTwo)
      }else{
        dds <- DESeqDataSetFromMatrix(countData = ceiling(input_table[,input_columns]),
                                      colData = colData,
                                      design = ~ factorOne*factorTwo)
      }
    }
  }
}}


dds


####3. differential expression analysis ####
{{
if (input_cores == 1) {
  dds <- DESeq(dds, test = input_test, fitType = input_fitType) # the standard differntial expression analysis steps
}else{
  register(MulticoreParam(input_cores)) # for experiments with many samples (eg.100 samples) to speed up the analysis
  dds <- DESeq(dds,  test = input_test, fitType = input_fitType, parallel=TRUE)
}
}}

resultnames <- resultsNames(dds)


{{
if (input_contrast == "") {
  res <- results(dds, addMLE=TRUE)
}else{
  input_contrast2 = eval(parse(text = paste("c(",input_contrast,")")))
  res <- results(dds, contrast = input_contrast2)
}
}}

# results tables, ie. log2 fold changes, p values, adjusted p values
# add to the result the unshrunken maximum likelihood estimate (MLE) for the log2 fold change

{{
  if(input_ruvseq) {
    # To estimate the factors of unwanted variation, we need a set of negative control genes,
    # e.g., least significantly DE genes based on a first-pass DE analysis performed prior to RUVg normalization
    if (!input_advanced) {
      set <- newSeqExpressionSet(as.matrix(ceiling(input_table[,input_columns])),
                                 phenoData = data.frame(input_factor1, row.names=colnames(input_table[,input_columns])))
    }else{
      if (input_prefactor2 == ""){
        set <- newSeqExpressionSet(as.matrix(ceiling(input_table[,input_columns])),
                                   phenoData = data.frame(input_factor1, row.names=colnames(input_table[,input_columns])))
      }else{
        set <- newSeqExpressionSet(as.matrix(ceiling(input_table[,input_columns])),
                                   phenoData = data.frame(input_factor1, input_factor2, row.names=colnames(input_table[,input_columns])))
      }
    }

    res  = cbind(ceiling(input_table[,input_columns]), pval = res$pvalue, padj = res$padj)
    res2  = res[ order(res$padj), ]
    # we consider all but the top 5000 genes as ranked by edgeR p-values.
    ruvempirical <- rownames(set)[which(!(rownames(set) %in% rownames(res2)[1:input_ruvn]))]
    set2 <- RUVg(set, ruvempirical, k=1)
    pData(set2)

    if (!input_advanced) {
      colData = data.frame( row.names = colnames(input_table[,input_columns]), factorOne = input_factor1, variation = pData(set2)$W_1)
      dds <- DESeqDataSetFromMatrix(countData = ceiling(input_table[,input_columns]),
                                    colData = colData,
                                    design = ~ factorOne + variation)
    }else{
      if (input_prefactor2 == ""){
        colData = data.frame( row.names = colnames(input_table[,input_columns]), factorOne = input_factor1, variation = pData(set2)$W_1)
        dds <- DESeqDataSetFromMatrix(countData = ceiling(input_table[,input_columns]),
                                      colData = colData,
                                      design = ~ factorOne  + variation)
      }else{
        colData = data.frame( row.names = colnames(input_table[,input_columns]), factorOne = input_factor1, factorTwo = input_factor2, variation = pData(set2)$W_1)
        if (input_interaction == FALSE){
          dds <- DESeqDataSetFromMatrix(countData = ceiling(input_table[,input_columns]),
                                        colData = colData,
                                        design = ~ factorOne + factorTwo + variation)
        }else{
          dds <- DESeqDataSetFromMatrix(countData = ceiling(input_table[,input_columns]),
                                        colData = colData,
                                        design = ~ factorOne*factorTwo + variation)
        }
      }
    }
    if (input_cores == 1) {
      dds <- DESeq(dds, test = input_test, fitType = input_fitType) # the standard differntial expression analysis steps
    }else{
      register(MulticoreParam(input_cores)) # for experiments with many samples (eg.100 samples) to speed up the analysis
      dds <- DESeq(dds,  test = input_test, fitType = input_fitType, parallel=TRUE)
    }
    resultnames <- resultsNames(dds)
    if (input_contrast == "") {
      res <- results(dds, addMLE=TRUE)
    }else{
      input_contrast2 = eval(parse(text = paste("c(",input_contrast,")")))
      res <- results(dds, contrast = input_contrast2)
    }

  }
}}



resOrdered <- res[order(res$padj),] # order results table by smallest adjusted p value
head(resOrdered)
summaryres <- summary(res) # summarize some basic tallies


## dispersion plot
#The dispersion estimate plot shows the gene-wise estimates (black), the
#fitted values (red), and the final maximum a posteriori estimates used in testing (blue).
plotDispEsts(dds)

#####4 explore the result #####
## Plot of counts for one gene with the smallest p value
#plotCounts(dds, gene=which.min(res$padj), intgroup="factorOne")
## info about which variables and tests were used
resDes <- mcols(res)$description
## All row-wise calculated values
rowWise <- mcols(dds,use.names=TRUE)
rowWiseDes <- mcols(mcols(dds), use.names=TRUE)


# resSig <- subset(resOrdered, padj < 0.1)

resSig = subset(resOrdered, padj < input_fdr)  # filter for significant genes, according to some chosen threshold for the false dicovery rate (FDR),
mostsig <- head( resSig[ order(resSig$padj ), ] ,input_n) # the most significantly differentially expressed genes
mostdown <- head( resSig[ order( resSig$log2FoldChange, -resSig$baseMean ), ] ,input_n) # the most strongly down-regulated of the significant genes
mostup <- head( resSig[ order( -resSig$log2FoldChange, -resSig$baseMean ), ],input_n ) # the most strongly up-regulated ones

output_baseMean <- res$baseMean
output_log2FoldChange <- res$log2FoldChange
output_lfcSE <- res$lfcSE
output_stat <- res$stat
output_pvalue <- res$pvalue
output_padj <- res$padj
output_lfcMLE <- res$lfcMLE #  "unshrunken" maximum likelihood estimates (MLE) of log2 fold change
output_cluster <- ifelse(output_padj < input_fdr, ifelse(output_log2FoldChange < 0, -1, 1), 0)
output_cluster[is.na(output_cluster)] = 0

# Histogram of p-values
hist(res$padj, breaks=100, col="skyblue", border="slateblue", main="")

## MA plot
plotMA(res, main="DESeq2",alpha = input_fdr, ylim=c(-2,2)) # the log2 fold changes attributable to a given variable over the mean of normalized counts
abline(h=c(-.5,.5),col="dodgerblue",lwd=2)

## barplot
counts<-prop.table(table((output_cluster)))*length(output_cluster)
bplt<-barplot(counts, main="DESeq2 Summary")
text(y= counts, x= bplt, labels=as.character(counts), xpd=TRUE, adj=c(0.5,0))



#####2 data explore #####
{{
  if (input_trans == "Regularized log transformation"){
    rld <- rlog(dds)
    rlogMat <- assay(rld)
    # The regularized logarithm or rlog incorporates a prior on the sample differences [1],
    # and the other uses the concept of variance stabilizing transformations (VST)
  }else{
    vsd <- varianceStabilizingTransformation(dds,fitType="local")
    vstMat <- assay(vsd)
  }
}}


select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

## heatmap of the count matrix
{{
  if (input_plot == "heatmap for the 30 most highly expressed genes") {
    if (input_trans == "Regularized log transformation"){
      heatmap.2(assay(rld)[select,], col = hmcol,
                Rowv = FALSE, Colv = FALSE, scale="none",
                dendrogram="none", trace="none", margin=c(10, 6))
    }else{
      heatmap.2(assay(vsd)[select,], col = hmcol,
                Rowv = FALSE, Colv = FALSE, scale="none",
                dendrogram="none", trace="none", margin=c(10, 6))
    }
  }
}}
## heatmap of the sample to sample distance
{{
  if (input_plot == "Heatmap of the sample-to-sample distances"){
    if (input_trans == "Regularized log transformation"){
      distsRL <- dist(t(assay(rld)))
    }else{
      distsRL <- dist(t(assay(vsd)))
    }

    mat <- as.matrix(distsRL)

    if (input_prefactor2 == ""){
      rownames(mat) <- colnames(mat) <- with(colData(dds), paste(factorOne, sep=" : "))
    }else{
      rownames(mat) <- colnames(mat) <- with(colData(dds), paste(factorOne, factorTwo, sep=" : "))
    }

    hc <- hclust(distsRL)
    heatmap.2(mat, Rowv=as.dendrogram(hc),
              symm=TRUE, trace="none",
              col = rev(hmcol), margin=c(13, 13))
  }
}}
## pca plot of the samples
{{
  if (input_plot == "Principal component plot of the samples"){
    if (input_prefactor2 == ""){
      if (input_trans == "Regularized log transformation"){
        plotPCA(rld,intgroup=c("factorOne"))
      }else{
        plotPCA(vsd,intgroup=c("factorOne"))
      }
    }else{
      if (input_trans == "Regularized log transformation"){
        plotPCA(rld,intgroup=c("factorOne", "factorTwo"))
      }else{
        plotPCA(vsd, intgroup=c("factorOne", "factorTwo"))
      }
    }
  }
}}


print(capture.output(cat("\n1.input\n")),collapse='\n')
print(capture.output(cat("\ndata\n")),collapse='\n')
print(capture.output(print(colData)),collapse='\n')

print(capture.output(cat("\n2.explore the data\n")),collapse='\n')

print(capture.output(cat("\n3.differential expression analysis\n")),collapse='\n')
print(capture.output(cat("\nterms in the design\n")),collapse='\n')
print(capture.output(print(resultnames)),collapse='\n')
print(capture.output(cat("\nsummary of basic tallies\n")),collapse='\n')
print(capture.output(print(summary(res))),collapse='\n')

print(capture.output(cat("\n4.explore the result\n")),collapse='\n')
print(capture.output(cat("\ninfo about which variables and tests were used \n")),collapse='\n')
print(capture.output(print(resDes)),collapse='\n')
print(capture.output(cat("\nthe most significantly differentially expressed genes:\n")),collapse='\n')
print(capture.output(print(mostsig)),collapse='\n')
print(capture.output(cat("\nthe most strongly down-regulated of the significant genes:\n")),collapse='\n')
print(capture.output(print(mostdown)),collapse='\n')
print(capture.output(cat("\nthe most strongly up-regulated ones:\n")),collapse='\n')
print(capture.output(print(mostup)),collapse='\n')

print(capture.output(print(counts)),collapse='\n')
