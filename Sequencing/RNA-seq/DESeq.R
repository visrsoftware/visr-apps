# rm(list=ls())
# assign("last.warning", NULL, envir = baseenv())
# setwd('/Users/hyounesy/Research/svn/Projects/GeneCalc/visr/srcR')
# options(repos='http://cran.rstudio.com/')
# .libPaths(c("/Users/hyounesy/VisRseq/RLibs", .libPaths()))

# input_table <- read.table("/Users/hyounesy/Research/Data/VisRseq/error/041616-Ref-seq.xls",sep="\t",header=TRUE)
# input_group1<-c("WT1count","WT2count","WT3count")
# input_group2<-c("KO1count","KO2count","KO3count")
# input_advanced<-FALSE
# input_columns<-c("WT1count","WT2count","WT3count","KO1count","KO2count","KO3count")
# input_GLM<-FALSE
# input_prefactor1<-"1,1,1,2,2,2"
# input_prefactor2<-""
# input_ruvseq<-TRUE
# input_ruvn<-5000
# input_fdr<-0.1
# input_dispersionmethod<-"pooled"
# input_dispersionsharingMode<-"maximum"
# input_dispersionfitType<-"local"
# input_plot<-""
# input_plotx<-"baseMean"
# input_ploty<-"foldChange"
# input_plotlog<-"xy"

source("visrutils.R")
visr.biocLite("DESeq")
visr.biocLite("RUVSeq")
visr.library("RColorBrewer")
visr.library("gplots")

visr.applyParameters()

if (!input_advanced && length(input_group1) == 1) {{ 
  visr.message("There are no replicates in group 1.", type="warning")}}

if (!input_advanced && length(input_group2) == 1) {{
  visr.message("There are no replicates in group 2.", type="warning")}}

if ((length(input_group1) == 1 || length(input_group2) == 1) && input_dispersionmethod != "blind") {{
  visr.message("None of your conditions is replicated.\nUse method='blind' to estimate across conditions, or 'pooled-CR', if you have crossed factors.", type="error")}}

if ((length(input_group1) == 1 || length(input_group2) == 1) && input_dispersionsharingMode!="fit-only") {{
  visr.message("None of your conditions is replicated.\nUse sharing-model='fit-only'.", type="error")}}

if (input_advanced) {{
  local_pre1 <- eval(parse(text = paste("c(",input_prefactor1,")")))}}

if (input_advanced && length(input_columns) != length(local_pre1)) {{
  visr.message("The number of samples is different from the number of treatment conditions specified in factor 1", type="error")}}

if (input_advanced && identical(local_pre1[duplicated(local_pre1)], numeric(0))) {{
  visr.message("There are no replicates in the data for factor 1. To proceed: choose method -> blind; sharing mode -> fit-only.", type="warning")}}

if (input_advanced && input_GLM) {{
  local_pre2 <- eval(parse(text = paste("c(",input_prefactor2,")")))}}

if (input_advanced && input_GLM && length(input_columns) != length(local_pre2)) {{
  visr.message("The number of samples is different from the number of treatment conditions specified in factor 2", type="error")}}

if (input_advanced && input_GLM && identical(local_pre2[duplicated(local_pre2)], numeric(0))) {{
  visr.message("There are no replicates in the data for factor 2. To proceed: choose method -> blind; sharing mode -> fit-only.", type="warning")}}

# basic
if (!input_advanced) {{
  par(mfrow = c(1,1))
  input_columns <- c(input_group1, input_group2)
  group <- factor(c(rep(1,length(input_group1)), rep(2,length(input_group2))))

  #rownames(countsTable) <- countsTable$ID
  countsTable <- subset(input_table, select = input_columns)
  #countsTable <- countsTable / 50
  #head(countsTable)
  cds <- newCountDataSet(ceiling(countsTable), group)
  cds <- estimateSizeFactors(cds)
  
  #head( counts(cds) )
  sizeFactors(cds)
  cds <- estimateDispersions(cds, method=input_dispersionmethod, sharingMode=input_dispersionsharingMode, fitType = input_dispersionfitType)
  normalizedCounts <- t( t(counts(cds)) / sizeFactors(cds))
  #output_rle <- generateRLE(normalizedCounts)
  res <- nbinomTest( cds,  1, 2)
    
  output_baseMean <- res$baseMean
  output_baseMeanA <- res$baseMeanA
  output_baseMeanB <- res$baseMeanB
  output_foldChange <- res$foldChange
  output_log2FoldChange <- res$log2FoldChange
  output_pval <- res$pval
  output_padj <- res$padj
  output_cluster <- ifelse(output_padj < input_fdr, ifelse(output_log2FoldChange < 0, -1, 1), 0)
  output_cluster[is.na(output_cluster)] = 0

  par(mfrow = c(2,1))
  cluster_color <- ifelse(output_cluster == 0, "darkgrey", ifelse(output_cluster==1, "red", "blue"))
  plot(res[,input_plotx], res[,input_ploty], xlab=input_plotx, ylab=input_ploty, log = input_plotlog, col=cluster_color, cex=.3, pch=20)

  counts<-prop.table(table((output_cluster)))*length(output_cluster)
  bplt<-barplot(counts, main="DESeq Summary")
  text(y= counts, x= bplt, labels=as.character(counts), xpd=TRUE, adj=c(0.5,0))
}}

{{
  if (input_advanced && !input_GLM) {
    par(mfrow = c(2,2))

    #### design ####
    countTable <- input_table[,input_columns]

    input_factor1 = as.character(eval(parse(text = paste("c(",input_prefactor1,")"))))
    design = factor(input_factor1)
    cds = newCountDataSet(ceiling(countTable), design)

    ##### normalization (estimate the effective library size)
    cds = estimateSizeFactors(cds)
    sfcds <- sizeFactors(cds)
    # head(counts(cds, normalized=TRUE))

    ##### variance estimation
    cds = estimateDispersions(cds,method = input_dispersionmethod, sharingMode = input_dispersionsharingMode, fitType = input_dispersionfitType)
    # str(fitInfo(cds))
    plotDispEsts(cds) # dispersion plot: plot the per-gene estimates against the mean normalized counts per gene and overlay the fitted curve though estimates
    # head(fData(cds)) # the dispersion values used by the subsequent testing
    vsd = varianceStabilizingTransformation(cds)

    #### standard
    normalizedCounts <- t( t(counts(cds)) / sizeFactors(cds))
    output_rle <- generateRLE(normalizedCounts)
    res = nbinomTest(cds, "1","2")
    # head(res)
    # plotMA(res,col = ifelse(res$padj>=input_fdr, "gray32", "red3")) # plot the log2 fold changes against the mean normalised counts, colouring in red those genes that are significant at 10% FDR
    hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="") # Histogram of p-values from the call to nbinomTest.

    resSig = res[ res$padj < input_fdr, ] # filter for significant genes, according to some chosen threshold for the false dicovery rate (FDR),
    mostsig <- head( resSig[ order(resSig$pval), ] ) # the most significantly differentially expressed genes
    mostdown <- head( resSig[ order( resSig$foldChange, -resSig$baseMean ), ] ) # the most strongly down-regulated of the significant genes
    mostup <- head( resSig[ order( -resSig$foldChange, -resSig$baseMean ), ] ) # the most strongly up-regulated ones

    output_baseMean <- res$baseMean
    output_baseMeanA <- res$baseMeanA
    output_baseMeanB <- res$baseMeanB
    output_foldChange <- res$foldChange
    output_log2FoldChange <- res$log2FoldChange
    output_pval <- res$pval
    output_padj <- res$padj
    output_cluster <- ifelse(output_padj < input_fdr, ifelse(output_log2FoldChange < 0, -1, 1), 0)
    output_cluster[is.na(output_cluster)] = 0

    # draw fold change
    cluster_color <- ifelse(output_cluster == 0, "darkgrey", ifelse(output_cluster==1, "red", "blue"))
    plot(res[,input_plotx], res[,input_ploty], xlab=input_plotx, ylab=input_ploty, log = input_plotlog, col=cluster_color, cex=.3, pch=20)

    counts<-prop.table(table((output_cluster)))*length(output_cluster)
    bplt<-barplot(counts, main="DESeq Summary")
    text(y= counts, x= bplt, labels=as.character(counts), xpd=TRUE, adj=c(0.5,0))

    #### explore the data
    if(input_plot == "heatmap for the 30 most highly expressed genes") {
      select = order(rowMeans(counts(cds)), decreasing=TRUE)[1:30]
      hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
      heatmap.2(exprs(vsd)[select,], col = hmcol, trace="none", margin=c(10, 6)) # heatmap for the 30 most highly expressed genes
    }
    if(input_plot == "Heatmap of the sample-to-sample distances") {
      dists = dist(t(exprs(vsd))) #A heatmap of this distance matrix gives us an overview over similarities and dissimilarities between samples
      mat = as.matrix(dists)
      rownames(mat) = colnames(mat) = with(pData(cds), paste(condition,sep=" : "))
      heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
    }
    if(input_plot == "Principal component plot of the samples") {
      #principal component plot of the samples (PCA plot)
      print(plotPCA(vsd, intgroup=c("condition")))
    }
  } #if (input_advanced && !input_GLM)
}}

{{
  if (input_advanced && input_GLM)
  {
    par(mfrow = c(2,2))
    ################## glm ################

    #### design
    countTable <- input_table[,input_columns]

    #       if(input_singlefactor == TRUE) {
    #         input_factor1 = as.character(eval(parse(text = paste("c(",input_prefactor1,")"))))
    #         design = data.frame(
    #           row.names = colnames(input_table[,input_columns]),
    #           condition = input_factor1)
    #       }else{
    input_factor1 = as.character(eval(parse(text = paste("c(",input_prefactor1,")"))))
    input_factor2 = as.character(eval(parse(text = paste("c(",input_prefactor2,")"))))

    design = data.frame(
      row.names = colnames(input_table[,input_columns]),
      condition = input_factor1,
      libType = input_factor2)
    #       }

    cdsFull = newCountDataSet(ceiling(countTable), design ) # creating a count data set with multiple factors

    #### normalization (estimate the effective library size)
    cdsFull = estimateSizeFactors( cdsFull )
    sfcds <- sizeFactors(cdsFull)

    ##### variance estimation
    cdsFull = estimateDispersions(cdsFull, method = input_dispersionmethod, sharingMode = input_dispersionsharingMode, fitType = input_dispersionfitType)
    plotDispEsts(cdsFull)
    vsdFull = varianceStabilizingTransformation(cdsFull)
    normalizedCounts <- t( t(counts(cdsFull)) / sizeFactors(cdsFull))
    #output_rle <- generateRLE(normalizedCounts)
    #### model
    # model: compare the full model and reduced model in order to infer whether the additional specification of the treatment improves the fit and hence,
    # whether the treatment has significant effect
    #       if(input_singlefactor == TRUE) {
    #         fit1 = fitNbinomGLMs( cdsFull, count ~ condition )
    #         fit0 = fitNbinomGLMs( cdsFull, count ~ 1 )
    #       }else{
    fit1 = fitNbinomGLMs( cdsFull, count ~ libType + condition )
    fit0 = fitNbinomGLMs( cdsFull, count ~ libType )
    #       }
    #str(fit1)

    #### perform the test:
    pvalsGLM = nbinomGLMTest( fit1, fit0 )
    padjGLM = p.adjust( pvalsGLM, method="BH" )

    if (input_ruvseq) {
      # To estimate the factors of unwanted variation, we need a set of negative control genes,
      # e.g., least significantly DE genes based on a first-pass DE analysis performed prior to RUVg normalization

      set <- newSeqExpressionSet(as.matrix(ceiling(countTable)),
                                 phenoData = data.frame(input_factor1, row.names=colnames(countTable)))
      #set
      res  = cbind(fit1, pval = pvalsGLM, padj = padjGLM)
      res2  = res[ order(res$padj), ]
      # we consider all but the top 5000 genes as ranked by p-values.
      ruvempirical <- rownames(set)[which(!(rownames(set) %in% rownames(res2)[1:input_ruvn]))]
      set2 <- RUVg(set, ruvempirical, k=1)
      pData(set2)

      #design <- model.matrix(~input_group + W_1, data=pData(set2))
      #         if(input_singlefactor == TRUE) {
      #           design = data.frame(
      #             row.names = colnames(input_table[,input_columns]),
      #             condition = input_factor1,
      #             variation = pData(set2)$W_1)
      #         }else{
      design = data.frame(
        row.names = colnames(input_table[,input_columns]),
        condition = input_factor1,
        libType = input_factor2,
        variation = pData(set2)$W_1)
      #         }

      cdsFull = newCountDataSet(ceiling(countTable), design ) # creating a count data set with multiple factors

      #### normalization (estimate the effective library size)
      cdsFull = estimateSizeFactors( cdsFull )
      sfcds <- sizeFactors(cdsFull)

      ##### variance estimation
      cdsFull = estimateDispersions(cdsFull, method = input_dispersionmethod, sharingMode = input_dispersionsharingMode, fitType = input_dispersionfitType)
      normalizedCounts <- t( t(counts(cdsFull)) / sizeFactors(cdsFull))
      #output_rle <- generateRLE(normalizedCounts)
      
      #### model
      # model: compare the full model and reduced model in order to infer whether the additional specification of the treatment improves the fit and hence,
      # whether the treatment has significant effect
      #         if(input_singlefactor == TRUE) {
      #           fit1 = fitNbinomGLMs( cdsFull, count ~ condition + variation)
      #           fit0 = fitNbinomGLMs( cdsFull, count ~ variation)
      #         }else{
      fit1 = fitNbinomGLMs( cdsFull, count ~ libType + condition + variation)
      fit0 = fitNbinomGLMs( cdsFull, count ~ libType +variation)
      #         }
      #### perform the test:
      pvalsGLM = nbinomGLMTest( fit1, fit0 )
      padjGLM = p.adjust( pvalsGLM, method="BH" )
    } #if(input_ruvseq)

    #extract the significant genes from the vector padjGLM
    res  = cbind(fit1, pval = pvalsGLM, padj = padjGLM)
    resSig  = res[which(res$padj  < input_fdr),] # filter for significant genes, according to some chosen threshold for the false dicovery rate (FDR),
    hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="") # Histogram of p-values from the call to the test.

    mostsig<-head( resSig[ order(resSig$padj), ] ) # the most significantly differentially expressed genes:
    mostdown<- head( resSig[ order(resSig[,3], decreasing = FALSE), ] )
    mostup<- head( resSig[ order(resSig[,3], decreasing = TRUE), ] )

    # head(fit1) # the corresponding fold changes

    output_pval <- res$pval
    output_padj <- res$padj
    output_cluster <- ifelse(output_padj < input_fdr, ifelse(res[,3] < 0, -1, 1), 0)
    output_cluster[is.na(output_cluster)] = 0

    # draw fold change
    cluster_color <- ifelse(output_cluster == 0, "darkgrey", ifelse(output_cluster==1, "red", "blue"))
    #plot(res[,1], res[,3], log = "xy", xlab="base", ylab="treatment", col=cluster_color, cex=.3, pch=20)

    counts<-prop.table(table((output_cluster)))*length(output_cluster)
    bplt<-barplot(counts, main="DESeq Summary")
    text(y= counts, x= bplt, labels=as.character(counts), xpd=TRUE, adj=c(0.5,0))

    #### explore the data
    if(input_plot == "heatmap for the 30 most highly expressed genes") {
      select = order(rowMeans(counts(cdsFull)), decreasing=TRUE)[1:30]
      hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
      heatmap.2(exprs(vsdFull)[select,], col = hmcol, trace="none", margin=c(10, 6)) # heatmap for the 30 most highly expressed genes
    }

    if(input_plot == "Heatmap of the sample-to-sample distances") {
      dists = dist( t( exprs(vsdFull) ) ) #A heatmap of this distance matrix gives us an overview over similarities and dissimilarities between samples
      mat = as.matrix( dists )
      rownames(mat) = colnames(mat) = with(pData(cdsFull), paste(condition, libType, sep=" : "))
      heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
    }

    if(input_plot == "Principal component plot of the samples") {
      #principal component plot of the samples (PCA plot)
      print(plotPCA(vsdFull, intgroup=c("condition", "libType")))
    }
  } # if (input_advanced && input_GLM)
}}

{{
  if (input_advanced) {
    print(capture.output(cat("\n1.design\n")),collapse='\n')
    print(capture.output(print(design)),collapse='\n')

    print(capture.output(cat("\n2.normalization (estimate the effective library size)\n")),collapse='\n')
    print(capture.output(cat("\nsize factors\n")),collapse='\n')
    print(capture.output(print(sfcds)),collapse='\n')

    print(capture.output(cat("\n3.variance estimation\n")),collapse='\n')

    print(capture.output(cat("\n4.explore the data\n")),collapse='\n')

    print(capture.output(cat("\n5.perform the test\n")),collapse='\n')

    print(capture.output(cat("\nthe most significantly differentially expressed genes:\n")),collapse='\n')
    print(capture.output(print(mostsig)),collapse='\n')
    print(capture.output(cat("\nthe most strongly down-regulated of the significant genes:\n")),collapse='\n')
    print(capture.output(print(mostdown)),collapse='\n')
    print(capture.output(cat("\nthe most strongly up-regulated ones:\n")),collapse='\n')
    print(capture.output(print(mostup)),collapse='\n')

    print(capture.output(print(counts)),collapse='\n')
  } #if (input_advanced)
}}
