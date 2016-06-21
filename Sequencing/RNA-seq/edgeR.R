source("visrutils.R")
visr.biocLite("edgeR")
visr.biocLite("RUVSeq")

visr.applyParameters()

if (!input_advanced && length(input_group1) == 1)
  visr.message("There are no replicates in group 1.", type="warning")

if (!input_advanced && length(input_group2) == 1)
  visr.message("There are no replicates in group 2.", type="warning")

if (input_advanced)
  input_group_c <- eval(parse(text = paste("c(",input_group,")")))

if (input_advanced  && length(input_columns) != length(input_group_c))
  visr.message("The number of samples is different from the number of treatment conditions specified in group", type="error")
if (input_advanced  && identical(input_group_c[duplicated(input_group_c)], numeric(0)))
  visr.message("There are no replicates in the data.", type="warning")

if (input_advanced) {
  if (input_experimentaldesign == "full interaction" || input_experimentaldesign == "blocking") {
    if (length(input_columns) != length(eval(parse(text = paste("c(",input_prefactor1,")")))))
      visr.message("The number of samples is different from the number of treatment conditions specified in factor 2", type="error")
    if (length(input_columns) != length(eval(parse(text = paste("c(",input_prefactor2,")")))))
      visr.message("The number of samples is different from the number of treatment conditions specified in factor 2", type="error")
  }
}


#' RLE normalization
generateRLE <- function (object) {
  nARR = dim(object)[2]
  nGEN = dim(object)[1]
  y = apply(object, 1, median)
  rle = matrix(nrow = nGEN, ncol = nARR)
  for (i in 1:nARR) {
    x = object[, i]
    rle[, i] = log(x/y)
  }
  colnames(rle) = colnames(object)
  return(rle)
}

if (!input_advanced) {
  par(mfrow = c(1,1))

  if (length(input_group1) == 1) {
    input_group1 <- c(input_group1, input_group1)
  }

  if (length(input_group2) == 1) {
    input_group2 <- c(input_group2, input_group2)
  }

  input_columns <- c(input_group1, input_group2)
  group <- factor(c(rep(1,length(input_group1)), rep(2,length(input_group2))))

  x <- subset(input_table, select = input_columns)
  #x <- x / input_read_length

  y <- DGEList(counts=x, group=group)

  #filtering
  keep <- rowSums(cpm(y) >= input_cpm_cutoff) >= 1
  y <- y[keep,]
  dim(y)


  if(input_method == "exact test") {
    y <- calcNormFactors(y)
    y <- estimateCommonDisp(y)
    y <- estimateTagwiseDisp(y)
    #output_rle <- generateRLE(y$pseudo.counts)
    output_normalized <- y$pseudo.counts
    et <- exactTest(y)
    de <- decideTestsDGE(et, p.value=input_pvalue, adjust.method=input_adjustmethod)
  } else { #input_method == "glm"
    design <- model.matrix(~group)
    y <- estimateGLMCommonDisp(y, design)#, method="deviance", robust=TRUE, subset=NULL)
    y <- estimateGLMTrendedDisp(y,design)
    y <- estimateGLMTagwiseDisp(y,design)
    #output_rle <- generateRLE(cpm(y))
    output_normalized <- cpm(y)
    fit <- glmFit(y,design)
    lrt <- glmLRT(fit,coef=2)
    de <- decideTestsDGE(lrt, p.value=input_pvalue, adjust.method=input_adjustmethod)
    topTags(lrt)
  }

  summary(de)

  output_decide<-rep(0, length(keep))
  output_decide[keep]<-de

  #hist(output_decide, col="gray", labels = TRUE)
  counts<-prop.table(table((output_decide)))*length(output_decide)
  bplt<-barplot(counts, main="EdgeR Summary")
  text(y= counts, x= bplt, labels=as.character(counts), xpd=TRUE, adj=c(0.5,0))
}

if (input_advanced) {
  input_verbose <- TRUE
  par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
  design <- NULL

  ####1 read data ####
  print(capture.output(cat("\n1.read data\n")),collapse='\n')


  if (anyDuplicated(eval(parse(text = paste("c(",input_group,")")))) == 0) {
    input_group = c(input_group, input_group)
    input_columns = c(input_columns, input_columns)
  }


  input_group = as.character(eval(parse(text = paste("c(",input_group,")"))))
  y <- DGEList(counts=input_table[,input_columns], group=input_group)
  # y$samples
  # head(y$counts)
  # summary(y$counts)
  print(capture.output(cat("\ndimension of the data\n")),collapse='\n')
  print(capture.output(dim(y)),collapse='\n')
  # dim(y)

  ####2 filtering ####
  print(capture.output(cat("\n2.filter low expression tags\n")),collapse='\n')

  keep <- rowSums(cpm(y) >= input_cpm_cutoff) >= min(table(input_group)) # filter out very lowly expressed tags, keeping genes that are expressed at a reasonable level in at least one treatment condition.
  y <- y[keep,,keep.lib.sizes=FALSE] # recompute the library sizes
  print(capture.output(cat("\ndimension of the data after filtering\n")),collapse='\n')
  print(capture.output(print(dim(y))),collapse='\n')
  # dim(y)

  ####3 normalization ####
  print(capture.output(cat("\n3.normalize data\n")),collapse='\n')

  y <- calcNormFactors(y, method= input_normalizationmethod)
  # d$samples

  ####4 data exploration ####
  print(capture.output(cat("\n4.explore the data\n")),collapse='\n')

  plotMDS(y, method=input_mdsmethod) # shows distances between samples
}


##################### exact test #####################
if (input_advanced && input_method == "exact test") {

  ####5 estimate the dispersions ####

  y <- estimateCommonDisp(y, verbose=input_verbose)
  y <- estimateTagwiseDisp(y, trend=input_trend) # estimate the dispersion parameter for each tag, a measure of the degree of inter-library variation for that tag
  #output_rle <- generateRLE(y$pseudo.counts)
  output_normalized <- y$pseudo.counts
  bcv<- sqrt(y$common.disp) ## BCV

  plotBCV(y) # plots the tagwise dispersions against log2-CPM

  ####6 tests ####
  #print(capture.output(cat("\nperform tests\n")),collapse='\n')

  et <- exactTest(y,pair = eval(parse(text = paste("c(",input_pair,")"))))

  ptoptags <- topTags(et, n=input_n) # The test results for the n most significant tags

  detags <- rownames(topTags(et, n=input_n)$table)

  pcpm <- cpm(y)[detags, order(y$samples$group)]  # counts per million for the tags that edgeR has identified as the most differentially expressed

  psummary <- summary(de <- decideTestsDGE(et, p=input_pvalue, adjust.method=input_adjustmethod)) # The total number of differentially expressed genes at FDR< 0.05

  detags <- rownames(y)[as.logical(de)]
  plotSmear(et, de.tags=detags) # A smearplot displays the tagwise log-fold changes against log-cpm with the DE genes highlighted
  abline(h = c(-input_foldchanges, input_foldchanges), col = "blue") # 4-fold changes
}

##################### glm ####################
if (input_advanced && input_method != "exact test") { #input_method == "glm"
  ####5 design matrix ####
  #print(capture.output(cat("\ndesign matrix\n")),collapse='\n')

  input_factor1 = as.character(eval(parse(text = paste("c(",input_prefactor1,")"))))
  input_factor2 = as.character(eval(parse(text = paste("c(",input_prefactor2,")"))))

  if(input_experimentaldesign == "simple") {design <- model.matrix(~input_group, data=y$samples)}
  if(input_experimentaldesign == "full interaction") {design <- model.matrix(~input_factor1*input_factor2, data=y$samples)}
  if(input_experimentaldesign == "blocking") {design <- model.matrix(~input_factor1+input_factor2, data=y$samples)}
  design


  ####6 estimate the dispersions ####
  #print(capture.output(cat("\nestimate the dispersions\n")),collapse='\n')

  y <- estimateGLMCommonDisp(y, design, verbose=input_verbose)
  y <- estimateGLMTrendedDisp(y, design) # a possible trend with averge count size
  y <- estimateGLMTagwiseDisp(y, design)
  #output_rle <- generateRLE(cpm(y))
  output_normalized <- cpm(y)

  bcv<- sqrt(y$common.disp)
  plotBCV(y)

  ####7 tests ####
  #paste(capture.output(cat("\nperform tests\n")),collapse='\n')

  #colnames(design)
  fit <- glmFit(y, design)


  if(input_coef == "" & input_contrast =="") {
    lrt <- glmLRT(fit)
  }else{
    if(input_coef == "" & input_contrast !=""){
      input_contrast2 = eval(parse(text = paste("c(",input_contrast,")")))
      lrt <- glmLRT(fit,contrast=input_contrast2)
    }else{
      lrt <- glmLRT(fit,coef= input_coef)
    }
  }


  if(input_ruvseq == TRUE) {
    # To estimate the factors of unwanted variation, we need a set of negative control genes,
    # e.g., least significantly DE genes based on a first-pass DE analysis performed prior to RUVg normalization

    set <- newSeqExpressionSet(as.matrix(y$counts),
                               phenoData = data.frame(input_group, row.names=colnames(y$counts)))
    #set
    ruvtop <- topTags(lrt, n=nrow(set))$table
    # we consider all but the top 5000 genes as ranked by edgeR p-values.
    ruvempirical <- rownames(set)[which(!(rownames(set) %in% rownames(ruvtop)[1:input_ruvn]))]

    set2 <- RUVg(set, ruvempirical, k=1)
    pData(set2)

    #plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[input_group])
    #plotPCA(set2, col=colors[input_group], cex=1.2)

    if(input_experimentaldesign == "simple") {design <- model.matrix(~input_group + W_1, data=pData(set2))}
    if(input_experimentaldesign == "full interaction") {design <- model.matrix(~input_factor1*input_factor2 + W_1, data=pData(set2))}
    if(input_experimentaldesign == "blocking") {design <- model.matrix(~input_factor1+input_factor2 + W_1, data=pData(set2))}
    #design <- model.matrix(~input_group + W_1, data=pData(set2))
    y <- DGEList(counts=counts(set), group=input_group)
    y <- calcNormFactors(y, method="upperquartile")
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    #output_rle <- generateRLE(cpm(y))
    output_normalized <- cpm(y)

    fit <- glmFit(y, design)

    if(input_coef == "" & input_contrast =="") {
      lrt <- glmLRT(fit)
    }else{
      if(input_coef == "" & input_contrast !=""){
        input_contrast2 = eval(parse(text = paste("c(",input_contrast,")")))
        lrt <- glmLRT(fit,contrast=input_contrast2)
      }else{
        lrt <- glmLRT(fit,coef= input_coef)
      }
    }
  }

  # likelihood ratio tests for tumour vs normal tissue differences
  # Note that glmLFT has conducted a test for the last coefficient in the linear model, which we can see is the tumor vs normal tissue effect
  # conduct likelihood ratio tests for the pathogen effect and show the top genes
  # the test is for the last coefficient in the design matrix, which in this case is the treatment effect

  ptoptags <- topTags(lrt, n=input_n) # show the top genes

  # The top DE tags have tiny p-values and FDR values, as well as large fold changes.
  top <- rownames(topTags(lrt, n=input_n))
  pcpm <- cpm(y)[top,] # counts-per-million in individual samples for the top genes

  psummary <- summary(de <- decideTestsDGE(lrt, p=input_pvalue, adjust.method=input_adjustmethod)) # The total number of differentially expressed genes at 5% FDR

  isDE <- as.logical(de) # pick out DE genes
  detags <- rownames(y)[isDE]
  plotSmear(lrt, de.tags=detags) # plot all the logFCs against average count size, highlighting the DE genes
  abline(h=c(-input_foldchanges, input_foldchanges), col="blue")
}

if (input_advanced) {
  output_decide<-rep(0, length(keep))
  output_decide[keep]<-de
  counts<-prop.table(table((output_decide)))*length(output_decide)
  bplt<-barplot(counts, main="EdgeR Summary")
  text(y= counts, x= bplt, labels=as.character(counts), xpd=TRUE, adj=c(0.5,0))


  print(capture.output(cat("\n5.estimate the dispersions\n")),collapse='\n')
  print(capture.output(cat("\nbiological coefficient of variation\n")),collapse='\n')
  print(capture.output(print(bcv)),collapse='\n')

  print(capture.output(cat("\n6.design\n")),collapse='\n')
  print(capture.output(print(design)),collapse='\n')

  print(capture.output(cat("\n7.perform tests\n")),collapse='\n')
  print(capture.output(cat("\ntest results for the most significant tags\n")),collapse='\n')
  print(capture.output(print(ptoptags)),collapse='\n')

  print(capture.output(cat("\ncounts per million for the tags that edgeR has identified as the most differentially expressed\n")),collapse='\n')
  print(capture.output(print(pcpm)),collapse='\n')

  #print(capture.output(print(psummary)),collapse='\n')
  print(capture.output(print(counts)),collapse='\n')
}
