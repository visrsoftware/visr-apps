biocPackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {source("http://bioconductor.org/biocLite.R"); biocLite(pkg);require(pkg, character.only = TRUE)}}
biocPackage("edgeR")

#input_group1 <- c("TT2_count")
#input_group2 <- c("SETDB1KO_count")
#input_adjust_method <- "BH"
input_pvalue <- 0.01
input_cpm_cutoff <- 1
#input_isGLM <- TRUE

visr.applyParameters()

if (length(input_group1) == 1) input_group1 <- c(input_group1, input_group1)

if (length(input_group2) == 1) input_group2 <- c(input_group2, input_group2)

input_columns <- c(input_group1, input_group2)
group <- factor(c(rep(1,length(input_group1)), rep(2,length(input_group2))))

x <- subset(visr.input, select = input_columns)
#x <- x / input_read_length

y <- DGEList(counts=x, group=group)

#filtering
keep <- rowSums(cpm(y) > input_cpm_cutoff) >= 1
y <- y[keep,]
dim(y)

{{
  if (input_isGLM)
  {
    design <- model.matrix(~group)
    y <- estimateGLMCommonDisp(y, design)#, method="deviance", robust=TRUE, subset=NULL)
    y <- estimateGLMTrendedDisp(y,design)
    y <- estimateGLMTagwiseDisp(y,design)
    fit <- glmFit(y,design)
    lrt <- glmLRT(fit,coef=2)
    de <- decideTestsDGE(lrt, p.value=input_pvalue, adjust.method=input_adjust_method)
    topTags(lrt)
  } else {
    y <- calcNormFactors(y)
    y <- estimateCommonDisp(y)
    y <- estimateTagwiseDisp(y)
    et <- exactTest(y)
    de <- decideTestsDGE(et, p.value=input_pvalue, adjust.method=input_adjust_method)
  }
}}
summary(de)

output_decide<-rep(0, length(keep))
output_decide[keep]<-de

#hist(output_decide, col="gray", labels = TRUE)
counts<-prop.table(table((output_decide)))*length(output_decide)
#bplt<-barplot(counts, main="EdgeR Summary")
#text(y= counts, x= bplt, labels=as.character(counts), xpd=TRUE, adj=c(0.5,0))
#cluster_color <- ifelse(output_decide == 0, "darkgrey", ifelse(output_decide==1, "red", "blue"))
#plot(res[[input_plot_x]], res[[input_plot_y]], xlab=input_plot_x, ylab=input_plot_y, log = "xy", col=cluster_color, cex=.3, pch=20)
detags <- rownames(y)[as.logical(de)]
plotSmear(et, de.tags=detags)