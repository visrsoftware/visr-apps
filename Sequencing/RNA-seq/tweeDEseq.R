biocPackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {source("http://bioconductor.org/biocLite.R"); biocLite(pkg);require(pkg, character.only = TRUE)}}
biocPackage("tweeDEseq")
biocPackage("edgeR")

# input_table <- read.csv("C:/Users/hexu/Desktop/VisRSeq/VisRSeq/data/counts_tab.txt",sep = "\t")
# input_table <- input_table[1:100,]
# input_columns <- c("CT.PA.1","CT.PA.2","KD.PA.3","KD.PA.4") #,"CT.SI.5","KD.SI.6","CT.SI.7")
# # input_sample1 <- c("CT.PA.1","CT.PA.2")
# # input_sample2 <- c("KD.PA.3","KD.PA.4")
# #input_precon <- "1,1,2,2"
# input_pregroup = "1,2,1,2"

# input_method = c("TMM", "cqn")[1]
# input_commondisp = FALSE
# input_priordf = 8
# 
# input_cpmcutoff = 0.5
# input_nsamplescutoff = 2
# input_meancpmcutoff = 0
# 
# input_log2fc = log2(1.5)
# input_fdr = 0.05


visr.applyParameters()

par(mfrow=c(2,2), mar=c(4, 5, 3, 2))

input_group = as.factor(as.character(eval(parse(text = paste("c(",input_pregroup,")")))))
if(length(input_columns) != length(input_group)) error_message<-"the number of samples is different from the number of treatment conditions specified in group"
test <- eval(parse(text = paste("c(",input_pregroup,")")))
if(identical(test[duplicated(test)], numeric(0))) error_message <- "There are no replicates in the data."

data <- normalizeCounts(input_table[,input_columns], input_group, method = input_method, common.disp = input_commondisp, prior.df = input_priordf)
print(capture.output(cat("\nnormalize the initial table to remove technical biases\n")),collapse='\n') 
print(capture.output(cat("\ndataset dimension\n")),collapse='\n') 
print(capture.output(dim(data)),collapse='\n')

data <- filterCounts(data, cpm.cutoff = input_cpmcutoff, n.samples.cutoff = input_nsamplescutoff, mean.cpm.cutoff = input_meancpmcutoff)
print(capture.output(cat("\nfilter out genes with very low espression\n")),collapse='\n') 
print(capture.output(cat("\ndataset dimension\n")),collapse='\n') 
print(capture.output(dim(data)),collapse='\n')

# show us by default the top 6 genes ranked by P-value including information on
# the magnitude of the fold-change in log2 scale (log2fc), overall mean expression in counts
# (overallMean), mean expression in counts for each sample group, raw P-value (pval) and the
# Benjamini-Hochberg (FDR) adjusted P-value (pval.adjust)

resPT <- tweeDE(data, group = input_group)

print(resPT)

output_overallMean <- resPT$overallMean
output_mean1 <- resPT$X1
output_mean2 <- resPT$X2
output_log2FoldChange <- resPT$log2fc
output_stat <- resPT$stat
output_pvalue <- resPT$pval
output_padj <- resPT$pval.adjust
output_cluster <- ifelse(output_padj < input_fdr, ifelse(output_log2FoldChange < 0, -1, 1), 0)
output_cluster[is.na(output_cluster)] = 0

# Histogram of p-values 
hist(output_padj, breaks=100, col="skyblue", border="slateblue", main="Histogram of p-values") 
## barplot
counts<-prop.table(table((output_cluster)))*length(output_cluster)
bplt<-barplot(counts, main="tweeDEseq Summary")
text(y= counts, x= bplt, labels=as.character(counts), xpd=TRUE, adj=c(0.5,0))

print(capture.output(cat("\nsummary\n")),collapse='\n') 
print(capture.output(print(counts)),collapse='\n')

# allows us to call differentially expressed a subset of gene meeting
# cutoffs on the minimum magnitude of the fold-change and maximum FDR and store that subset
# in a data.frame object

deGenes <- print(resPT, n=Inf, log2f.cutoff=input_log2fc, pval.adjust.cutoff=input_fdr, print=FALSE)

dim(deGenes)

 hl <- list(list(genes=rownames(deGenes), pch=1, col="red", lwd=2, cex=1.5), list(genes=rownames(resPT), pch=21, col="blue", bg="blue", cex=0.7))
#   list(genes=XiEgenes, pch=21, col="skyblue", bg="skyblue", cex=0.7),
         

MAplot(resPT, cex=0.7, highlight=hl, log2fc.cutoff=input_log2fc, main="MA-plot")
Vplot(resPT, cex=0.7, highlight=hl,log2fc.cutoff=input_log2fc,pval.adjust.cutoff=input_fdr, main="Volcano plot")

# mod <- glmPT(data ~ input_group)
# mod
# summary(mod)
# anova(mod)
# 
# 
# resPTglm <- tweeDEglm( ~ input_group, counts = data)
# head(resPTglm[sort(resPTglm$pval.adjust, index.return = TRUE)$ix,])
# 
# # performs conditional quantile normalization in order to remove possibly existing bias arising from differences in GC content or gene lengths.
# tweeDEglm( ~ input_group, counts = data,offset = cqn.subset$offset)
