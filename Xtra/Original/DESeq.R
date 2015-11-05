biocPackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {source("http://bioconductor.org/biocLite.R"); biocLite(pkg);require(pkg, character.only = TRUE)}}
biocPackage("DESeq")

input_plot_drawFoldChange<-FALSE

visr.applyParameters()

input_columns <- c(input_group1, input_group2)
group <- factor(c(rep(1,length(input_group1)), rep(2,length(input_group2))))

#rownames(countsTable) <- countsTable$ID
countsTable <- subset(visr.input, select = input_columns)
#countsTable <- countsTable / 50

#head(countsTable)
cds <- newCountDataSet(ceiling(countsTable), group)
cds <- estimateSizeFactors(cds)
#head( counts(cds) )
sizeFactors(cds)
cds <- estimateDispersions(cds, method=input_method, sharingMode=input_sharingMode, fitType=input_fitType)
res <- nbinomTest( cds,  1, 2)

output_baseMean <- res$baseMean
output_baseMeanA <- res$baseMeanA
output_baseMeanB <- res$baseMeanB
output_foldChange <- res$foldChange
output_log2FoldChange <- res$log2FoldChange
output_pval <- res$pval
output_padj <- res$padj
output_cluster <- ifelse(output_padj < input_FDRThreshold, ifelse(output_log2FoldChange < 0, -1, 1), 0)
output_cluster[is.na(output_cluster)] = 0

{{
if (input_plot_drawFoldChange) {
  cluster_color <- ifelse(output_cluster == 0, "darkgrey", ifelse(output_cluster==1, "red", "blue"))
  plot(res[[input_plot_x]], res[[input_plot_y]], xlab=input_plot_x, ylab=input_plot_y, log = input_plot_log, col=cluster_color, cex=.3, pch=20)
} else {
  counts<-prop.table(table((output_cluster)))*length(output_cluster)
  bplt<-barplot(counts, main="DESeq Summary")
  text(y= counts, x= bplt, labels=as.character(counts), xpd=TRUE, adj=c(0.5,0))
}
}}

