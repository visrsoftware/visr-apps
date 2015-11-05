biocPackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {source("http://bioconductor.org/biocLite.R"); biocLite(pkg);require(pkg, character.only = TRUE)}}
biocPackage("sSeq")

# input_table <- read.csv("C:/Users/hexu/Desktop/old versions/VisRSeq.0.72.7/VisRSeq.0.72.7/data/counts_tab.txt",sep = "\t")
# input_table <- read.table("C:/Users/hexu/Desktop/old versions/VisRSeq.0.72.7/VisRSeq.0.72.7/data/test.txt",sep="\t",header=TRUE)
# #input_columns <- c("CT.PA.1","CT.PA.2","KD.PA.3","KD.PA.4","CT.SI.5","KD.SI.6","CT.SI.7")
# input_sample1 <- c("CT.PA.1","CT.PA.2")
# input_sample2 <- c("KD.PA.3","KD.PA.4")
# #input_precon <- "1,1,2,2"
# input_precol = "1,2,1,2"

# input_pairedDesign = FALSE
# input_pairedDesigndispMethod = c("per-pair","pooled")[2]
# input_useFisher = FALSE
# 
# input_cutoff = 0.05


visr.applyParameters()


if (length(input_sample1) == 1) error_message <- "Warning: there are no replicates in group 1."  
if (length(input_sample2) == 1) error_message <- "Warning: there are no replicates in group 2."

par(mfrow = c(2,3))

countsTable1 = input_table[,input_sample1]
countsTable2 = input_table[,input_sample2]
countsTable = cbind(countsTable1, countsTable2)
#input_con = as.character(eval(parse(text = paste("c(",input_precon,")"))))
input_con = c(rep("A",length(input_sample1)),rep("B",length(input_sample2)))
#exact test to get differential expressed genes
{{
  if( input_precol == "") {
  res1 = nbTestSH(countsTable, 
                conds = input_con, "A", "B",
                useFisher = input_useFisher)
}else{
  input_coLevels=data.frame(subjects=(eval(parse(text = paste("c(",input_precol,")")))))
  res1 = nbTestSH(countsTable, 
                  conds = input_con, "A", "B",
                  coLevels = input_coLevels,
                  pairedDesign = input_pairedDesign, 
                  pairedDesign.dispMethod = input_pairedDesigndispMethod,
                  useFisher = input_useFisher)
}
}}

print(capture.output(cat("\nexact test to get differential expressed genes\n")),collapse='\n') 
print(capture.output(head(res1)),collapse='\n')

# head(res1)
output_mean <- res1$Mean  # The row per-gene averages over the values in countsTable
output_log2FoldChange <- res1$rawLog2FoldChange
output_pvalue <- res1$pval
output_cluster <- ifelse(output_pvalue < 0.1, ifelse(output_log2FoldChange < 0, -1, 1), 0)
output_cluster[is.na(output_cluster)] = 0

counts<-prop.table(table((output_cluster)))*length(output_cluster)
print(capture.output(cat("\nsummary\n")),collapse='\n') 
print(capture.output(print(counts)),collapse='\n')

# a plot of ASD values when varying the shrinkage targets
{{
  if( input_precol == "") {
    disp1 = nbTestSH(countsTable, conds = input_con, "A", "B",SHonly=TRUE, plotASD=TRUE)
  }else{
    disp1 = nbTestSH(countsTable, 
                     conds = input_con, "A", "B",
                     coLevels = input_coLevels,
                     pairedDesign = input_pairedDesign, 
                     pairedDesign.dispMethod = input_pairedDesigndispMethod,
                     SHonly=TRUE, 
                     plotASD=TRUE)
  }
}}
# a scatter plot visualizing the relationship between the dispersion estimates and the mean estimates
plotDispersion(disp1, legPos="bottomright")

# variance plot
rV = rowVars(countsTable);
mu = rowMeans(countsTable); 
SH.var =  (disp1$SH * mu^2 + mu)
smoothScatter(log(rV)~log(mu), main="Variance Plot", ylab='log(variance)', xlab='log(mean)', col=blues9[5], cex.axis=1.8)
points( log(SH.var)~log(mu), col="black", pch=16)
leg1 =  expression(paste("log(", hat("V")[g]^"sSeq", ")", sep=''));
leg2 =  expression(paste("log(", hat("V")[g]^"MM", ")", sep=''));
legend("bottomright", legend=c(leg1,leg2), col=c("black",blues9[5]), pch=c(16, 1), cex=0.8)


#obtain the p-values for the comparison AvsA.
conds2.Hammer = c("A","B");
res1.2 =  nbTestSH( countsTable[,1:2], conds2.Hammer, "A", "B");
#draw the ECDF plot;
dd = data.frame(AvsA=res1.2$pval, AvsB=res1$pval);
ecdf <- ecdfAUC(dd, col.line=c("green", "red"), main = "ECDF, Hammer", drawRef = TRUE, rm1=TRUE)
print(capture.output(cat("\ndraw the ECDF plot\n")),collapse='\n') 
print(capture.output(print(ecdf)),collapse='\n')

# MV plot helps to detect any dependent structures between the means and the differences in condition A and B. 
# In a MV plot, we expect to see that the dots are roughly distributed on the two sides of the zero horizontal line without any
# dependent patter between M and V. 
# A volcano plot is a scatter plot that visualizes the linear dependence 
# between the statistical changes (e.g. -log2(p-value)) and the biological changes (e.g. log2(fold change)). 
# We expect to see that the dots are linearly and evenly distributed on the two sides of the zero vertical line.
{{
  drawMA_vol.2 <- function (y, groups2, pv, cutoff = NULL, xlab1 = "(log2(A)+log2(B))/2", 
                            ylab1 = "log2(A)-log2(B)", tt1 = "MA plot", tt2 = "volcano plot", 
                            log2FoldChange = NULL, col1 = c("black", "red")) {
    AB = unique(groups2)
    A = AB[1]
    B = AB[2]
    if (length(pv) != dim(y)[[1]]) 
      stop("check if pv is matched to y for drawMAfunction.")
    if (ncol(y) != length(groups2) & length(unique(groups2)) != 
          2) 
      stop("check y and groups2 for drawMA function.")
    muA = rowMeans(as.matrix(y[, groups2 %in% A]))
    muA[muA == 0] = 1
    muB = rowMeans(as.matrix(y[, groups2 %in% B]))
    muB[muB == 0] = 1
    if (is.null(log2FoldChange)) 
      log2FoldChange = log2(muA/muB)
    baseMean = log2(muA * muB)/2
    # op <- par(mfrow = c(1, 2))
    if (is.null(xlab1) | is.null(ylab1)) {
      cond1 = unique(groups2)
      leg.A = cond1[1]
      leg.B = cond1[2]
      xlab1 = paste("(log2(", leg.A, ")+log2(", leg.B, "))/2", 
                    sep = "")
      ylab1 = paste("(log2(", leg.A, ")-log2(", leg.B, "))/2", 
                    sep = "")
    }
    plot(y = log2FoldChange, x = baseMean, pch = 16, cex = 0.3, 
         col = col1[1], xlab = xlab1, ylab = ylab1, main = tt1)
    if (is.null(cutoff)) {
      cutoff = min(pv[order(pv)][floor(nrow(y) * 0.05)], 0.05)
      print(paste("cutoff", cutoff))
    }
    sel = pv < cutoff
    points(x = baseMean[sel], y = log2FoldChange[sel], col = col1[2], 
           cex = 0.6, pch = 16)
    abline(h = 0, col = "blue")
    log2pv = -log2(pv)
    plot(y = log2pv, x = log2FoldChange, pch = 16, cex = 0.3, 
         col = col1[1], ylab = "-log2(pvalue)", xlab = "log2(fold change)", 
         main = tt2)
    points(y = log2pv[sel], x = log2FoldChange[sel], col = col1[2], 
           cex = 0.3, pch = 16)
    abline(v = 0, col = "blue")
  }
}}
drawMA_vol.2(countsTable, input_con, res1$pval, cutoff=input_cutoff);


