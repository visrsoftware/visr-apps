usePackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {install.packages(pkg, repos = "http://cran.us.r-project.org");require(pkg, character.only = TRUE)}}
usePackage("data.table")
usePackage("devtools")
usePackage("graphics")
biocPackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {source("http://bioconductor.org/biocLite.R"); biocLite(pkg);require(pkg, character.only = TRUE)}}
biocPackage("GenomicRanges")
biocPackage("IRanges")

if(!require(methylKit)){install_github("al2na/methylKit",build_vignettes=FALSE);require(methylKit)}

visr.applyParameters()

### 1.read methylation data

{{
  if(input_ann == FALSE & input_cpgann == FALSE) {
    par(mfrow=c(2,3))
  }else{
    if(input_ann == TRUE & input_cpgann == TRUE) {
      par(mfrow=c(2,4))
    }else{
      par(mfrow=c(2,4))
    }
  }
}}


file.list <- list(input_file1,input_file2,input_file3,input_file4)#,input_file5,input_file6,input_file7,input_file8)
ind<- c()
for(i in 1:length(file.list)){if (file.list[i] == "") ind <- append(ind,i)}
if(!is.null(ind)) file.list <- file.list[-ind]

# read the files to a methylRawList object: myobj
in_treatment <- eval(parse(text = paste("c(",input_treatment,")")))
if(length(file.list) != length(in_treatment)) error_message <-"the number of samples is different from the number of treatment conditions specified"

{{
myobj <- read(file.list, 
              sample.id = as.list(paste("condition",in_treatment,sep="")),
              assembly = "",# input_assembly,         #https://groups.google.com/forum/#!topic/methylkit_discussion/keXZoQlbups
              treatment = in_treatment, 
              context = input_context,
              pipeline = input_pipeline)          
}}

### 2.quality check and basic features of the data
getMethylationStats(myobj[[2]], plot = T, both.strands = F)  # plots the histogram for percent methylation distribution
# Typically, percent methylation histogram should have two peaks on both ends.

getCoverageStats(myobj[[2]], plot = T, both.strands = F)   # plot the read coverage per base information
# Experiments that are highly suffering from PCR duplication bias will have a secondary peak towards the right hand side of the histogram.

print(capture.output(cat("\nquality check: percent methylation statistics for second sample\n")),collapse='\n')
print(capture.output(getMethylationStats(myobj[[2]], plot = F, both.strands = F)),collapse='\n')

### 3.filter samples based on read coverage
# filters a methylRawList and discards bases that havecoverage below 10X,
# and also discards the bases that have more than 99.9th percentile of coverage ineach sample.
dimbefore <- NULL
for(i in 1:length(in_treatment)) dimbefore[i] <- paste(dim(myobj[[i]]),collapse=" ")
print(capture.output(cat("\nthe dimension of the data:\n")),collapse='\n')
print(capture.output(dimbefore),collapse='\n')
{{
  if(input_filter == TRUE) {
    if (input_locount == 0) input_locount = NULL
    if (input_loperc == 0) input_loperc = NULL
    if (input_hicount == 0) input_hicount = NULL
    if (input_hiperc == 0) input_hiperc = NULL  
    myobj <- filterByCoverage(myobj, lo.count = input_locount,lo.perc = input_loperc, hi.count = input_hicount, hi.perc = input_hiperc)
    dimafter <- NULL
    for(i in 1:length(in_treatment)) dimafter[i] <- paste(dim(myobj[[i]]),collapse=" ")
    print(capture.output(cat("\nafter filtering, the dimension of the data:\n")),collapse='\n')
    print(capture.output(dimafter),collapse='\n')
  }
}}
## (optional) get differentially methylated regions
# For some situations, it might be desirable to summarize methylation information over tiling windows
# or over a set of predefined regions (promoters, CpG islands, introns, etc.)
# tiles the genome with windows 1000bp length and 1000bp step-size and summarizes the methylation information on those tiles
{{
  if (input_tilingwindow == TRUE) {
    myobj <- tileMethylCounts(myobj, win.size = input_winsize, step.size = input_stepsize, cov.bases = input_covbase)
  }
}}

### 4.sample correlation
# merge all samples to one object for base-pair locations that are covered in all samples. 
# return a methylBase object
# Setting destrand=TRUE will merge reads on both strands of a CpG dinucleotide. 
# This provides better coverage, but only advised when looking at CpG methylation. 
# In addition, setting destrand=TRUE will only work when operating on base-pair resolution
meth <- unite(myobj, destrand = input_destrand)
print(capture.output(cat("\nsample correlation:\n")),collapse='\n')
print(capture.output(print(getCorrelation(meth, method = input_cormethod, plot = F))),collapse='\n')


### 5.cluster samples based on similarity of methylation profiles
# cluster (hierarchical) the samples and draw a dendrogram
clusterSamples(meth, dist = input_dist, method = input_clustermethod, sd.filter=input_sdfilter,sd.threshold=input_sdthreshold, filterByQuantile=input_filterbyquantile, plot = TRUE)
# plot a scree plot for importance of components.
#PCASamples(meth, screeplot = TRUE)
# plot PC1 and PC2 axis and a scatter plot of our samples on those axis which will reveal how they cluster.
PCASamples(meth, scale = input_scale, center = input_center, transpose = input_transpose, sd.filter=input_sdfilter,sd.threshold=input_sdthreshold, filterByQuantile=input_filterbyquantile)


### 6.get differentially methylated bases
# calculate differential methylation 
# Depending on the sample size per each set it will either use Fisher's exact or logistic regression to calculate Pvalues. 
# P-values will be adjusted to Q-values.
# for machines with multiple cores, use multiple-cores to increse speed for both fisher"s easct test and logistic regression based test
myDiff <- calculateDiffMeth(meth, slim = input_slim, weighted.mean = input_weightedmean) #, num.cores = 2)  
# select the differentially methylated regions/bases based on q-value and percent methylation difference cutoffs. 
# here selects the bases that have qvalue<0.01 and percent methylation difference larger than 25%. 
# If you specify type="hyper" or type="hypo" options, you will get hyper-methylated or hypo-methylated regions/bases.
#
# get hyper methylated bases
myDiff.hyper <- get.methylDiff(myDiff, difference = input_difference, qvalue = input_qvalue, type = "hyper")
#
# get hypo methylated bases
myDiff.hypo <- get.methylDiff(myDiff, difference = input_difference, qvalue = input_qvalue, type = "hypo")
#
# get all differentially methylated bases
myDiff.de <- get.methylDiff(myDiff, difference = input_difference,qvalue = input_qvalue, type = "all")

output_all <- cbind(getData(meth),getData(myDiff)[,c("pvalue","qvalue","meth.diff")])
output_hyper <- cbind("hyper",merge(output_all,getData(myDiff.hyper),by = c("chr","start","end","strand","pvalue","qvalue","meth.diff")))
colnames(output_hyper)[1] <- "type"
output_hypo <- cbind("hypo",merge(output_all,getData(myDiff.hypo),by = c("chr","start","end","strand","pvalue","qvalue","meth.diff")))
colnames(output_hypo)[1] <- "type"
output_diff <- rbind(output_hyper,output_hypo)
output_all <- merge(output_all,output_diff[,c("chr","start","end","strand","pvalue","qvalue","meth.diff","type")],by = c("chr","start","end","strand","pvalue","qvalue","meth.diff"), all.x = TRUE,stringsAsFactors=FALSE)

type1 <- rep("",dim(output_all)[1])
output_all <- cbind(output_all, type1,stringsAsFactors=FALSE)
{{
  for(i in 1:dim(output_all)[1]) {
    if (is.na(output_all[i,"type"])) {
      output_all[i,"type1"] = "none"
    }else{
      if(output_all[i,"type"] == "hyper") {
        output_all[i,"type1"] = "hyper"
      }else{
        if(output_all[i,"type"] == "hypo") output_all[i,"type1"] = "hypo"
      }
    }
  }
}}
output_all <- subset(output_all, select = -c(type))
names(output_all)[names(output_all)=="type1"] <- "type"


### 7. differential methylation events per chr
#  visualize the distribution of hypo/hyper-methylated bases/regions per chr
# The list shows percentages of hypo/hyper methylated bases over all the covered bases in a given chromosome.
# if plot=FALSE a list having per chromosome differentially methylation events will be returned
pctperchr <- diffMethPerChr(myDiff, plot = FALSE, qvalue.cutoff = input_qvalue, meth.cutoff = input_difference)
# if plot=TRUE a barplot will be plotted
diffMethPerChr(myDiff, plot = TRUE, qvalue.cutoff = input_qvalue, meth.cutoff = input_difference)
print(capture.output(cat("\npercentages of hypo/hyper methylated bases over all the covered bases in a given chromosome:\n")),collapse='\n')
print(capture.output(print(pctperchr)),collapse='\n')

### 8. annoate differntial methylation events
# what percentage of our differentially methylated regions are on promoters/introns/exons/intergenic region.
if (input_ann == TRUE) print(capture.output(cat("\nannotate differentially methylated Cs with promoter/exon/intron using annotation data:\n")),collapse='\n')
{{
  if(input_ann == TRUE){
    gene.obj <- read.transcript.features(input_annotation)
    # annotate differentially methylated Cs with promoter/exon/intron using annotation data
    ann <- annotate.WithGenicParts(myDiff.de, gene.obj)
    #print(capture.output(cat("\nannotate differentially methylated Cs with promoter/exon/intron using annotation data:\n")),collapse='\n')
    print(capture.output(print(ann)),collapse='\n')
  }
}}

{{
  if(input_ann == TRUE){
    # get the distance to TSS (transcription start site) and nearest gene name
    diffAnn <- annotate.WithGenicParts(myDiff.de, gene.obj)
    # target.row is the row number in myDiff
    #head(getAssociationWithTSS(diffAnn), 3)
    # the annotation table which shows overlap status of the methylation events with the annotation such as promoter/exon/intron
    #head(getMembers(diffAnn))
    output_anndiff <- cbind(output_diff,getAssociationWithTSS(diffAnn),getMembers(diffAnn))
    
    # get percentage/number of differentially methylated regions that overlap with intron/exon/promoters
    pctregion <- getTargetAnnotationStats(diffAnn, percentage = TRUE,precedence = TRUE)
    
    # plot the percentage of differentially methylated bases overlapping with exon/intron/promoters
    plotTargetAnnotation(diffAnn, precedence = TRUE, main = "differential methylation annotation")
    print(capture.output(cat("\npercentage/number of differentially methylated regions that overlap with intron/exon/promoters:\n")),collapse='\n')
  }
}}
if(input_ann == TRUE) print(capture.output(print(pctperchr)),collapse='\n')

# read the CpG island annotation and annotate our differentially methylated bases/regions with them.
# read the shores and flanking regions and name the flanks as shores and CpG islands as CpGi
if(input_cpgann == TRUE) print(capture.output(cat("\nread the CpG island annotation and annotate our differentially methylated bases/regions.\nread the shores and flanking regions and name the flanks as shores and CpG islands as CpGi:\n")),collapse='\n')
{{
  if(input_cpgann == TRUE){
    cpg.obj <- read.feature.flank(input_cpgannotation, feature.flank.name = c("CpGi", "shores"),remove.unsual=input_removeunsual,flank=input_flank,clean=input_clean)
    diffCpGann <- annotate.WithFeature.Flank(myDiff.de, cpg.obj$CpGi,cpg.obj$shores, feature.name = "CpGi", flank.name = "shores", strand=input_strand)
    print(capture.output(print(diffCpGann)),collapse='\n')    
  }
}}

{{
  if(input_cpgann == TRUE){
    # plot the percentage of differentially methsylated bases are on CpG islands, CpG island shores and other regions.
    plotTargetAnnotation(diffCpGann, col = c("green", "gray","white"), main = "differential methylation annotation")
    # percentage of intron/exon/promoters that overlap with differentially methylated bases.
    pctint <- getFeatsWithTargetsStats(diffAnn, percentage = TRUE)
    print(capture.output(cat("\npercentage of intron/exon/promoters that overlap with differentially methylated bases:\n")),collapse='\n')
  }
}}
if(input_cpgann == TRUE) print(capture.output(print(pctint)),collapse='\n')



