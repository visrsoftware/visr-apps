source("visrutils.R")

visr.biocLite("chromPlot")
visr.biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
visr.biocLite("TxDb.Mmusculus.UCSC.mm10.knownGene")
visr.biocLite("GenomicFeatures")

visr.app.start("Chrom Plot 4.0")

visr.category("Data Visibility")

visr.param("plottitle", type = "char", label = "Plot Title", default = "Chromosome Plot")
visr.param("humangenome", type = "char", default = "hg19 Human Genome",
           items = c("hg19 Human Genome", "mm10 Mouse Genome"), label = "Genome")
visr.param("bandsorideogram", type = "char", default = "Bands",
           items = c("Bands", "Ideogram", "None"), label = "Data On Chromosomes")
visr.param("chromside", type = "char", default = "Gene Transcripts",
           items = c("Gene Transcripts", "Statistics", "Statistics & Gene Transcripts", "Gene IDs", "None"),
           label = "Data On Sides Of Chromosomes",
           info = "Please select Statistics data column first if choosing Statistics option. For Gene IDs, small dataset is preferable.")
visr.param("datacolor", type = "color", default = "seagreen", label = "Gene Transcripts Color")
visr.param("textsizepercentage", type = "int", default = 100, min = 50, max = 200, label = "Text Size Percentage")

visr.category("Chromosomes")

visr.param("chromosome1", type = "logical", default = TRUE, label = "Chromosome 1")
visr.param("chromosome2", type = "logical", default = TRUE, label = "Chromosome 2")
visr.param("chromosome3", type = "logical", default = TRUE, label = "Chromosome 3")
visr.param("chromosome4", type = "logical", default = TRUE, label = "Chromosome 4")
visr.param("chromosome5", type = "logical", default = TRUE, label = "Chromosome 5")
visr.param("chromosome6", type = "logical", default = TRUE, label = "Chromosome 6")
visr.param("chromosome7", type = "logical", default = TRUE, label = "Chromosome 7")
visr.param("chromosome8", type = "logical", default = TRUE, label = "Chromosome 8")
visr.param("chromosome9", type = "logical", default = TRUE, label = "Chromosome 9")
visr.param("chromosome10", type = "logical", default = TRUE, label = "Chromosome 10")
visr.param("chromosome11", type = "logical", default = TRUE, label = "Chromosome 11")
visr.param("chromosome12", type = "logical", default = TRUE, label = "Chromosome 12")
visr.param("chromosome13", type = "logical", default = TRUE, label = "Chromosome 13")
visr.param("chromosome14", type = "logical", default = TRUE, label = "Chromosome 14")
visr.param("chromosome15", type = "logical", default = TRUE, label = "Chromosome 15")
visr.param("chromosome16", type = "logical", default = TRUE, label = "Chromosome 16")
visr.param("chromosome17", type = "logical", default = TRUE, label = "Chromosome 17")
visr.param("chromosome18", type = "logical", default = TRUE, label = "Chromosome 18")
visr.param("chromosome19", type = "logical", default = TRUE, label = "Chromosome 19")
visr.param("chromosome20", type = "logical", default = TRUE, label = "Chromosome 20 (Human Only)")
visr.param("chromosome21", type = "logical", default = TRUE, label = "Chromosome 21 (Human Only)")
visr.param("chromosome22", type = "logical", default = TRUE, label = "Chromosome 22 (Human Only)")
visr.param("chromosomeX", type = "logical", default = TRUE, label = "Chromosome X")
visr.param("chromosomeY", type = "logical", default = TRUE, label = "Chromosome Y")

visr.category("Bands Input")

visr.param("chromno", type = "column", label = "Chromosome Number")
visr.param("start", type = "column-numerical", label = "Gene Start")
visr.param("end", type = "column-numerical", label = "Gene End")
visr.param("id", type = "column", label = "Gene ID")
visr.param("bandcolor", type = "color", default = "firebrick1", label = "Band Color")
visr.param("strandsvisible", type = "logical", default = FALSE, label = "Bands Colored By Strand +/-")
visr.param("strandsinput", type = "column", label = "Strands")
visr.param("strandscolor", type = "color", default = "blue3", label = "Band Color (-)")

visr.category("Statistics")

visr.param("statline", type = "logical", default = FALSE, label = "Lines")
visr.param("statmean", type = "logical", default = FALSE, label = "Mean")
visr.param("statcol", type = "column-numerical", label = "Statistics Data")
visr.param("statcolor", type = "color", default = "darkorchid", label = "Statistics Color")

visr.app.end(printjson = TRUE, writefile = TRUE)
visr.applyParameters()

# Loading Gaps & Ideogram Data
data(hg_cytoBandIdeo)
data(mm10_cytoBandIdeo)
data(hg_gap)
data(mm10_gap)

if (visr.param.humangenome == "hg19 Human Genome") {
  
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  txgr <- transcripts(txdb)
  
  whichgap <- hg_gap
  whichcyto <- hg_cytoBandIdeo
  
} else {
  
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  txgr <- transcripts(txdb)
  
  whichgap <- mm10_gap
  whichcyto <- mm10_cytoBandIdeo
  
}

# Loading Bands Input Data Into Data Frame
our.frame <- data.frame(Chrom=visr.input[[visr.param.chromno]],
                        Start=visr.input[[visr.param.start]],
                        End=visr.input[[visr.param.end]],
                        Name=visr.input[[visr.param.id]])

# Adding Statistics Column To Data Frame
if (visr.param.chromside == "Statistics" || visr.param.chromside == "Statistics & Gene Transcripts") {
  our.frame$StatCol <- visr.input[[visr.param.statcol]]
  if (visr.param.statline == TRUE) {
    statisticstype <- "l"
  } else {
    statisticstype <- "p"
  }
  if (visr.param.statmean == TRUE) {
    statisticsmean <- "mean"
  } else {
    statisticsmean <- "none"
  }
}

# Adding Band Colors To Data Frame
our.frame$Colors <- visr.param.bandcolor
if (visr.param.strandsvisible == TRUE) {
  strandscolumn <- visr.input[[visr.param.strandsinput]]
  for (i in 1:length(strandscolumn)) {
    if (strandscolumn[i] == "+") {
      strandscolumn[i] <- visr.param.bandcolor
    } else {
      strandscolumn[i] <- visr.param.strandscolor
    }
  }
  our.frame$Colors <- strandscolumn
}

# Adding ID Column To Data Frame
if (visr.param.chromside == "Gene IDs") {
  our.frame$ID <- visr.input[[visr.param.id]]
}

# Creating Chromosomes Visibility Vector
chromosome.parameters = c(visr.param.chromosome1, visr.param.chromosome2, visr.param.chromosome3,
                          visr.param.chromosome4, visr.param.chromosome5, visr.param.chromosome6,
                          visr.param.chromosome7, visr.param.chromosome8, visr.param.chromosome9,
                          visr.param.chromosome10, visr.param.chromosome11, visr.param.chromosome12,
                          visr.param.chromosome13, visr.param.chromosome14, visr.param.chromosome15,
                          visr.param.chromosome16, visr.param.chromosome17, visr.param.chromosome18,
                          visr.param.chromosome19, visr.param.chromosome20, visr.param.chromosome21,
                          visr.param.chromosome22, visr.param.chromosomeX, visr.param.chromosomeY)

chromosomes <- vector(mode = "character", length = 0)

for (i in 1:24) {
  if (chromosome.parameters[i] == TRUE &&
      ((visr.param.humangenome == "hg19 Human Genome" && i < 23)
       || (visr.param.humangenome == "mm10 Mouse Genome" && i < 20))) {
    chromosomes <- cbind(chromosomes, as.character(i))
  }
  else if (chromosome.parameters[i] == TRUE && i == 23) {
    chromosomes <- cbind(chromosomes, "X")
  }
  else if (chromosome.parameters[i] == TRUE && i == 24) {
    chromosomes <- cbind(chromosomes, "Y")
  }
}

# Creating Positioning Vector
side.vector <- c(-1, -1, -1, -1, -1, -1, -1, -1)
if (visr.param.chromside == "Statistics & Gene Transcripts") {
  side.vector <- c(-1, -1, -1, -1, -1, -1, 1, -1)
}

# Columns Number
if (length(chromosomes) > 12) {
  no.columns <- 12
} else {
  no.columns <- length(chromosomes)
}

# Bands Or Ideogram
if (visr.param.bandsorideogram == "Bands") {
  whichcyto <- our.frame
}

# Text Size
textsize <- visr.param.textsizepercentage / 100

# Calling ChromPlot
if (!(visr.param.bandsorideogram == "None")) {
  if (visr.param.chromside == "Gene Transcripts") {
    chromPlot(gaps=whichgap, bands=whichcyto, annot1=txgr, chr=chromosomes,
              figCols=no.columns, chrSide=side.vector, colAnnot1=visr.param.datacolor,
              title=visr.param.plottitle, cex=textsize)
  } else if (visr.param.chromside == "Statistics") {
    chromPlot(gaps=whichgap, bands=whichcyto, chr=chromosomes,
              figCols=no.columns, chrSide=side.vector,
              stat=our.frame, statCol="StatCol", statTyp=statisticstype, statSumm=statisticsmean,
              colStat=visr.param.statcolor, title=visr.param.plottitle, cex=textsize)
  } else if (visr.param.chromside == "Statistics & Gene Transcripts") {
    chromPlot(gaps=whichgap, bands=whichcyto, annot1=txgr, chr=chromosomes,
              stat=our.frame, statCol="StatCol", statTyp=statisticstype, statSumm=statisticsmean, noHist=TRUE,
              statName="StatCol", colStat=visr.param.statcolor,
              figCols=no.columns, chrSide=side.vector, colAnnot1=visr.param.datacolor,
              title=visr.param.plottitle, cex=textsize)
  } else if (visr.param.chromside == "Gene IDs") {
    chromPlot(gaps=whichgap, bands=whichcyto, chr=chromosomes,
              figCols=no.columns, chrSide=side.vector, cex=textsize,
              stat=our.frame, statCol="Value", statName="Value", noHist=TRUE, statTyp="n",
              colStat=visr.param.datacolor)
  } else {
    chromPlot(gaps=whichgap, bands=whichcyto, chr=chromosomes,
              figCols=no.columns, chrSide=side.vector,
              title=visr.param.plottitle, cex=textsize)
  }
} else {
  if (visr.param.chromside == "Gene Transcripts") {
    chromPlot(gaps=whichgap, annot1=txgr, chr=chromosomes,
              figCols=no.columns, chrSide=side.vector, colAnnot1=visr.param.datacolor,
              title=visr.param.plottitle, cex=textsize)
  } else if (visr.param.chromside == "Statistics") {
    chromPlot(gaps=whichgap, chr=chromosomes,
              figCols=no.columns, chrSide=side.vector,
              stat=our.frame, statCol="StatCol", statTyp=statisticstype, statSumm=statisticsmean,
              colStat=visr.param.statcolor, title=visr.param.plottitle, cex=textsize)
  } else if (visr.param.chromside == "Statistics & Gene Transcripts") {
    chromPlot(gaps=whichgap, annot1=txgr, chr=chromosomes,
              stat=our.frame, statCol="StatCol", statTyp=statisticstype, statSumm=statisticsmean, noHist=TRUE,
              statName="StatCol", colStat=visr.param.statcolor,
              figCols=no.columns, chrSide=side.vector, colAnnot1=visr.param.datacolor,
              title=visr.param.plottitle, cex=textsize)
  } else if (visr.param.chromside == "Gene IDs") {
    chromPlot(gaps=whichgap, chr=chromosomes,
              figCols=no.columns, chrSide=side.vector, cex=textsize,
              stat=our.frame, statCol="Value", statName="Value", noHist=TRUE, statTyp="n",
              colStat=visr.param.datacolor)
  } else {
    chromPlot(gaps=whichgap, chr=chromosomes,
              figCols=no.columns, chrSide=side.vector,
              title=visr.param.plottitle, cex=textsize)
  }
}