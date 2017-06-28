source("visrutils.R")

visr.biocLite("DMRcaller")

visr.app.start("DMRcaller")

visr.category("Data Files")

visr.param("cxreportfilewt", type = "filename", label = "CX Report File For Wild Type")
visr.param("cxreportfilemutant", type = "filename", label = "CX Report File For Mutant")

visr.category("Plot Options")

visr.param(
  "plot",
  type = "char",
  label = "Plot Type",
  items = c(
    "Low Resolution Profile",
    "Coverage",
    "DMR Distribution",
    "DMR Profiles"
  )
)

visr.param("plottitle",
           type = "char",
           label = "Plot Title",
           default = "Methylation")

visr.param("titlewt",
           type = "char",
           label = "Label For WT",
           default = "WT")
visr.param("titlemutant",
           type = "char",
           label = "Label For Mutant",
           default = "Mutant")

visr.param(
  "context",
  type = "char",
  label = "Context",
  items = c("CG", "CHH", "CHG")
)

visr.category("Selected Data Options")

visr.param(
  "chromosome",
  type = "char",
  default = "Chr3",
  label = "Region of Interest Chromosome",
  info = "Must correspond to CX Report file chromosome number column."
)
visr.param(
  "locationstart",
  type = "integer",
  default = 0,
  min = 0,
  label = "Region of Interest Start Location"
)
visr.param(
  "locationend",
  type = "integer",
  default = 1000000,
  min = 0,
  label = "Region of Interest End Location",
  info = "Warning: May be slow or crash (CPU intensive) with large >300Kb region when plotting DMR Profiles."
)

visr.category("Visual Options")

visr.param("color1",
           type = "color",
           label = "Plot Color For WT",
           default = "#006A6A")
visr.param("color2",
           type = "color",
           label = "Plot Color For Mutant",
           default = "#E69F00")

visr.param(
  "point1",
  type = "integer",
  min = 0,
  max = 25,
  label = "Point Type For WT",
  default = 1
)
visr.param(
  "point2",
  type = "integer",
  min = 0,
  max = 25,
  label = "Point Type For Mutant",
  default = 0
)

visr.param(
  "line1",
  type = "integer",
  min = 1,
  max = 6,
  label = "Line Type For WT",
  default = 4
)
visr.param(
  "line2",
  type = "integer",
  min = 1,
  max = 6,
  label = "Line Type For Mutant",
  default = 1
)

visr.app.end(printjson = TRUE, writefile = TRUE)
visr.applyParameters()

#####

# load presaved data
data(methylationDataList)

# load your data
methylationDataWT <- readBismark(visr.param.cxreportfilewt)
methylationDataMutant <- readBismark(visr.param.cxreportfilemutant)
methylationDataList <- GRangesList(WT = methylationDataWT,
                                   Mutant = methylationDataMutant)

#####

# plots

if (visr.param.plot == "Low Resolution Profile") {
  chr_local <- GRanges(
    seqnames = Rle(visr.param.chromosome),
    ranges = IRanges(visr.param.locationstart, visr.param.locationend)
  )
  
  par(mar = c(4, 4, 3, 1) + 0.1)
  plotMethylationProfileFromData(
    methylationDataList[["WT"]],
    methylationDataList[["Mutant"]],
    conditionsNames = c(visr.param.titlewt, visr.param.titlemutant),
    regions = chr_local,
    windowSize = 10000,
    autoscale = FALSE,
    context = c(visr.param.context),
    col = c(visr.param.color1, visr.param.color2),
    pch = c(visr.param.point1, visr.param.point2),
    lty = c(visr.param.line1, visr.param.line2)
  )
  
} else if (visr.param.plot == "Coverage") {
  chr_local <- GRanges(
    seqnames = Rle(visr.param.chromosome),
    ranges = IRanges(visr.param.locationstart, visr.param.locationend)
  )
  
  # plot the coverage in the two contexts
  par(mar = c(4, 4, 3, 1) + 0.1)
  plotMethylationDataCoverage(
    methylationDataList[["WT"]],
    methylationDataList[["Mutant"]],
    breaks = c(1, 5, 10, 15),
    regions = chr_local,
    conditionsNames = c(visr.param.titlewt, visr.param.titlemutant),
    context = c(visr.param.context),
    proportion = TRUE,
    labels = LETTERS,
    contextPerRow = FALSE,
    col = c(visr.param.color1, visr.param.color2),
    pch = c(visr.param.point1, visr.param.point2),
    lty = c(visr.param.line1, visr.param.line2)
  )
  
} else if (visr.param.plot == "DMR Distribution") {
  chr_local <- GRanges(
    seqnames = Rle(visr.param.chromosome),
    ranges = IRanges(visr.param.locationstart, visr.param.locationend)
  )
  
  # compute the DMRs in CG context with noise_filter method
  DMRsNoiseFilterCG <- computeDMRs(
    methylationDataList[["WT"]],
    methylationDataList[["Mutant"]],
    regions = chr_local,
    context = visr.param.context,
    method = "noise_filter",
    windowSize = 100,
    kernelFunction = "triangular",
    test = "score",
    pValueThreshold = 0.01,
    minCytosinesCount = 4,
    minProportionDifference = 0.4,
    minGap = 0,
    minSize = 50,
    minReadsPerCytosine = 4,
    cores = 1
  )
  
  DMRsNoiseFilterCGMerged <- mergeDMRsIteratively(
    DMRsNoiseFilterCG,
    minGap = 200,
    respectSigns = TRUE,
    methylationDataList[["WT"]],
    methylationDataList[["Mutant"]],
    context = visr.param.context,
    minProportionDifference = 0.4,
    minReadsPerCytosine = 4,
    pValueThreshold = 0.01,
    test = "score"
  )
  # compute the distribution of DMRs
  hotspots <-
    computeOverlapProfile(DMRsNoiseFilterCGMerged,
                          chr_local,
                          windowSize = 5000,
                          binary = TRUE)
  # plot the distribution of DMRs
  plotOverlapProfile(
    GRangesList("Chromosome" = hotspots),
    col = c(visr.param.color1, visr.param.color2),
    title = visr.param.plottitle
  )
  
} else if (visr.param.plot == "DMR Profiles") {
  chr_local <- GRanges(
    seqnames = Rle(visr.param.chromosome),
    ranges = IRanges(visr.param.locationstart, visr.param.locationend)
  )
  
  # compute the DMRs in CG context with noise_filter method
  DMRsNoiseFilterCG <- computeDMRs(
    methylationDataList[["WT"]],
    methylationDataList[["Mutant"]],
    regions = chr_local,
    context = visr.param.context,
    method = "noise_filter",
    windowSize = 100,
    kernelFunction = "triangular",
    test = "score",
    pValueThreshold = 0.01,
    minCytosinesCount = 4,
    minProportionDifference = 0.4,
    minGap = 0,
    minSize = 50,
    minReadsPerCytosine = 4,
    cores = 1
  )
  
  # compute the DMRs in CG context with neighbourhood method
  DMRsNeighbourhoodCG <- computeDMRs(
    methylationDataList[["WT"]],
    methylationDataList[["Mutant"]],
    regions = chr_local,
    context = visr.param.context,
    method = "neighbourhood",
    test = "score",
    pValueThreshold = 0.01,
    minCytosinesCount = 4,
    minProportionDifference = 0.4,
    minGap = 200,
    minSize = 1,
    minReadsPerCytosine = 4,
    cores = 1
  )
  
  # compute the DMRs in CG context with bins method
  DMRsBinsCG <- computeDMRs(
    methylationDataList[["WT"]],
    methylationDataList[["Mutant"]],
    regions = chr_local,
    context = visr.param.context,
    method = "bins",
    binSize = 100,
    test = "score",
    pValueThreshold = 0.01,
    minCytosinesCount = 4,
    minProportionDifference = 0.4,
    minGap = 200,
    minSize = 50,
    minReadsPerCytosine = 4,
    cores = 1
  )
  
  # load the gene annotation data
  data(GEs)
  #select the genes
  genes <- GEs[which(GEs$type == "gene")]
  # compute the DMRs in CG context over genes
  DMRsGenesCG <- filterDMRs(
    methylationDataList[["WT"]],
    methylationDataList[["Mutant"]],
    potentialDMRs = genes[overlapsAny(genes, chr_local)],
    context = visr.param.context,
    test = "score",
    pValueThreshold = 0.01,
    minCytosinesCount = 4,
    minProportionDifference = 0.4,
    minReadsPerCytosine = 3,
    cores = 1
  )
  
  DMRsNoiseFilterCGMerged <- mergeDMRsIteratively(
    DMRsNoiseFilterCG,
    minGap = 200,
    respectSigns = TRUE,
    methylationDataList[["WT"]],
    methylationDataList[["Mutant"]],
    context = visr.param.context,
    minProportionDifference = 0.4,
    minReadsPerCytosine = 4,
    pValueThreshold = 0.01,
    test = "score"
  )
  
  # select a location on the chromosome
  chr3Reg <- GRanges(
    seqnames = Rle(visr.param.chromosome),
    ranges = IRanges(visr.param.locationstart, visr.param.locationend)
  )
  
  # create a list with all DMRs
  DMRsCGList <- list(
    "noise filter" = DMRsNoiseFilterCGMerged,
    "neighbourhood" = DMRsNeighbourhoodCG,
    "bins" = DMRsBinsCG,
    "genes" = DMRsGenesCG
  )
  
  # plot the local profile
  par(cex = 0.9)
  par(mar = c(4, 4, 3, 1) + 0.1)
  plotLocalMethylationProfile(
    methylationDataList[["WT"]],
    methylationDataList[["Mutant"]],
    chr3Reg,
    DMRsCGList,
    conditionsNames = c(visr.param.titlewt, visr.param.titlemutant),
    GEs,
    windowSize = 300,
    main = visr.param.plottitle,
    context = visr.param.context
  )
}