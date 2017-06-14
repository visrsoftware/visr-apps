source("visrutils.R")

visr.biocLite("methylKit")
visr.biocLite("genomation")

visr.app.start("Methyl Kit")

visr.category("Data Files")

visr.param(
  "methcallorbismark",
  type = "char",
  label = "File Type",
  items = c("Methylation Call", "Bismark")
)
visr.param(
  "onefileormulti",
  type = "char",
  label = "Number of Files",
  items = c("Multiple Files", "One File")
)
visr.param("file1",
           type = "filename",
           filename.mode = "load",
           label = "File 1")
visr.param("file2",
           type = "filename",
           filename.mode = "load",
           label = "File 2")
visr.param("file3",
           type = "filename",
           filename.mode = "load",
           label = "File 3")
visr.param("file4",
           type = "filename",
           filename.mode = "load",
           label = "File 4")
visr.param("filename1",
           type = "char",
           label = "File 1 Title",
           default = "Test 1")
visr.param("filename2",
           type = "char",
           label = "File 2 Title",
           default = "Test 2")
visr.param("filename3",
           type = "char",
           label = "File 3 Title",
           default = "Ctrl 1")
visr.param("filename4",
           type = "char",
           label = "File 4 Title",
           default = "Ctrl 2")

visr.category("Output Type")

visr.param(
  "outputtype",
  type = "char",
  label = "Output Plot Type",
  items = c(
    "Methylation",
    "Coverage",
    "Correlation",
    "Clustering",
    "PCA",
    "Differential Methylation Annotation"
  ),
  item.labels = c(
    "Methylation",
    "Coverage",
    "Correlation (MultiFile Only)",
    "Clustering (MultiFile MethCall Only)",
    "PCA (MultiFile MethCall Only)",
    "Differential Methylation Annotation (MultiFile MethCall Only)"
  )
)

visr.category("Methylation/Coverage Options")

visr.param(
  "filechoose",
  type = "int",
  min = 1,
  max = 4,
  default = 1,
  label = "Choose file number for plots"
)
visr.param("bothstrands",
           type = "logical",
           default = FALSE,
           label = "Both strands")

visr.category("Correlation/Clustering/PCA Options")

visr.param("methmin",
           type = "logical",
           label = "CpGs with >=1 sample per group only",
           default = FALSE)
visr.param(
  "clusteringdist",
  type = "char",
  default = "correlation",
  label = "Clustering Distance Measure",
  items = c(
    "correlation",
    "euclidean",
    "maximum",
    "manhattan",
    "binary",
    "minkowski"
  )
)
visr.param(
  "clusteringmethod",
  type = "char",
  default = "ward",
  label = "Clustering Method",
  items = c(
    "ward",
    "single",
    "complete",
    "average",
    "mcquitty",
    "median",
    "centroid"
  )
)
visr.param(
  "screeoraxis",
  type = "char",
  default = "PCA Scree Plot",
  label = "PCA Plot Type",
  items = c("PCA Scree Plot", "PCA Axis Plot")
)

visr.category("Differential Methylation Annotation Options")

visr.param(
  "overlap",
  type = "char",
  default = "Exons/Introns/Promoters",
  label = "Overlap With",
  items = c("Exons/Introns/Promoters", "CpG Islands")
)
visr.param(
  "difference",
  type = "int",
  min = 0,
  max = 100,
  default = 25,
  label = "Percent Methylation Difference"
)
visr.param(
  "qvalue",
  type = "double",
  min = 0,
  max = 1,
  default = 0.01,
  label = "q-value"
)

visr.category("Filter By Coverage")

visr.param("filter",
           type = "logical",
           label = "Filter samples based on read coverage",
           default = FALSE)
visr.param(
  "locount",
  type = "int",
  min = 0,
  max = 100,
  default = 10,
  label = "Coverage lower limit in count"
)
visr.param(
  "loperc",
  type = "double",
  min = 0,
  max = 100,
  default = 0,
  label = "Coverage lower limit in percentile"
)
visr.param(
  "hiperc",
  type = "double",
  min = 0,
  max = 100,
  default = 99.9,
  label = "Coverage upper limit in percentile"
)

visr.app.end(printjson = TRUE, writefile = TRUE)
visr.applyParameters()

# Reading files into list
if (visr.param.methcallorbismark == "Methylation Call" &&
    visr.param.onefileormulti == "Multiple Files") {
  file.list = list(visr.param.file1,
                   visr.param.file2,
                   visr.param.file3,
                   visr.param.file4)
  
  # Read the files to a methylRawList object: myobj
  myobj = methRead(
    file.list,
    sample.id = list(
      visr.param.filename1,
      visr.param.filename2,
      visr.param.filename3,
      visr.param.filename4
    ),
    assembly = "hg18",
    treatment = c(1, 1, 0, 0),
    context = "CpG"
  )
  
} else if (visr.param.methcallorbismark == "Methylation Call" &&
           visr.param.onefileormulti == "One File") {
  # Read the files to a methylRaw object: myobj
  myobj = methRead(
    visr.param.file1,
    sample.id = visr.param.filename1,
    assembly = "hg18",
    treatment = c(1, 1, 0, 0),
    context = "CpG"
  )
  
} else if (visr.param.methcallorbismark == "Bismark" &&
           visr.param.onefileormulti == "One File") {
  # Read the files to a methylRaw object: myobj
  myobj = processBismarkAln(
    location = visr.param.file1,
    sample.id = visr.param.filename1,
    assembly = "hg18",
    read.context = "CpG",
    save.folder = getwd()
  )
  
} else if (visr.param.methcallorbismark == "Bismark" &&
           visr.param.onefileormulti == "Multiple Files") {
  file.list = list(visr.param.file1,
                   visr.param.file2,
                   visr.param.file3,
                   visr.param.file4)
  
  # Read the files to a methylRawList object: myobj
  myobj = processBismarkAln(
    location = file.list,
    sample.id = list(
      visr.param.filename1,
      visr.param.filename2,
      visr.param.filename3,
      visr.param.filename4
    ),
    assembly = "hg18",
    read.context = "CpG",
    save.folder = getwd(),
    treatment = c(1, 1, 0, 0)
  )
}

if (visr.param.filter == TRUE) {
  myobj = filterByCoverage(
    myobj,
    lo.count = visr.param.locount,
    lo.perc = visr.param.loperc,
    hi.count = NULL,
    hi.perc = visr.param.hiperc
  )
}

if (visr.param.onefileormulti == "Multiple Files") {
  meth = unite(myobj, destrand = FALSE)
  if (visr.param.methmin == TRUE) {
    meth = unite(myobj, min.per.group = 1L)
  }
} else {
  meth = myobj
}

if (visr.param.screeoraxis == "PCA Scree Plot") {
  scree = TRUE
} else {
  scree = FALSE
}

if (visr.param.onefileormulti == "Multiple Files") {
  if (visr.param.outputtype == "Methylation") {
    getMethylationStats(myobj[[visr.param.filechoose]], plot = TRUE, both.strands =
                          visr.param.bothstrands)
  } else if (visr.param.outputtype == "Coverage") {
    getCoverageStats(myobj[[visr.param.filechoose]], plot = TRUE, both.strands =
                       visr.param.bothstrands)
  } else if (visr.param.outputtype == "Correlation") {
    getCorrelation(meth, plot = TRUE)
  } else if (visr.param.outputtype == "Clustering" &&
             visr.param.methcallorbismark == "Methylation Call") {
    clusterSamples(meth,
                   dist = visr.param.clusteringdist,
                   method = visr.param.clusteringmethod,
                   plot = TRUE)
  } else if (visr.param.outputtype == "PCA" &&
             visr.param.methcallorbismark == "Methylation Call") {
    PCASamples(meth, screeplot = scree)
  } else if (visr.param.outputtype == "Differential Methylation Annotation" &&
             visr.param.methcallorbismark == "Methylation Call") {
    myDiff = calculateDiffMeth(meth)
    # get hyper methylated bases
    myDiff25p.hyper = getMethylDiff(myDiff,
                                    difference = visr.param.difference,
                                    qvalue = visr.param.qvalue,
                                    type = "hyper") # add difference, qvalue parameters
    # get hypo methylated bases
    myDiff25p.hypo = getMethylDiff(myDiff,
                                   difference = visr.param.difference,
                                   qvalue = visr.param.qvalue,
                                   type = "hypo")
    # get all differentially methylated bases
    myDiff25p = getMethylDiff(myDiff, difference = visr.param.difference, qvalue =
                                visr.param.qvalue)
    
    if (visr.param.overlap == "Exons/Introns/Promoters") {
      # read the gene BED file
      gene.obj = readTranscriptFeatures(system.file("extdata", "refseq.hg18.bed.txt",
                                                    package = "methylKit"))
      # annotate differentially methylated CpGs with promoter/exon/intron using annotation data
      annotateWithGeneParts(as(myDiff25p, "GRanges"), gene.obj)
      
      diffAnn = annotateWithGeneParts(as(myDiff25p, "GRanges"), gene.obj)
      # target.row is the row number in myDiff25p
      head(getAssociationWithTSS(diffAnn))
      
      plotTargetAnnotation(diffAnn, precedence = TRUE,
                           main = "differential methylation annotation")
    } else {
      # read the shores and flanking regions and name the flanks as shores and CpG islands as CpGi
      cpg.obj = readFeatureFlank(
        system.file("extdata", "cpgi.hg18.bed.txt",
                    package = "methylKit"),
        feature.flank.name = c("CpGi", "shores")
      )
      
      # convert methylDiff object to GRanges and annotate
      diffCpGann = annotateWithFeatureFlank(
        as(myDiff25p, "GRanges"),
        cpg.obj$CpGi,
        cpg.obj$shores,
        feature.name = "CpGi",
        flank.name = "shores"
      )
      
      plotTargetAnnotation(diffCpGann,
                           col = c("green", "gray", "white"),
                           main = "differential methylation annotation")
    }
  }
} else if (visr.param.onefileormulti == "One File") {
  if (visr.param.outputtype == "Methylation") {
    getMethylationStats(myobj, plot = TRUE, both.strands = visr.param.bothstrands)
  } else if (visr.param.outputtype == "Coverage") {
    getCoverageStats(myobj, plot = TRUE, both.strands = visr.param.bothstrands)
  }
}
