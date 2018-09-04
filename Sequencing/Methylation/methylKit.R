source("visrutils.R")

visr.biocLite("methylKit")
visr.biocLite("genomation")

####################################################################################################
####################################################################################################
visr.app.start("Methyl Kit")
####################################################################################################
####################################################################################################

####################################################################################################
visr.category("Input Options")
####################################################################################################
visr.param("methcallorbismark", label = "Input Type", items = c("Methylation Call", "Bismark"))

param.max_input = 10L
visr.param("input_count", label = "Number of input files", default = 4L, min = 1L, max = param.max_input)

####################################################################################################
visr.category("Data Files")
####################################################################################################

for (input_index in seq(param.max_input)) {
  visr.param(paste0("file", input_index), label=paste0("File ", input_index), type = "filename", filename.mode = "load",
             info = "location of input file",
             active.condition = paste0("visr.param.input_count >= ", input_index))

  visr.param(paste0("sample.id", input_index), label = paste0("Title for Sample ", input_index),
             info = "The id of the sample for the specified input file. e.g. test1 or ctrl1",
             default = paste0("Sample ", input_index),
             active.condition = paste0("visr.param.input_count >= ", input_index))
}

####################################################################################################
visr.category("Sample Options")
####################################################################################################

visr.param("treatment", info = "list of 0 and 1 denoting which samples are control which samples are test",
           default = "1, 1, 0, 0")

visr.param("assembly", info = "a string that defines the genome assembly such as hg18, mm9. this is just a string for book keeping. It can be any string. Although, when using multiple files from the same assembly, this string should be consistent in each object.",
           default = "hg18")

visr.param("context", label="methylation context string", info = "Determines what type of methylation context will be read-in to the memory which can be immediately used for analysis.",
           items = c('CpG','CHG','CHH','none'))

####################################################################################################
visr.category("Output Options")
####################################################################################################

visr.param("save.folder", label = "Save Folder",
           type="filename", filename.mode = "dir",
           info = "The folder which will be used to save methylation call files. if not specified no methylation call file will be saved as a text file. The files saved can be read in less time using 'Methylation Call' input type option",
           active.condition = "visr.param.methcallorbismark == 'Bismark'")

visr.param(
  "outputtype",
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

####################################################################################################
visr.category("Methylation: Coverage Options")
####################################################################################################

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

####################################################################################################
visr.category("Differential Methylation Annotation Options")
####################################################################################################

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

####################################################################################################
visr.category("Filter By Coverage")
####################################################################################################

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

####################################################################################################
####################################################################################################
visr.app.end(printjson = TRUE, writefile = TRUE)
visr.applyParameters()
####################################################################################################
####################################################################################################

param.treatment = eval(parse(text=paste0("c(", visr.param.treatment, ")")))
param.location = visr.param.file1
param.sample.id = visr.param.sample.id1

param.save.folder = NULL
if (visr.param.save.folder != '')
  param.save.folder = visr.param.save.folder

if (visr.param.input_count > 1) {
  for (input_index in seq(2, visr.param.input_count)) {
    param.location = c(param.location, eval(parse(text=paste0("visr.param.file",input_index))))
    param.sample.id = c(param.sample.id, eval(parse(text=paste0("visr.param.sample.id",input_index))))
  }
  param.location = as.list(param.location)
  param.sample.id = as.list(param.sample.id)
}

# Reading files into list
if (visr.param.methcallorbismark == "Methylation Call") {
  # Read the files to a methylRawList object: myobj
  myobj = methRead(
    param.location,
    param.sample.id,
    assembly = visr.param.assembly,
    treatment = param.treatment,
    context = visr.param.context
  )
} else if (visr.param.methcallorbismark == "Bismark" &&
           visr.param.input_count == 1) {
  # Read the files to a methylRaw object: myobj
  myobj = processBismarkAln(
    location = param.location,
    sample.id = param.sample.id,
    assembly = visr.param.assembly,
    read.context = visr.param.context,
    save.folder = param.save.folder
  )
} else if (visr.param.methcallorbismark == "Bismark" &&
           visr.param.input_count > 1) {
  # Read the files to a methylRawList object: myobj
  myobj = processBismarkAln(
    location = param.location,
    sample.id = param.sample.id,
    assembly = visr.param.assembly,
    read.context = visr.param.context,
    save.folder = param.save.folder,
    treatment = param.treatment
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

if (visr.param.input_count > 1) {
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

if (visr.param.input_count > 1) {
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
} else if (visr.param.input_count == 1) {
  if (visr.param.outputtype == "Methylation") {
    getMethylationStats(myobj, plot = TRUE, both.strands = visr.param.bothstrands)
  } else if (visr.param.outputtype == "Coverage") {
    getCoverageStats(myobj, plot = TRUE, both.strands = visr.param.bothstrands)
  }
}
