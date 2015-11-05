#setwd('/Users/hyounesy/Research/svn/Projects/GeneCalc/visr//srcR')
#options(repos='http://cran.rstudio.com/')
#.libPaths(c("/Users/hyounesy/VisRseq/RLibs", .libPaths()))
source("visrutils.R")

visr.biocLite("S4Vectors")
visr.biocLite("org.Mm.eg.db")
visr.biocLite("graphite")
visr.library("car")
visr.library("Rcpp")
visr.library("igraph")
visr.biocLite("RBGL")
visr.library("gRbase")
visr.biocLite("qpgraph")
visr.biocLite("KEGGgraph")
visr.library("corpcor")
visr.library("rrcov")
visr.libraryURL("timeClip","http://romualdi.bio.unipd.it/wp-uploads/2014/11/timeClip_0.99.3.tar.gz")
visr.library("graph")
visr.library("nlme")

# loading these due to the local modifications required to run the packages
source("lib/timeClip/R/getExpression.R")
source("lib/timeClip/R/removeNArows.R")
source("lib/timeClip/R/prunePaths.R")
source("lib/timeClip/R/clipperBackbone.R")
source("lib/timeClip/R/timeClip.R")
source("lib/timeClip/R/plotTimeInCytoscape.R")
source("lib/timeClip/R/getJunctionTreePaths.R")
source("lib/timeClip/R/cliqueTimeCourseTest.R")
source("lib/timeClip/R/covariance.R")
source("lib/timeClip/R/dag-graph.R")
source("lib/timeClip/R/deleteEdge.R")
source("lib/timeClip/R/getGraphEntryGenes.R")
source("lib/timeClip/R/pathwayTimeCourseTest.R")
source("lib/timeClip/R/removeSelfLoops.R")
source("lib/timeClip/R/myTimeClipSpaced.R",chdir=TRUE)

######## Paramters ############

visr.applyParameters()

###### Imporint data ##########

col_num <- length(visr.param.columns)

if (visr.param.pathwaydb == "kegg") {
  pathwaynames <- c(visr.param.kegg.pathways)
} else if (visr.param.pathwaydb == "reactome") {
  pathwaynames <- c(visr.param.reactome.pathways)
}
Path_num <- length(pathwaynames)

mydata <- input_table[,visr.param.columns]
colnames(mydata) <- c(1:length(visr.param.columns))
rownames(mydata) <- as.vector(input_table[,visr.param.columnTranscriptID])

# awk '{print $1"\t"$2+1-1"\t"$3+1-1"\t"$4+1-1"\t"$5+1-1"\t"$6+1-1"\t"$7+1-1"\t"$8+1-1}' Endothelium.txt > x
#as.numeric(as.character(mydata))
#apply(mydata, 2, is.numeric)
#allpaths = list()

allpathsmatrix <- matrix(nrow=Path_num, ncol=col_num + 1) # additional column is for p-value

#library(S4Vectors)
#library(org.Mm.eg.db)
#library(graph)
#library(graphite)
#library(timeClip)
#library(nlme)

i <- 0
pathwaydb <- pathways(visr.param.species, visr.param.pathwaydb)

for (pathwayname in pathwaynames ) {

  wasSuccessfull = T
  i <- i + 1
  if (is.null(pathwaydb[[pathwayname]])) { # skip the pathway
    wasSuccessfull = F
    #print(paste("Skipping NULL pathway: ", pathwayname))
  }
  else { # generate the pathway
    #print(paste("Processing pathway: ", pathwayname))
    paths <-convertIdentifiers(pathwaydb[pathwayname], visr.param.genedb)
    pathgraph <-convertIdentifiers(pathwaydb[[pathwayname]], visr.param.genedb)
    genes <- unlist(lapply(paths, nodes))
    g <- pathwayGraph(pathgraph)

    #genes <- as.vector(genes[grep('NM_', as.vector(genes),'r')])
    #g <- subGraph(genes, g)
    #length(unique(genes))
    #as.vector(genes)
    #write(as.vector(genes), file = "191.xls")

    ######### Converting to appropriate format ########

    timedata <- as.matrix(mydata)
    times <- as.numeric(colnames(timedata))
    genessubset <- intersect(as.vector(genes), rownames(timedata))
    #grep('NM_', as.vector(genes),'r')
    data.subset <- timedata[genessubset,]
    typeof(data.subset)
    if (length(genessubset)==0)
      wasSuccessfull = F
  }

  if(wasSuccessfull) {
    ######### Runnign TimeClip stage 1 ########
    #library(graph)
    g_temp <- subGraph(genessubset, g)
    if (class(try(pathwayTimeCourse(data.subset, times, g_temp, npc=1, eqids=c(1:col_num)),silent = TRUE)) == "try-error")
    {
      allpathsmatrix[i, 1 : col_num] <- matrix(0, col_num)
      allpathsmatrix[i, col_num + 1] <- c(1)
    } else {
      pathwayAnalysis <- pathwayTimeCourse(data.subset, times, g_temp, npc=1, eqids=c(1 : col_num))
      pathwayAnalysis$alpha
      allpathsmatrix[i,] <- as.matrix(c(pathwayAnalysis$pc, pathwayAnalysis$alpha))
      ##########################################

      #allpaths[[i]] <- c(pathwayAnalysis$pc,pathwayAnalysis$alpha)
    }
  } else {
    allpathsmatrix[i, 1 : col_num] <- matrix(0,col_num)
    allpathsmatrix[i, col_num + 1] <- c(1)
  }
}

#remove(allpathsmatrix)
#allpathsmatrix <- matrix(nrow=240, ncol=10)
#for (i in 1:5 ) {

#allpathsmatrix[i,] <- as.matrix(allpaths[[i]])
#}

rownames(allpathsmatrix) <- pathwaynames #names(pathwaydb[1:Path_num])
colnames(allpathsmatrix) <- c(visr.param.columns, "p-value")

visr.output.pathwaymatrix<-allpathsmatrix
#library(Hmisc)
#write.table(allpathsmatrix, file = "~/test_timeclip.txt",row.names = TRUE,col.names=FALSE,sep = "\t")
