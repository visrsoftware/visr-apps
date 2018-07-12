source("visrutils.R")

visr.app.start("ModEx", info = "A general purpose system for exploration of computational models")
visr.category("Input")
visr.param("directory", label = "Runs directory", type = "filename", filename.mode = "dir",
           debugvalue = "~/SFU/TFPlaygroundPSA/output/gauss_25_1")


DERIVATION_NONE <- "No derived output"
DERIVATION_FIRST_ROW <- "Take first row"
DERIVATION_COUNT_PER_CLASS <- "Aggregate: Count per class label"

visr.category("Derived output", active.condition = "visr.param.directory != ''")
visr.param("derivationMethod", label =  "Derivation Method",
           items = c(
             DERIVATION_NONE,
             DERIVATION_FIRST_ROW,
             "Aggregate: Number of class labels", # TODO
             "Aggregate: Mode of class labels", #TODO
             DERIVATION_COUNT_PER_CLASS,
             "Aggregate: Average", #TODO
             "Aggregate: Sum",#TODO
             "Aggregate: Median",#TODO
             "Aggregate: Min",#TODO
             "Aggregate: Max",#TODO
             "Dimensionality Reduction: PCA",#TODO
             "Dimensionality Reduction: MDS",#TODO
             "Dimensionality Reduction: tSNE"#TODO
           ),
           debugvalue = DERIVATION_NONE)
visr.param("tableForDerivatives", label = "Table to Use for Derivatives",
           items = c("", "quality_criteria"), #TODO: these items should be populated with proper table names in PSA java code
           active.condition = sprintf("visr.param.derivationMethod != '%s'", DERIVATION_NONE))

visr.param("columnsForDerivatives", label = "Column to Use for Derivatives",
           active.condition = sprintf("visr.param.derivationMethod != '%s' && visr.param.derivationMethod != '%s'", DERIVATION_NONE, DERIVATION_FIRST_ROW),
           debugvalue = "edgeR_is.de")

# visr.param("derivativePCA", label = "Calculate PCA on Derivatives", default = FALSE, # TODO
#           active.condition = paste0("visr.param.derivationMethod != '", DERIVATION_NONE, "'"))
#visr.param("derivativeMDA", label = "Calculate MDS on Derivatives", default = FALSE, # TODO
#           active.condition = paste0("visr.param.derivationMethod != '", DERIVATION_NONE, "'"))
#visr.param("derivativeTSNE", label = "Calculate tSNE on Derivatives", default = FALSE, # TODO
#           active.condition = paste0("visr.param.derivationMethod != '", DERIVATION_NONE, "'"))

visr.category("View Options", active.condition = "visr.param.directory != ''")
visr.param("launchExplorer", label = "Start Exploration", default = TRUE)
visr.param("output_showOutputDist", label = "Show Derived Output Filters", default = TRUE)
visr.param("impactSort", label = "Sort Parameters by impact", default = FALSE)

'
#TODO: remove, since this is handled through specifying derivationMethod
if (TRUE) {
  visr.param.output_mds <- FALSE
  visr.param.mdsColumnIndex <- 1
  visr.param.mdsMethod <- "manhattan"
  visr.param.summaryMDS <- "mds_"
} else {
  #TODO: remove
  visr.param("output_mds", label = "Calculate MDS for all runs", default = FALSE)
  visr.param("mdsColumnIndex", label = "Column index for MDS",
             info = "Index of output column to compute MDS and summary for. (1 = first column of output)",
             default = 1L, min = 1L,
             active.condition = "visr.param.output_mds == TRUE")
  #TODO: remove
  visr.param("mdsMethod", label = "distance metric for MDS",
             default = "manhattan",
             items = c(
               "euclidean",
               "maximum",
               "manhattan",
               "canberra",
               "binary",
               "minkowski",
               "hamming"
             ),
             active.condition = "visr.param.output_mds == TRUE"
             )
  visr.param("summaryMDS", label = "MDS column name",
             type = "output-multi-column",
             default = "mds_",
             active.condition = "visr.param.output_mds == TRUE"
  )

}
'

visr.app.end(printjson = TRUE, writefile = TRUE)
visr.applyParameters()

visr.param.recalc <- FALSE # not controlled by user for now
'
if (!visr.isGUI()) {
  ### FOR TESTING in RStudio ONLY - apply variables manually
  # Hamids DEBUG:
  .libPaths(c("~/VisR/RLibs", .libPaths()))
  #visr.param.directory <- "/Users/hyounesy/Downloads/spiral_25"
  visr.param.directory <- "~/SFU/visrseq-prototypes/Data/joy_clustering_all_300"
  #visr.param.directory <- "~/Research/Data/VisRseq/Runs/DESeq_EdgeR_LinSca"
  #visr.param.directory <- "~/Research/git/cradle-of-visrseq/Data/edgeR[500]_simulated_12492_up624_down625"
  #visr.param.directory <- "~/Research/Data/edgeR_2016-3-30(18.6)[500]_mydata_5spc"
  #visr.param.directory <- "~/Research/Data/PSA/edgeR_2016-4-14(7.51)[1000]_Jafar_SNwt_ESC_vs_SNcDKO_ESC"
  #visr.param.directory <- "~/SFU/visrseq-prototypes/Data/GSE49712/edgeR[4000]_GSE49712_SEQC_Count"

  visr.param.impactSort <- FALSE
  visr.param.recalc <- FALSE
  visr.param.output_mds <- FALSE
  visr.param.mdsColumnIndex <- 1
  visr.param.mdsMethod <- "manhattan"
  visr.param.summaryMDS <- "mds_"

  visr.param.tableForDerivatives <- "" # "quality_criteria" # "" #TODO: user input: string
  visr.param.derivationMethod <- DERIVATION_COUNT_PER_CLASS #DERIVATION_FIRST_ROW # "Count Levels"
}
'

path <- visr.param.directory

# initializing the file paths
pathInput    <- paste0(path, "/input.txt")
pathRunsInfo <- paste0(path, "/runsInfo.txt")
if (!file.exists(pathRunsInfo)) {
  pathIndexOld <- paste0(path, "/index.txt") # previous name for runsInfo
  if (file.exists(pathIndexOld)) {
    file.rename(pathIndexOld, pathRunsInfo)
  }
}

pathParamInfo           <- paste0(path, "/paramInfo.txt")
pathParamViewInfo       <- paste0(path, "/paramViewInfo.txt")
pathAllRunsMatrixRData  <- paste0(path, "/allRunsMatrix.RData")
pathAllRunsMatrixTxt    <- paste0(path, "/allRunsMatrix.txt")
pathDistanceMatrix      <- paste0(path, "/distanceMatrix.txt")
pathDistanceMatrixRData <- paste0(path, "/distanceMatrix.RData")
pathImpact              <- paste0(path, "/impact.txt")

allRunsMatrix <- NULL
distanceMatrix <- NULL


exportRunFilesFromAllRunsMatrix <- function() {
  numRuns = nrow(allRunsMatrix)
  for (i in seq_len(numRuns)) {
    visr.print(paste("saving",fileNames[i]))
    write.table(allRunsMatrix[i,], fileNames[i], row.names = FALSE, sep="\t")
  }
}

loadAllRunsMatrixIfNull<-function() {
  visr.logProgress("loading all runs matrix")
  if (is.null(allRunsMatrix)) {
    if (file.exists(pathAllRunsMatrixRData)) {
      # first try to load the binary version
      visr.print(paste("loading", pathAllRunsMatrixRData))
      load(file = pathAllRunsMatrixRData, .GlobalEnv)
    }
  }

  if (is.null(allRunsMatrix)) {
    # could not load the binary, try the txt one
    visr.print(paste("loading", pathAllRunsMatrixTxt))
    allRunsMatrix <<- read.table(pathAllRunsMatrixTxt, header=FALSE, sep="\t", check.names = F)
    # write the binary
    visr.print(paste("saving", pathAllRunsMatrixRData))
    save(allRunsMatrix, file = pathAllRunsMatrixRData)
  }
}

loadDistanceMatrixIfNull<-function() {
  visr.logProgress("loading distance matrix")
  if (is.null(distanceMatrix)) {
    if (file.exists(pathDistanceMatrixRData)) {
      # first try to load the binary version
      visr.print(paste("loading", pathDistanceMatrixRData))
      load(file = pathDistanceMatrixRData, .GlobalEnv)
    }
  }

  if (is.null(distanceMatrix)) {
    # could not load the binary, try the txt one
    visr.print(paste("loading", pathDistanceMatrix))
    distanceMatrix <<- read.table(pathDistanceMatrix, sep="\t", check.names = F)
    # write the binary
    visr.print(paste("saving", pathDistanceMatrixRData))
    save(distanceMatrix, file = pathDistanceMatrixRData)
  }
}


calculateMode=function(x){
  names(tail(sort(table(x)),1))
}

# whether to store the minimum distance from a runs with one parameter values to nearset run with a different parameter value
# It will be one column for each parameter
storeMinDistColumns <- FALSE

###############################################################################
# Setting paths and loading needed files
###############################################################################

visr.logProgress("Loading the runs files")
visr.assert_file_exists(pathInput, "input")
visr.assert_file_exists(pathParamInfo, "paramInfo")
visr.assert_file_exists(pathRunsInfo, "runsInfo")

inputTable <- read.table(pathInput, header=TRUE, sep="\t", check.names = F)
paramInfo  <- read.table(pathParamInfo, header=TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = F)
runsInfoTable <- read.table(pathRunsInfo, header=TRUE, sep = "\t", check.names = F)
if ("succeded" %in% colnames(runsInfoTable)) {
  # only take the runs that succeeded.
  # TODO: better way of handling the failed runs. e.g. to show which parameter combinations were problematic.
  runsInfoTable <- runsInfoTable[which(runsInfoTable$succeded == "true"),]
}

if (visr.param.tableForDerivatives == "") {
  fileNames <- paste0(path, "/runs/", (runsInfoTable[,"ID"]), ".txt")
} else {
  fileNames <- paste0(path, "/runs/", visr.param.tableForDerivatives, "_", (runsInfoTable[,"ID"]), ".txt")
}

numRuns <- length(fileNames)

for (name in names(inputTable)) {
  if (any(is.na(inputTable[,name])) || any(grep("NA", inputTable[,name])) ||  any(grep("NaN", inputTable[,name]))) {
    inputTable[,name] <- NULL
  }
}
visr.print("[DONE] loading the runs files")

###############################################################################
# Generate paramViewInfo if doesn't exist
###############################################################################
if (!file.exists(pathParamViewInfo)) {
  visr.print("creating paramViewInfo")
  viewParams <- c()
  viewTypes <- c()
  for (i in 1: nrow(paramInfo)) {
    #if( !grepl("output", paramInfo$type[i]) )
    { # skip the output parameters
      if (!paramInfo$label[i] %in% c("ID", "imagePath")) {
        viewParams <- c(viewParams, paramInfo$label[i])
        viewTypes <- c(viewTypes, if (grepl("output", paramInfo$type[i])) "Output" else "Parameter")
      }
    }
  }

  paramViewInfo = data.frame(c1 = c(0:(length(viewParams)-1)), c2=viewParams, c3=viewTypes, c4=-1)

  write.table(paramViewInfo, file=pathParamViewInfo, row.names = FALSE, col.names = FALSE, sep = "\t", quote=FALSE)
}

###############################################################################
# Calculate derived output for all Runs (dOutput)
###############################################################################

if (visr.param.derivationMethod != DERIVATION_NONE) {
  visr.print(paste("Calculating derivative results:", visr.param.derivationMethod))

  #if (visr.param.recalc || (!file.exists(pathAllRunsMatrixRData) && !file.exists(pathAllRunsMatrixTxt)))

  if (numRuns > 0) {
    for (i in 1:numRuns) {
      visr.logProgress(paste("reading",fileNames[i]))
      tt <- read.table(fileNames[i], header=TRUE,sep="\t", check.names = F)
      if (visr.param.derivationMethod == DERIVATION_FIRST_ROW) {
        if (i == 1) {
          dOutput <- tt[1,] # take the first row
        } else {
          dOutput <- rbind(dOutput, tt[1,]) # append the first row
        }
      } else if (visr.param.derivationMethod == DERIVATION_COUNT_PER_CLASS) {
        if (i == 1) {
          dOutput <- matrix(nrow=numRuns, ncol=0) # matrix of derivations
        }
        for (derivCol in visr.param.columnsForDerivatives) {
          tt_col <- tt[,derivCol]
          tt_table <- table(tt_col)
          for (n1 in names(tt_table)) {
            colName <- paste0(derivCol, "(", n1, ")") # e.g. "is.de(-1)"
            if (!colName %in% colnames(dOutput)) {
              # seeing the factor level for the first time
              newColnames <- c(colnames(dOutput), colName)
              dOutput <- cbind(dOutput, rep(0, numRuns))
              colnames(dOutput) <- newColnames
            }
            dOutput[i, colName] <- tt_table[n1]
          }
        }
        '
        if (i == 1) {
          allRunsMatrix <- matrix(nrow = numRuns, ncol=nrow(tt))
        }
        allRunsMatrix[i,] <- tt[, visr.param.mdsColumnIndex]
        '
      } else {
        visr.message(msg = sprintf("The derivation method '%s' cannot be computed", visr.param.derivationMethod), type = "error")
      }
      #else if (visr.param.derivationMethod == DERIVATION_MDS) { #TODO
      #}
    }

    # post-process
    #if (visr.param.derivationMethod == DERIVATION_FIRST_ROW) {
      # take the numeric columns and transpose
    allRunsMatrix <- as.matrix(dOutput[, which(apply(dOutput, 2, is.numeric))])
    #}
  }

  save(allRunsMatrix, file = pathAllRunsMatrixRData)

  if (exists("dOutput") && ncol(dOutput) > 0) {
    for (c1 in colnames(dOutput)) {
      runsInfoTable[[c1]] <- dOutput[,c1]
    }
    # writing the table, in case of a crash, etc...
    write.table(runsInfoTable, file=pathRunsInfo, row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE)

    # Add missing meta info to the parameter information file
    for (c1 in colnames(dOutput)) {
      if(! c1 %in% paramInfo$label) {
        paramInfo <- rbind(paramInfo, c(c1, c1, "output_generated"))
      }
    }
  }
  visr.print("[DONE] Calculating derivative results")
}

###############################################################################
# Calculate Distance matrix (distance between runs)
###############################################################################
if (visr.param.derivationMethod != DERIVATION_NONE) {
  loadAllRunsMatrixIfNull()
  visr.logProgress("computing the distance matrix. may take some time...")
  distanceMeasure <- "euclidean" # , "maximum", "manhattan", "canberra", "binary" or "minkowski"
  distanceMatrix <- as.matrix(dist(allRunsMatrix, method = distanceMeasure)) #  distances between the rows
  # save distanceMatrix to a file
  # write.table(distanceMatrix, file=pathDistanceMatrix, col.names = F, row.names = F, sep = "\t")
  visr.print(paste("saving", pathDistanceMatrixRData))
  save(distanceMatrix, file = pathDistanceMatrixRData)

  '
  if ((!file.exists(pathDistanceMatrixRData) && !file.exists(pathDistanceMatrix)) || # whether the distance matrix was already computed and saved
      visr.param.recalc) # forcing the recalculation
  {
    loadAllRunsMatrixIfNull()
    visr.logProgress("computing the distance matrix. may take some time...")
    # recompute the distance matrix
    if (visr.param.mdsMethod == "hamming")
    { # hamming ditance not supported in dist() function. need to calculate manually
      n <- nrow(allRunsMatrix)
      distanceMatrix <- matrix(nrow=n, ncol=n)
      for(i in seq_len(n)) {
        visr.print(paste("computing hamming dist for row ", i))
        for(j in seq(i, n)) {
          distanceMatrix[j, i] <- distanceMatrix[i, j] <- sum(allRunsMatrix[i,] != allRunsMatrix[j,])
        }
      }
    } else {
      distanceMatrix <- as.matrix(dist(allRunsMatrix, method = visr.param.mdsMethod)) #  distances between the rows
      # save distanceMatrix to a file
      # write.table(distanceMatrix, file=pathDistanceMatrix, col.names = F, row.names = F, sep = "\t")
      visr.print(paste("saving", pathDistanceMatrixRData))
      save(distanceMatrix, file = pathDistanceMatrixRData)
    }
    visr.print("[DONE] recompute the distance matrix.")
  } else {
    # load the existing distance matrix from file.
    # delay loading this until needed (performance optimization):
    #distanceMatrix <- read.table(pathDistanceMatrix, sep="\t", check.names = F)
  }
  '
  ###############################################################################
  # perform MDS on the distance matrix and add it to the runsInfo table
  ###############################################################################
  '
  if (visr.param.output_mds) {
    mdsDims <- 2
    mdsColNames <- paste(visr.param.summaryMDS, "coord" , 1:mdsDims, sep = "")
    c1ColName   <- paste(visr.param.summaryMDS, "coord1", sep = "")
    c2ColName <- paste(visr.param.summaryMDS, "coord2", sep = "")
    if ( !all(mdsColNames %in% colnames(runsInfoTable)) || visr.param.recalc == TRUE) {
      visr.logProgress("computing the MDS. may take some time...")
      visr.print("perform MDS on the distance matrix")
      loadDistanceMatrixIfNull()
      visr.logProgress("computing the MDS. may take some time...")
      fit <- cmdscale(distanceMatrix, eig = TRUE, k = mdsDims) # k is the number of dim
      visr.param.summaryMDS <- fit$points[, c(1:mdsDims)] # If subscript out of bounds error: means most likely, that all columns of the input data were identical and therefore removed
      # colnames(visr.param.summaryMDS)<-c("coord1","coord2")

      for(i in 1:mdsDims) {
        runsInfoTable[[mdsColNames[i]]] <- fit$points[,i]

        # Add missing mds_coordinate meta info to the parameter information file
        if(! (mdsColNames[i] %in% paramInfo$label) ) {
          paramInfo <- rbind(paramInfo, c(mdsColNames[i], mdsColNames[i], "output_generated"))
        }
      }
      visr.print("[Done] perform MDS on the distance matrix")
    }
  }
  '
}
###############################################################################
# Calculate Parameters' significance/impact score
###############################################################################

if(visr.param.recalc ||
   !file.exists(pathImpact)) # whether the impact file exists already
{
  if (visr.param.derivationMethod != DERIVATION_NONE) {
    loadAllRunsMatrixIfNull()
    loadDistanceMatrixIfNull()
  }
  visr.logProgress("Calculating Parameters' significance/impact score")
  visr.print("Calculate Parameters' significance/impact score")
  # find columns of the runsInfo table, that are input parameters with more than 1 value
  numlevels<-function(x) {
    length(levels(as.factor(x)))
  }

  coolParameters <- match(names(which(lapply(runsInfoTable, numlevels) > 1)), paramInfo$label)
  coolParameters <- coolParameters[!is.na(coolParameters)]
  coolParameters <- coolParameters[order(coolParameters)]
  inputParams <- c()
  # visr.print( coolParameters)
  # TODO make this more R-like
  for (i in 1: length(coolParameters)) {
    #if( !grepl("output", paramInfo$type[coolParameters[i]]) )
    { # skip the output parameters
      inputParams <- c(inputParams, coolParameters[i])
    }
  }

  # data frame containing the significance/impact values per parameters
  sigFrame <- data.frame(row.names = paramInfo$label[inputParams]) #error row names contain missing values

  if (length(inputParams) > 0) {
    for (vp in 1: length(inputParams)) {
      # iterate through input parameters to calculate impact score
      #paramName<-colnames(runsInfoTable)[inputParams[vp]]
      visr.logProgress(paste("Calculating Parameters' significance/impact score:", vp, "of", length(inputParams)))
      paramName<-paramInfo$label[inputParams[vp]]
      visr.print(paste("processing", paramName))
      paramValues <- c()

      if (storeMinDistColumns) {
        minDistColName <- paste("mindist(",paramName,")", sep = "")
      }

      # Binning for numerical parameters that have more than 10 values
      binned <- FALSE
      paramValues<-levels(as.factor(runsInfoTable[, paramName]))
      numParamValues<-length(paramValues)
      MAX_NUMERIC_PARAM_CATEGORIES = 10
      NUMERIC_PARAM_BIN_COUNTS = 10
      if (numParamValues > MAX_NUMERIC_PARAM_CATEGORIES &&
          (grepl("double", paramInfo$type[inputParams[vp]]) ||
           grepl("int", paramInfo$type[inputParams[vp]]))) {
        binned <- TRUE
        bins <- tapply(runsInfoTable[, paramName], cut(runsInfoTable[, paramName], NUMERIC_PARAM_BIN_COUNTS))
        paramValues <- levels(as.factor(bins))
        numParamValues <- length(paramValues)
      }

      isOutputParam = grepl("output", paramInfo$type[inputParams[vp]])

      if (!isOutputParam && visr.param.derivationMethod != DERIVATION_NONE) {
        sigV <- c()
        if (FALSE && numParamValues > 1) {
          breakdown <- data.frame(row.names = paramValues)

          paramsDist<-matrix(0,nrow=numParamValues,ncol=numParamValues)

          # now compare each pair of parameter values and find the minimum distance between rows
          for (v1 in 1:(numParamValues)) {
            if (binned) {
              v1Rows <- which(bins == paramValues[v1]) # index of rows with their parameter value in paramValues[v1] bin
            } else {
              v1Rows <- which(runsInfoTable[, paramName]==paramValues[v1]) # index of rows with their parameter value == paramValues[v1]
            }

            if (visr.param.derivationMethod == DERIVATION_COUNT_PER_CLASS) {
              # Calculate the mode for each value of each input parameter
              vRuns = allRunsMatrix[v1Rows,]
              if(is.null(dim(vRuns))) {
                modes = data.frame(sapply(vRuns, calculateMode))
              } else {
                modes = data.frame(apply(vRuns, 2, calculateMode))
              }
              inputTable[,paste(paramName,paramValues[v1],sep = "-")] <- modes
            }


            if (v1 >= numParamValues) {
              break# Only calculate impatct for numParamValues -1
            }

            for (v2 in (v1+1) : numParamValues) {
              visr.print(paste("comparing", paramValues[v1], "with", paramValues[v2]))
              if (binned) {
                v2Rows <- which(bins == paramValues[v2]) # index of rows with their parameter value in paramValues[v2] bin
              } else {
                v2Rows <- which(runsInfoTable[, paramName]==paramValues[v2]) # index of rows with their parameter value == paramValues[v2]
              }
              mm <- distanceMatrix[v1Rows, v2Rows]

              if (length(mm) > 1 && length(dim(mm) > 1)) {
                v1MinDist <- apply(mm,1,min)
                v2MinDist <- apply(mm,2,min)
                minDist <- mean(c(mean(v1MinDist), mean(v2MinDist)))
                if (storeMinDistColumns) {
                  runsInfoTable[v1Rows, minDistColName] <- v1MinDist
                  runsInfoTable[v2Rows, minDistColName] <- v2MinDist
                }
              } else {
                minDist <- min(mm)
              }
              paramsDist[v1,v2] <- paramsDist[v2,v1] <- minDist;
              visr.print(paste("average min dist = ", minDist))
              sigV <- c(sigV, minDist)

              # Breakdown
              breakdown[v1, paramValues[v2]] <- minDist
              breakdown[paramValues[v2], v1] <- minDist
            }
            breakdown[v1, v1] <- 0
          }

          # parameterImpact = breakdown[]
          breakdown = apply(breakdown, 1, mean)
          write.table(breakdown, file=paste(path, "/breakdown-", paramName, ".txt", sep = ""), row.names = TRUE, col.names = NA, sep = "\t", quote=TRUE)
        }
        sigFrame$sig[vp] <- mean(sigV)
      } else {
        # can't calculate actual significance without full output
        sigFrame$sig[vp] <- 0
      }
      sigFrame$index[vp] <- inputParams[vp]
    }
    visr.print("[DONE] Calculate Parameters' significance/impact score")
  }

  # Order and save the sigFrame
  sigFrame <- sigFrame[with(sigFrame, order(-sig)), ]
  write.table(as.matrix(sigFrame), file = pathImpact, col.names = NA, sep = "\t", quote=TRUE)
}

# visr.message("after impact")

# Generate the frequency / stability of the assigned clusters to the genes

ENABLE_GENERATE_FREQUENCY_TABLE = FALSE
if (ENABLE_GENERATE_FREQUENCY_TABLE) {
  source("frequency.R")
  pathFrequency           <- paste(path, "/frequency.txt",sep = "")
  pathStableGenes         <- paste(path, "/stableGenes.txt",sep = "")
  if(!file.exists(pathFrequency) || !file.exists(pathStableGenes) || visr.param.recalc) {
    visr.print("Generate the frequency / stability of the assigned clusters to the genes")

    loadAllRunsMatrixIfNull()
    loadDistanceMatrixIfNull()

    visr.logProgress("Generating the frequency / stability tables")
    f <- freq(allRunsMatrix)
    f$name <- inputTable$name # Add the names of the genes to the table for easier identification
    f <- f[with(f, order(-freq)), ] # Order them descending for frequency
    top <- f[(f$sym != "0") & (f$freq > 0.95),] # Assable the top (most stable) genes which are not unregulate
    write.table(f, file=pathFrequency, row.names = TRUE, col.names = NA, sep = "\t", quote=TRUE)
    write.table(top, file=pathStableGenes, row.names = TRUE, col.names = NA, sep = "\t", quote=TRUE)

    visr.print("[DONE] Generated the frequency / stability of the assigned clusters to the genes")
  } else {
    f <- read.table(pathFrequency, header=TRUE,sep="\t", check.names = F)
    top <- read.table(pathStableGenes, header=TRUE,sep="\t", check.names = F)
  }
}

# Save the possibly changed files
visr.print("Saving the (possibly) changed files ...")

write.table(runsInfoTable, file=pathRunsInfo,     row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE)
write.table(paramInfo,  file=pathParamInfo, row.names = FALSE, col.names = TRUE, sep = "\t", quote=TRUE)
visr.print("[DONE] Saving the (possibly) changed files")
