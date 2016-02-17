#setwd("~/Research/svn/Projects/GeneCalc/visr/srcR/")

# source("visrutils.R")
# source("frequency.R")
# visr.applyParameters()

# Hamid's DEBUG:
#visr.param.directory <- "/Users/hyounesy/Research/Data/VisRseq/Runs/DESeq_EdgeR_LinSca"

### FOR TESTING ONLY - apply variables manually
source("C:/Users/Joseph/Documents/visrseq/visr/srcR/frequency.R")
# source("B:/Diverses/Dropbox/Uni/BA/VisRSeq/visr/srcRfrequency.R")
# visr.param.directory <- "B:/Diverses/Dropbox/Uni/BA/VisRSeq/output/LinSca_Lin-Pos-Neg_edgeR_115-9-22(6.43)[550]"
visr.param.directory <- "C:/Users/Joseph/Documents/visrseq/output/LinSca_Lin-Pos-Neg_edgeR_115-9-22(6.43)[550]"
# visr.param.directory <- "C:/Users/Joseph/Documents/visrseq/output/test-method"
# visr.param.directory <- "C:/Users/Joseph/Documents/visrseq/output/k-means_2016-1-26(18.14)[6]_iris-sub"
visr.param.mdsColumnIndex <- 1
visr.param.summaryUpDown <- "sum_"
visr.param.recalc <- F
visr.param.mdsMethod <- "manhattan"
visr.output.summaryMDS <- "mds_"
###

### FUNCTIONS - should be put in an other file!
loadMDS<-function(path) {
  read.table(paste(path, "/allRunsMatrix.txt", sep=""), header=FALSE, sep="\t", check.names = F)
}

calculateMode=function(x){
	names(tail(sort(table(x)),1))
}

# whether to store the minimum distance from a runs with one parameter values to nearset run with a different parameter value
# It will be one column for each parameter
storeMinDistColumns <- FALSE

###############################################################################
# Setting paths and loading crucial index files
###############################################################################
path <- visr.param.directory
pathIndex <- paste(path, "/index.txt", sep="")
pathParamInfo <- paste(path, "/paramInfo.txt", sep="")
pathInput = paste(path, "/input.txt", sep = "")

inputTable <- read.table(pathInput, header=TRUE,sep="\t", check.names = F)
paramInfo <- read.table(pathParamInfo, header=TRUE,sep = "\t", stringsAsFactors = FALSE, check.names = F)
indexTable <- read.table(pathIndex, header=TRUE, sep = "\t", check.names = F)

fileNames <- paste(path, "/runs/", (indexTable[,"ID"]), ".txt", sep="")
numRuns <- length(fileNames)

for (name in names(inputTable)) {
	if (any(is.na(inputTable[,name])) || any(grep("NA", inputTable[,name])) ||  any(grep("NaN", inputTable[,name]))) {
		inputTable[,name] <- NULL
	}
}


###############################################################################
# Calculate sum of up and down regulated genes for all Runs
###############################################################################
upColName   <- paste(visr.param.summaryUpDown, "up", sep = "")
downColName <- paste(visr.param.summaryUpDown, "down", sep = "")
upDown <- matrix(nrow=numRuns, ncol=2)
colnames(upDown) <- c(upColName, downColName)

# !(upColName %in% colnames(indexTable)) ||    # check if index table has the number of up-regulated genes for each run
# 	!(downColName %in% colnames(indexTable)) ||  # check if index table has the number of down-regulated genes for each run
if (
	  visr.param.recalc || # force recalculating
    !file.exists(paste(path,"/allRunsMatrix.txt",sep="")) )
{
  # Recreate the matrix containing the up/down regulation state for all runs
  # TODO slow (not R enough?)
  if (numRuns > 0) {
    for (i in 1:numRuns) {
      tt <- read.table(fileNames[i], header=TRUE,sep="\t", check.names = F)

      upDown[i,1] <- sum(tt[,visr.param.mdsColumnIndex] > 0.1)
      upDown[i,2] <- sum(tt[,visr.param.mdsColumnIndex] < -0.1)

      if (i == 1) {
      	mdsMatrix <- matrix(nrow = numRuns, ncol=nrow(tt))
      }
      mdsMatrix[i,] <- tt[, visr.param.mdsColumnIndex]
    }
  }
  write.table(mdsMatrix, file=paste(path,"/allRunsMatrix.txt",sep=""),col.names = F, row.names = F, sep = "\t")
  indexTable[[upColName]] <- upDown[,1]
  indexTable[[downColName]] <- upDown[,2]

} else {
	mdsMatrix <- loadMDS(path)
# 	upDown[,1] <- indexTable[[upColName]]
# 	upDown[,2] <- indexTable[[downColName]]
}

# Add missing up/down regulated meta info to the parameter information file
if(! upColName %in% paramInfo$label) {
	paramInfo <- rbind(paramInfo, c(upColName, upColName, "output_generated"))
}
if(! downColName %in% paramInfo$label) {
	paramInfo <- rbind(paramInfo, c(downColName, downColName, "output_generated"))
}


#plot(upDown)

###############################################################################
# Calculate Distance matrix (distance between runs)
###############################################################################
if (!file.exists(paste(path,"/distanceMatrix.txt",sep="")) || # whether the distance matrix was already computed and saved
    visr.param.recalc) # forcing the recalculation
{
  # recompute the distance matrix
	if (visr.param.mdsMethod == "hamming")
  { # hamming ditance not supported in dist() function. need to calculate manually
		n <- nrow(mdsMatrix)
		distanceMatrix <- matrix(nrow=n, ncol=n)
		for(i in seq_len(n))
			for(j in seq(i, n))
			  distanceMatrix[j, i] <- distanceMatrix[i, j] <- sum(mdsMatrix[i,] != mdsMatrix[j,])
	} else {
	  distanceMatrix <- as.matrix(dist(mdsMatrix, method = visr.param.mdsMethod)) #  distances between the rows
		# save distanceMatrix to a file
		write.table(distanceMatrix, file=paste(path,"/distanceMatrix.txt",sep=""),col.names = F, row.names = F, sep = "\t")
	}
} else {
  # load the existing distance matrix from file.
  distanceMatrix <- read.table(paste(path, "/distanceMatrix.txt" ,sep=""), sep="\t", check.names = F)
}


###############################################################################
# perform MDS on the distance matrix and add it to the index table
###############################################################################
mdsDims <- 2
mdsColNames <- paste(visr.output.summaryMDS, "coord" , 1:mdsDims, sep = "")
c1ColName   <- paste(visr.output.summaryMDS, "coord1", sep = "")
c2ColName <- paste(visr.output.summaryMDS, "coord2", sep = "")
if ( !all(mdsColNames %in% colnames(indexTable)) || visr.param.recalc == TRUE) {
	fit <- cmdscale(distanceMatrix, eig = TRUE, k = mdsDims) # k is the number of dim
	visr.output.summaryMDS <- fit$points[, c(1:mdsDims)] # If subscript out of bounds error: means most likely, that all columns of the input data were identical and therefore removed
	# colnames(visr.output.summaryMDS)<-c("coord1","coord2")

	for(i in 1:mdsDims) {
		indexTable[[mdsColNames[i]]] <- fit$points[,i]

		# Add missing mds_coordinate meta info to the parameter information file
		if(! (mdsColNames[i] %in% paramInfo$label) ) {
			paramInfo <- rbind(paramInfo, c(mdsColNames[i], mdsColNames[i], "output_generated"))
		}
	}
}

###############################################################################
# Calculate Parameters' significance/impact score
###############################################################################
if(visr.param.recalc ||
   !file.exists(paste(path, "/impact.txt", sep=""))) # whether the impact file exists already
{
  # find columns of the index table, that are input parameters with more than 1 value
  numlevels<-function(x) { length(levels(as.factor(x))) }
	coolParameters <- match(names(which(lapply(indexTable, numlevels) > 1)), paramInfo$label)
	coolParameters <- coolParameters[!is.na(coolParameters)]
	inputParams <- c()
	# print( coolParameters)
	# TODO make this more R-like
	for (i in 1: length(coolParameters)) {
		if( !grepl("output", paramInfo$type[coolParameters[i]]) ) { # skip the output parameters
			inputParams <- c(inputParams, coolParameters[i])
		}
	}

  # data frame containing the significance/impact values per parameters
	sigFrame <- data.frame(row.names = paramInfo$label[inputParams]) #error row names contain missing values

	if(length(inputParams) > 0) {
		for (vp in 1: length(inputParams))
    { # iterate through input parameters to calculate impact score
			#paramName<-colnames(indexTable)[inputParams[vp]]
			paramName<-paramInfo$label[inputParams[vp]]
			print(paste("processing", paramName))
			paramValues <- c()

			if (storeMinDistColumns) {
			  minDistColName <- paste("mindist(",paramName,")", sep = "")
			}

			# Binning for numerical parameters that have more than 10 values
			binned <- FALSE
			paramValues<-levels(as.factor(indexTable[, paramName]))
			numParamValues<-length(paramValues)
			if ( numParamValues > 10 &&
           (grepl("double", paramInfo$type[inputParams[vp]]) ||
            grepl("int", paramInfo$type[inputParams[vp]])) ) {
				binned <- TRUE
				bins <- tapply(indexTable[, paramName], cut(indexTable[, paramName], 10))
				paramValues <- levels(as.factor(bins))
				numParamValues <- length(paramValues)
			}

			sigV <- c()
			if (numParamValues > 1) {
				breakdown <- data.frame(row.names = paramValues)

				paramsDist<-matrix(0,nrow=numParamValues,ncol=numParamValues)

        		# now compare each pair of parameter values and find the minimum distance between rows
				for (v1 in 1:(numParamValues)) {
					if(binned) {
						v1Rows <- which(bins == paramValues[v1]) # index of rows with their parameter value in paramValues[v1] bin
					}else {
						v1Rows <- which(indexTable[, paramName]==paramValues[v1]) # index of rows with their parameter value == paramValues[v1]
					}
					
					# Calculate the mode for each value of each input parameter
					vRuns = mdsMatrix[v1Rows,]
					modes = data.frame(apply(vRuns, 2, calculateMode))
					inputTable[,paste(paramName,paramValues[v1],sep = "-")] = modes
					if(v1 >= numParamValues) break; # Only calculate impatct for numParamValues -1
					
					for (v2 in (v1+1) : numParamValues) {
						print(paste("comparing", paramValues[v1], "to", paramValues[v2]))
						if(binned) {
						  v2Rows <- which(bins == paramValues[v2]) # index of rows with their parameter value in paramValues[v2] bin
						} else {
						  v2Rows <- which(indexTable[, paramName]==paramValues[v2]) # index of rows with their parameter value == paramValues[v2]
						}
						mm <- distanceMatrix[v1Rows, v2Rows]

						if (length(mm) > 1 && length(dim(mm) > 1)) {
						  v1MinDist <- apply(mm,1,min)
						  v2MinDist <- apply(mm,2,min)
						  minDist <- mean(c(mean(v1MinDist), mean(v2MinDist)))

              			if (storeMinDistColumns) {
  						  indexTable[v1Rows, minDistColName] <- v1MinDist
  						  indexTable[v2Rows, minDistColName] <- v2MinDist
              			}

						} else {
							minDist <- min(mm)
						}
						paramsDist[v1,v2] <- paramsDist[v2,v1] <- minDist;
						print(paste("average min dist = ", minDist))
						sigV <- c(sigV, minDist)

						# Breakdown
						breakdown[v1, paramValues[v2]] <- minDist
						breakdown[paramValues[v2], v1] <- minDist
					}
					breakdown[v1, v1] <- 0
				}

				# parameterImpact = breakdown[]
				breakdown = apply(breakdown, 1, mean)
				write.table(breakdown, file=paste(path, "/breakdown-", paramName, ".txt", sep = ""), row.names = TRUE, col.names = NA, sep = "\t")
			}
			sigFrame$sig[vp] <- mean(sigV)
			sigFrame$index[vp] <- inputParams[vp]

			# Normalize stuff
			#sigFrame$
		}
	}

	# Order and save the sigFrame
	sigFrame <- sigFrame[with(sigFrame, order(-sig)), ]
	write.table(as.matrix(sigFrame), file=paste(path,"/impact.txt",sep=""), col.names = NA, sep = "\t")
}

# visr.message("after impact")

# Generate the frequency / stability of the assigned clusters to the genes
if(!file.exists(paste(path,"/frequency.txt",sep = "")) || !file.exists(paste(path,"/stableGenes.txt",sep = "")) || visr.param.recalc == TRUE) {
	if(is.null(mdsMatrix)) loadMDS(path)
	f <- freq(mdsMatrix)
	f$name <- inputTable$name # Add the names of the genes to the table for easier identification
	f <- f[with(f, order(-freq)), ] # Order them descending for frequency
	top <- f[(f$sym != "0") & (f$freq > 0.95),] # Assable the top (most stable) genes which are not unregulate
	write.table(f, file=paste(path,"/frequency.txt",sep = ""), row.names = TRUE, col.names = NA, sep = "\t")
	write.table(top, file=paste(path,"/stableGenes.txt",sep = ""), row.names = TRUE, col.names = NA, sep = "\t")
} else {
	f <- read.table(paste(path, "/frequency.txt", sep = ""), header=TRUE,sep="\t", check.names = F)
	top <- read.table(paste(path, "/stableGenes.txt", sep = ""), header=TRUE,sep="\t", check.names = F)
}

# Save the possibly changed files
write.table(inputTable, file=pathInput,     row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(indexTable, file=pathIndex,     row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(paramInfo,  file=pathParamInfo, row.names = FALSE, col.names = TRUE, sep = "\t")
