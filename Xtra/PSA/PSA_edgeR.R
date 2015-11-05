source("visrutils.R")
source("frequency.R")

visr.applyParameters()


### FOR TESTING ONLY - apply variables manually
# #source("C:/Users/Joseph/Documents/visrseq/visr/srcR/frequency.R")
# source("B:/Diverses/Dropbox/Uni/BA/VisRSeq/visr/srcRfrequency.R")
# visr.param.directory = "B:/Diverses/Dropbox/Uni/BA/VisRSeq/output/LinSca_Lin-Pos-Neg_edgeR_115-9-22(6.43)[550]"
# #visr.param.directory = "C:/Users/Joseph/Documents/visrseq/output/LinSca_Lin-Pos-Neg_edgeR_115-9-22(6.43)[550]"
# visr.param.mdsColumnIndex = 1
# visr.param.summaryUpDown = "sum_"
# visr.param.recalc = F
# visr.param.mdsMethod = "manhattan"
# visr.output.summaryMDS = "mds_"
###

### FUNCTIONS - should be put in an other file!
numlevels<-function(x) {length(levels(as.factor(x)))}
loadMDS<-function(path) {read.table(paste(path, "/allRunsMatrix.txt", sep=""), header=FALSE, sep="\t", check.names = F)}

# Setting paths and loading crucial files
path = visr.param.directory
pathIndex = paste(path, "/index.txt" ,sep="")
pathParamInfo = paste(path, "/paramInfo.txt" ,sep="")


inputTable = read.table(paste(path, "/input.txt", sep = ""), header=TRUE,sep="\t", check.names = F)
paramInfo = read.table(pathParamInfo, header=TRUE,sep="\t", stringsAsFactors = FALSE, check.names = F)
indexTable <- read.table(pathIndex, header=TRUE, sep="\t", check.names = F)

fileNames <- paste(path, "/runs/", (indexTable[,"ID"]), ".txt",sep="")
numRuns <- length(fileNames)

mdsMatrix = NULL


# Calculate sum of up and down regulated genes
upColName 	=	paste(visr.param.summaryUpDown, "up", sep = "")
downColName = 	paste(visr.param.summaryUpDown, "down", sep = "")
upDown <- matrix(nrow=numRuns, ncol=2)
colnames(upDown) <- c(upColName, downColName)
if ( !(upColName %in% colnames(indexTable)) || !(downColName %in% colnames(indexTable)) || 
	 visr.param.recalc == TRUE || !file.exists(paste(path,"/allRunsMatrix.txt",sep="")) )
{
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
  indexTable[[upColName]] = upDown[,1]
  indexTable[[downColName]] = upDown[,2]
  
} else {
#   numrows <- nrow(input_table) # number of runs
#   numcols <- 8000 #hypotetical number of genes
#   mdsMatrix <- as.integer((sample.int(7,numrows * numcols,TRUE) - 4L) / 3)
#   dim(mdsMatrix) <- c(numrows,numcols)
	upDown[,1] = indexTable[[upColName]]
	upDown[,2] = indexTable[[downColName]]
}

# Add missing up/down regulated meta info to the parameter information file
if(! upColName %in% paramInfo$label) {
	paramInfo = rbind(paramInfo, c(upColName, upColName, "output_generated"))
}
if(! downColName %in% paramInfo$label) {
	paramInfo = rbind(paramInfo, c(downColName, downColName, "output_generated"))
}



plot(upDown)

#visr.output.summaryOneColumn <- mdsMatrix

#if (visr.param.mdsMethod == "PCA") {
#  output_princomp <- princomp(x = mdsMatrix, cor = TRUE, scores = TRUE)
#  ouput_scores <- output_princomp$scores
#  visr.output.summaryMDS <- data.matrix(ouput_scores[,c(1,2)])
#  colnames(visr.output.summaryMDS)<-c("pc1","pc2")
#} else 

if(!file.exists(paste(path,"/distanceMatrix.txt",sep="")) || visr.param.recalc == TRUE) {
	mdsMatrix = loadMDS(path)
	
	if (visr.param.mdsMethod == "hamming") {
		n <- nrow(mdsMatrix)
		d <- matrix(nrow=n, ncol=n)
		for(i in seq_len(n))
			for(j in seq(i, n))
				d[j, i] <- d[i, j] <- sum(mdsMatrix[i,] != mdsMatrix[j,])
	} else {
		d <- as.matrix(dist(mdsMatrix, method = visr.param.mdsMethod)) #  distances between the rows
		# save d to a file
		write.table(d, file=paste(path,"/distanceMatrix.txt",sep=""),col.names = F, row.names = F, sep = "\t")
	}
} else {
	d = read.table(paste(path, "/distanceMatrix.txt" ,sep=""), sep="\t", check.names = F)
}



# perform MDS on the distance matrix (d)
mdsDims = 4
mdsColNames = paste(visr.output.summaryMDS, "coord" , 1:mdsDims, sep = "")
c1ColName 	=	paste(visr.output.summaryMDS, "coord1", sep = "")
c2ColName = 	paste(visr.output.summaryMDS, "coord2", sep = "")
if ( !all(mdsColNames %in% colnames(indexTable)) || visr.param.recalc == TRUE) {
	fit <- cmdscale(d, eig = TRUE, k = mdsDims) # k is the number of dim
	visr.output.summaryMDS <- fit$points[,c(1:mdsDims)]
	# colnames(visr.output.summaryMDS)<-c("coord1","coord2")
	
	for(i in 1:mdsDims) {
		indexTable[[mdsColNames[i]]] = fit$points[,i]
		
		# Add missing mds_coordinate meta info to the parameter information file
		if(! (mdsColNames[i] %in% paramInfo$label) ) {
			paramInfo = rbind(paramInfo, c(mdsColNames[i], mdsColNames[i], "output_generated"))
		}
	}
}

# Calculate Parameter Significance
if(visr.param.recalc == TRUE || !file.exists(paste(path, "/impact.txt", sep=""))) {
	coolParameters<-match(names(which(lapply(indexTable, numlevels) > 1)), paramInfo$label)
	coolParameters = coolParameters[!is.na(coolParameters)]
	inputParams = c()
	#print( coolParameters)
	# TODO make this more r-like
	for (i in 1: length(coolParameters)){
		if( !grepl("output", paramInfo$type[coolParameters[i]]) ) {
			inputParams = c(inputParams, coolParameters[i])
		}
	}
	#print( paramInfo$label[inputParams])
	sigFrame = data.frame(row.names = paramInfo$label[inputParams]) #error row names contain missing values
	
	if(length(inputParams) > 0) {
		for (vp in 1: length(inputParams)) {
			#paramName<-colnames(indexTable)[inputParams[vp]]
			paramName<-paramInfo$label[inputParams[vp]]
			print(paste("processing", paramName))
			paramValues = c()
			
			# Binning
			if ( grepl("double", paramInfo$type[inputParams[vp]]) || grepl("int", paramInfo$type[inputParams[vp]]) ) {
				binned = TRUE
				bins = tapply(indexTable[, paramName], cut(indexTable[, paramName], 10))
				paramValues = levels(as.factor(bins))
			}else { 
				binned = FALSE
				paramValues<-levels(as.factor(indexTable[, paramName]))
			}
			
			numParamValues<-length(paramValues)
			
			sigV = c()
			if (numParamValues > 1) {
				breakdown = data.frame(row.names = paramValues)
				
				paramsDist<-matrix(0,nrow=numParamValues,ncol=numParamValues)
				for (v1 in 1:(numParamValues-1)) {
					for (v2 in (v1+1) : numParamValues) {
						print(paste("comparing", paramValues[v1], "to", paramValues[v2]))
							
						if(binned) {
							mm<-d[which(bins==paramValues[v1]), which(bins==paramValues[v2])] #error? - incorrect number of dimensions
						} else {
							mm<-d[which(indexTable[, paramName]==paramValues[v1]), which(indexTable[, paramName]==paramValues[v2])] #error? - incorrect number of dimensions
						}
							
						if(length(mm) > 1) {
							minDist <- mean(c(mean(apply(mm,1,min)), mean(apply(mm,2,min)))) 
						} else {
							minDist = mm
						}
						paramsDist[v1,v2] = paramsDist[v2,v1] = minDist;
						print(paste("average min dist = ", minDist))
						sigV = c(sigV, minDist)
						
						# Breakdown
						breakdown[v1,paramValues[v2]] = minDist
					}
				}
				
				write.table(breakdown, file=paste(path, "/breakdown-", paramName, ".txt", sep = ""), row.names = TRUE, col.names = NA, sep = "\t")
			}
			sigFrame$sig[vp] = mean(sigV)
			sigFrame$index[vp] = inputParams[vp]
			
			# Normalize stuff
			#sigFrame$ 
		}
	}
	
	# Order and save the sigFrame
	sigFrame = sigFrame[with(sigFrame, order(-sig)), ]
	write.table(as.matrix(sigFrame), file=paste(path,"/impact.txt",sep=""),col.names = NA, sep = "\t")
}

# Generate the frequency / stability of the assigned clusters to the genes
if(!file.exists(paste(path,"/frequency.txt",sep = "")) || !file.exists(paste(path,"/stableGenes.txt",sep = "")) || visr.param.recalc == TRUE) {
	if(is.null(mdsMatrix)) loadMDS(path)
	f = freq(mdsMatrix)
	f$name = inputTable$name # Add the names of the genes to the table for easier identification
	f = f[with(f, order(-freq)), ] # Order them descending for frequency
	top = f[(f$sym != "0") & (f$freq > 0.95),] # Assable the top (most stable) genes which are not unregulate
	write.table(f, file=paste(path,"/frequency.txt",sep = ""), row.names = TRUE, col.names = NA, sep = "\t")
	write.table(top, file=paste(path,"/stableGenes.txt",sep = ""), row.names = TRUE, col.names = NA, sep = "\t")
} else {
	f = read.table(paste(path, "/frequency.txt", sep = ""), header=TRUE,sep="\t", check.names = F)
	top = read.table(paste(path, "/stableGenes.txt", sep = ""), header=TRUE,sep="\t", check.names = F)
}

# Save the possibly changed files
write.table(indexTable, file=pathIndex, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(paramInfo, file=pathParamInfo, row.names = FALSE, col.names = TRUE, sep = "\t")
