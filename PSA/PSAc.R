source("visrutils.R")
source("frequency.R")
visr.applyParameters()

### FOR TESTING ONLY - apply variables manually
# source("C:/Users/Joseph/Documents/visrseq/visr/srcR/frequency.R")
# # source("B:/Diverses/Dropbox/Uni/BA/VisRSeq/visr/srcRfrequency.R")
# # visr.param.directory = "B:/Diverses/Dropbox/Uni/BA/VisRSeq/output/LinSca_Lin-Pos-Neg_edgeR_115-9-22(6.43)[550]"
# # visr.param.directory = "C:/Users/Joseph/Documents/visrseq/output/LinSca_Lin-Pos-Neg_edgeR_115-9-22(6.43)[550]"
# visr.param.directory = "C:/Users/Joseph/Documents/visrseq/output/k-means_2016-1-26(18.14)[6]_iris-sub"
# visr.param.mdsColumnIndex = 1
# visr.param.summaryUpDown = "sum_"
# visr.param.recalc = T
# visr.param.mdsMethod = "manhattan"
# visr.output.summaryMDS = "mds_"
###

### FUNCTIONS - should be put in an other file!
numlevels<-function(x) {length(levels(as.factor(x)))}
loadMDS<-function(path) {read.table(paste(path, "/allRunsMatrix.txt", sep=""), header=FALSE, sep="\t", check.names = F)}

paths = c(visr.param.directory, visr.param.directory2)
paths = unique(paths[paths!=""])
pathIndices = paste(paths, "/runsInfo.txt" ,sep="")
pathParamInfos = paste(paths, "/paramInfo.txt" ,sep="")
inputTables = c()
paramInfos = c()
indexTables = c()
runsList = list()
	
for (i in (1:length(paths))){
	inputTables[i] = read.table(paste(paths[i], "/input.txt", sep = ""), header=TRUE,sep="\t", check.names = F)
	paramInfos[i] = read.table(pathParamInfos[i], header=TRUE,sep="\t", stringsAsFactors = FALSE, check.names = F)
	indexTables[i] = read.table(pathIndices[i], header=TRUE, sep="\t", check.names = F)
	runsList[i] = paste(paths[i], "/runs/", (indexTable[,"ID"]), ".txt",sep="")
	
}
numRuns <- sum(sapply(runsList, nrow))

mdsMatrix = NULL

if(length(paths) <= 1){
	return
}

if (numRuns > 0) {
	for (j in 1:length(runsList)){
		for (i in 1:length(runsList[[j]])) {
			tt <- read.table(runsList[[j]][i], header=TRUE,sep="\t", check.names = F)
			
			if (i == 1) {
				mdsMatrix <- matrix(nrow = numRuns, ncol=nrow(tt))
			}
			mdsMatrix[i,] <- tt[, visr.param.mdsColumnIndex]
		}
	}
	
}
write.table(mdsMatrix, file=paste(path,"/allRunsMatrix.txt",sep=""),col.names = F, row.names = F, sep = "\t")
indexTable[[upColName]] = upDown[,1]
indexTable[[downColName]] = upDown[,2]

#if(!file.exists(paste(path,"/distanceMatrix.txt",sep="")) || visr.param.recalc == TRUE) {
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

# perform MDS on the distance matrix (d)
mdsDims = 2
mdsColNames = paste(visr.output.summaryMDS, "coord" , 1:mdsDims, sep = "")
c1ColName 	=	paste(visr.output.summaryMDS, "coord1", sep = "")
c2ColName = 	paste(visr.output.summaryMDS, "coord2", sep = "")
if ( !all(mdsColNames %in% colnames(indexTable)) || visr.param.recalc == TRUE) {
	fit <- cmdscale(d, eig = TRUE, k = mdsDims) # k is the number of dim
	visr.output.summaryMDS <- fit$points[, c(1:mdsDims)] # If subscript out of bounds error: means most likely, that all columns of the input data were identical and therefore removed
	# colnames(visr.output.summaryMDS)<-c("coord1","coord2")
	
	for(i in 1:mdsDims) {
		indexTable[[mdsColNames[i]]] = fit$points[,i]
		
		# Add missing mds_coordinate meta info to the parameter information file
		if(! (mdsColNames[i] %in% paramInfo$label) ) {
			paramInfo = rbind(paramInfo, c(mdsColNames[i], mdsColNames[i], "output_generated"))
		}
	}
}
