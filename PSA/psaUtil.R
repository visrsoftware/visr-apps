visr.message<-function(text, type=c("error","warning")) {
  
}

loadDistanceMatrixIfNull<-function(directory){
  # if (is.null(distanceMatrix)) {
    visr.print("loading the distanceMatrix.txt")
    distanceMatrix <<- read.table(paste(directory, "/distanceMatrix.txt" ,sep=""), sep="\t", check.names = F)
    return(distanceMatrix)
  # }
}

loadAllRunsMatrix<-function(directory) {
  # if (is.null(allRunsMatrix)) {
    visr.print("loading the allRunsMatrix.txt")
    allRunsMatrix <<- read.table(paste(directory, "/allRunsMatrix.txt", sep=""), header=FALSE, sep="\t", check.names = F)
    return(allRunsMatrix)
  # }
}