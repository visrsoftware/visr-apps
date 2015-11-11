source("visrutils.R")
visr.biocLite("SamSPECTRAL")
visr.applyParameters()

d <- subset(input_table, select = input_columns)
tableMatrix <- data.matrix(d)

# SamSPECTRAL uses =="NA" to check NA value
if(is_number_of_clusters_NA) {
	input_number_of_clusters = "NA"
} 

if(is_k_for_kmeans_NA) {
	input_k_for_kmeans = "NA"
}

if(is_min_eigenvalue_NA) {
	input_min_eigenvalue = "NA"
} 

# convert input string vector to actual vector
if(input_scale == "") {
	input_scale <- rep(1,dim(tableMatrix)[2])
} else {
	paste("c(",input_scale,")")
	input_scale <- eval(parse(text=paste("c(",input_scale,")")))
}



output_column <- SamSPECTRAL(tableMatrix,
     normal.sigma = input_normal_sigma,
     separation.factor = input_separation_factor,
     number.of.clusters = input_number_of_clusters,
     scale = input_scale,
     #talk = TRUE,
     precision = input_precision,
     #eigenvalues.num =NA,
     #return_only.labels=TRUE,
     do.sampling = input_do_sampling,
     beta = input_beta,
     stabilizer = input_stabilizer,
     k.for_kmeans = input_k_for_kmeans,
     maximum.number.of.clusters = 30,
     m = input_m,
     minimum.eigenvalue = input_min_eigenvalue,
     previous.result = NULL, 
     #replace.inf.with.extremum=TRUE,
     minimum.degree = input_minimum_degree)


