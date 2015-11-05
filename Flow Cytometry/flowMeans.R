source("visrutils.R")
visr.library("flowMeans")
visr.applyParameters()

d <- subset(input_table, select = input_columns)
tableMatrix <- data.matrix(d)

# flowMeans uses is.na() to check NA value
if(is_maxN_NA) {
	input_maxN = NA
} 

if(is_numC_NA) {
	input_numC = NA
} 

if(is_MaxCovN_NA) {
	input_MaxCovN = NA
}

if(is_MaxKernN_NA) {
	input_MaxKernN = NA
} 

output_strucutre <- flowMeans(tableMatrix,
    MaxN = input_maxN,
    NumC = input_numC,
    iter.max = input_iter_max,
    nstart = input_nstart,
    Mahalanobis = input_do_Mahalanobis,
    Standardize = input_do_Standardize,
    Update = input_update,
    OrthagonalResiduals = input_do_OrthagonalResiduals,
    MaxCovN = input_MaxCovN,
    MaxKernN = input_MaxKernN,
    addNoise = input_do_addNoise)

output_column <- output_strucutre@Label


