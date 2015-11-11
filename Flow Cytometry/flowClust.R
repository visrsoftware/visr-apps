source("visrutils.R")
visr.biocLite("flowMerge")
visr.biocLite("flowClust")
visr.applyParameters()

d <- subset(input_table, select = input_columns)
tableMatrix <- data.matrix(d)

if(input_nu_est_string == "no estimation") {
	input_nu_est <- 0
} else if (input_nu_est_string == "estimation") {
	input_nu_est <- 1
} else {  #input_nu_est_string=="cluster-specific estimation"
	input_nu_est <- 2
}

if(input_trans_string == "no estimation") {
	input_trans <- 0
} else if (input_trans_string == "estimation") {
	input_trans <- 1
} else {  #input_trans_string=="cluster-specific estimation"
	input_trans <- 2
}

if(is_min_count_NULL) {
	input_min_count <- NULL
}

if(is_max_count_NULL) {
	input_max_count <- NULL
}

if(input_min == "NULL") {
	input_min <- NULL
} else {
	paste("c(",input_min,")")
	input_min <- eval(parse(text=paste("c(",input_min,")")))
}

if(input_max == "NULL") {
	input_max <- NULL
} else {
	paste("c(",input_max,")")
	input_max <- eval(parse(text=paste("c(",input_max,")")))
}

if(is_u_cutoff_NULL) {
	input_u_cutoff <- NULL
}

output_strucutre <- flowClust(tableMatrix,
     expName="Flow Experiment",
     K = 1:input_K,
     B = input_B,
     tol = input_tol,
     nu = input_nu,
     lambda = input_lambda,
     nu.est= input_nu_est,
     trans = input_trans,
     min.count = input_min_count,
     max.count = input_max_count,
     min = input_min,
     max = input_max,
     level = input_level,
     u.cutoff = input_u_cutoff,
     z.cutoff = input_z_cutoff,
     randomStart = input_randomStart, 
     B.init = input_B_init,
     tol.init = input_tol_init,
     seed = input_seed,
     criterion = input_criterion,
     control = NULL,   #An argument reserved for internal use.
     prior = NULL,
     usePrior = "no")

if(use_flowMerge) {
	# create flowFrame data for creating flowObj 
	dta <-  tableMatrix
	meta <- data.frame(name=dimnames(dta)[[2]], desc=paste('this is column',dimnames(dta)[[2]],'from your table (.csv)'))
	meta$range <- apply(apply(dta,2,range),2,diff)
	meta$minRange <- apply(dta,2,min)
	meta$maxRange <- apply(dta,2,max)
	ff <- new("flowFrame", exprs=dta, parameters=AnnotatedDataFrame(meta))


	#extract the BIC of each solution with the internal flowMerge function BIC
	#flowMerge:::BIC(flowClust.res)

	# extracted the max BIC solution, found for K=?
	maxBIC <- output_strucutre[[which.max(flowMerge:::BIC(output_strucutre))]];
	
	# reate a flowObj object from the flowClust result and the flowFrame data
	flowobj <- flowObj(maxBIC,ff);
	
	# run the merge function on the flowObj and extract then entropy of clustering from the merged results using the ENT function
	mergeResult<- merge(flowobj,metric="entropy");
	
	# extract the merged solution with number of clusters equal to the position of the changepoint found by fitPiecewiseLinreg.
	# Note: fitPiecewiseLinreg fits a piecewise linear regression to the entropy vs number of clusters and 
	# 		locates the position of the changepoint, if appropriate. Model selection is done via the BIC criterion.
	i <- fitPiecewiseLinreg(mergeResult);
	
	# the optimal merged solution based on the entropy criterion
	mergeopt <- mergeResult[[i]];
	
	output_column <- mergeopt@label
} else {
	output_column <- output_strucutre@.Data[[input_K]]@label
}


