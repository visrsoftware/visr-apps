
usePackage<-function (pkg) {if (!require(pkg, character.only = TRUE)) {install.packages(pkg, repos = "http://cran.us.r-project.org");require(pkg, character.only = TRUE)}}

usePackage("cluster")


visr.applyParameters()

{{
  output_pam <- pam(x = input_table[,input_columns],
                    k = input_k,
                    #diss = input_diss,
                    metric = input_metric,
                    #medoids = input_medoids, 
                    stand = input_stand, 
                    cluster.only = input_clusteronly,
                    #keep.diss = input_keepdiss,
                    #keep.data = input_keepdata, 
                    trace.lev = input_tracelev)
}}
print(capture.output(print(summary(output_pam))),collapse='\n')

output_clusterid <- output_pam$clustering

if (input_whichplot == "clusplot") {input_whichplot2 = 1}
if (input_whichplot == "silhouette plot"){input_whichplot2 = 2}

if(input_main == "" & input_whichplot == "clusplot") {input_main = "Clusplot"}
if(input_main == "" & input_whichplot == "silhouette") {input_main = "Silhouette plot"}

{{
  if (input_whichplot == "clusters"){
    #plot clusters
    plot (input_table[,input_columns], col = output_pam$clustering)
    #add the medoids to the plot
    points(output_pam$medoids, col = output_pam$clustering, pch = 4)
  }else{
  plot(x = output_pam, 
       ask = FALSE,
       which.plot = input_whichplot2,
       main = input_main,
       sub = input_sub,
       #xlim = input_xlim,
       #ylim = input_ylim,
       nmax.lab = input_nmaxlab,
       max.strlen = input_maxstrlen)
  }
}}

