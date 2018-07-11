analysis_cond <-  sprintf("(visr.param.Output_directory != '' && %s)", load_data_valid)
analysis_cond2 <-  sprintf("(visr.param.Output_directory != '' && %s)", load_integrated_valid)

visr.app.category("Analysis steps",info = "Select the analysis steps to run",
                  active.condition = paste(analysis_cond,analysis_cond2,sep = "||"))

visr.param("Find_Variable_Genes", default = F, debugvalue = F, active.condition = analysis_cond)
visr.param("Dim_Reduction", label = "Dimensionality Reduction", default = F, debugvalue = F)
visr.param("Cluster_Cells", default = F, debugvalue = F)
visr.param("Find_Marker_Genes",label = "Differential Expression Analysis", default = F, debug = F)
