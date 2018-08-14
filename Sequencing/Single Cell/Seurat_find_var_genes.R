#Find variable genes
find_var_cond <- sprintf("visr.param.Find_Variable_Genes == T && %s",analysis_cond)
visr.app.category("Find Variable Genes", active.condition = find_var_cond)
visr.param("Normalization_scale_factor", min = 1, default = 10000, debugvalue = 10000,
           info = "Sets the scale factor for cell-level normalization")
visr.param("Mean_exp_low",label = "Minimum Mean Expression", min = 0, default = 0.0125,debugvalue = 0.0125,
           info = "Low cutoff on x-axis (mean expression) for identifying variable genes")
visr.param("Mean_exp_high", label = "Maximum Mean Expression", type = "double", min = 0, default = 3.0, debugvalue = 3.0, items = c("Inf"),item.labels = "No Upper Limit",
           info = "High cutoff on x-axis (mean expression) for identifying variable genes")
visr.param("Dispersion_low", label = "Minimum Dispersion", min = 0, default = 0.5, debugvalue = 0.5,
           info = "Low cutoff on y-axis (standard deviation) for identifying variable genes")
visr.param("Dispersion_high", label = "Maximum Dispersion", type = "double", min = 0, debugvalue = Inf, items = c("Inf"), item.labels = "No Upper Limit", default="Inf",
           info = "High cutoff on y-axis (standard deviation) for identifying variable genes")

plot_var_genes <- function(gbmData){
  par(mfrow = c(1,1), oma = c(0, 0, 1, 0))
  VariableGenePlot(gbmData,x.low.cutoff = visr.param.Mean_exp_low, x.high.cutoff = visr.param.Mean_exp_high, 
                   y.cutoff = visr.param.Dispersion_low, y.high.cutoff = visr.param.Dispersion_high)
  mtext("Variable Genes Selection", cex = 1.5)
}

find_variable_genes <- function(gbmData){
  if (is.null(gbmData@meta.data$percent.mito)){
    gbmData <- filter_Cells(gbmData)
  }
  # Normalization
  print(paste("Performing normalization"))
  normalization_scale_factor <- visr.param.Normalization_scale_factor
  raw_data <- gbmData@raw.data
  gbmData@raw.data <- gbmData@raw.data[rownames(gbmData@data),]
  gbmData <- NormalizeData(object = gbmData, normalization.method = "LogNormalize", scale.factor = normalization_scale_factor)
  gbmData@raw.data <- raw_data
  
  # Find variable genes
  print(paste("Finding variable genes"))
  x.low.cutoff <- visr.param.Mean_exp_low
  x.high.cutoff <- visr.param.Mean_exp_high
  y.cutoff <- visr.param.Dispersion_low
  y.high.cutoff <- visr.param.Dispersion_high
  print(x.low.cutoff)
  gbmData <- FindVariableGenes(object= gbmData, mean.function = ExpMean, dispersion.function = LogVMR, do.plot = F,
                               x.low.cutoff = x.low.cutoff, x.high.cutoff = x.high.cutoff, y.cutoff = y.cutoff, y.high.cutoff = y.high.cutoff)
  
  plot_var_genes(gbmData)
  switchPlotToScreen()
  plot_var_genes(gbmData)
  switchPlotToReport()
  
  print(paste(paste("Number of variable genes: ", toString(length(x = gbmData@var.genes)), sep = "")))
  if (length(gbmData@var.genes) == 0){
    visr.message("No variable genes found")
  }
  return(gbmData)
}

export_var_genes <- function(gbmData){
  if (length(gbmData@var.genes) == 0 || nrow(gbmData@hvg.info) == 0){return()}
  table <- gbmData@hvg.info
  table <- cbind(gene = rownames(table), table, is.variable=F)
  table[gbmData@var.genes,"is.variable"] <- T
  write.table(x = table, file = paste(output_folder,var_gene_output,sep = "/"), quote = F, 
              row.names = F, sep = "\t")
}
