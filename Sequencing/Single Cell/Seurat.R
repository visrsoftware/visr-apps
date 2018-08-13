source("visrutils.R")
curr_dir <- "Sequencing/Single Cell"

visr.app.start("Seurat", debugdata=mtcars, input.type = "none")
source(sprintf("%s/Utils.R",curr_dir))
source(sprintf("%s/Seurat_IO.R",curr_dir))
source(sprintf("%s/Seurat_filter.R",curr_dir))
source(sprintf("%s/Seurat_steps.R",curr_dir))
source(sprintf("%s/Seurat_find_var_genes.R",curr_dir))
source(sprintf("%s/Seurat_dim_reduction.R",curr_dir))
source(sprintf("%s/Seurat_cluster_cells.R",curr_dir))
source(sprintf("%s/Seurat_DE.R",curr_dir))
source(sprintf("%s/Seurat_additional.R",curr_dir))
visr.app.end(printjson=T, writefile=T)
visr.applyParameters()

# Check Input -------------------------------------------------------------
if (visr.param.workflow == "single"){
  if (visr.param.Import_method == "load_raw"){
    if (!dir.exists(sprintf("%s/outs",visr.param.Path_to_outs))){
      visr.message(paste("'",visr.param.Path_to_outs,"'"," is not a valid path to Cell Ranger pipeline output directory. It should contain another directory named \"outs\""))
    }
  } else {
    if (!file.exists(visr.param.Path_to_seurat_object)){visr.message(sprintf("'%s' is not a valid path to seurat object", visr.param.Path_to_seurat_object))}
  }
}else{
  if (visr.param.Import_method2 == "load_one"){
    assert_that(file.exists(visr.param.seurat_integrated), msg = sprintf("'%s' is not a valid path to seurat object",visr.param.seurat_integrated))
  }else{
    assert_that(file.exists(visr.param.seurat_obj_1), msg = sprintf("'%s' is not a valid path to seurat object 1",visr.param.seurat_obj_1))
    assert_that(file.exists(visr.param.seurat_obj_2), msg = sprintf("'%s' is not a valid path to seurat object 2",visr.param.seurat_obj_2))
  }
}

if (!dir.exists(visr.param.Output_directory)){visr.message("Path to output directory not found")}

assert_that(visr.param.max_nGenes > visr.param.min_nGenes,msg = "The high cutoff of the number of genes must be greather than the low cutoff" )
assert_that(visr.param.Mean_exp_high > visr.param.Mean_exp_low, msg = "The high cutoff of the mean expression level must be greather than the low cutoff")
assert_that(visr.param.Dispersion_high > visr.param.Dispersion_low, msg =  "The high cutoff of dispersion must be greather than the low cutoff" )
# Load Packages -----------------------------------------------------------

visr.library("Seurat")
visr.library("dplyr")
visr.library("reshape2")

# Analysis result file names ----------------------------------------------

plot_output <- "plots.pdf"
cell_info_output <- "Cells.tsv"
DE_output <- "Differential_Expression_Analysis_result.tsv"
var_gene_output <- "Variable_Genes.tsv"
selected_gene_output <- "Selected_Marker_Genes.tsv"

# Create PDF --------------------------------------------------------------

# save plots to this location
if (visr.param.create_subdir){
  subfolder <- Sys.time()
  subfolder <- gsub(" ","_",subfolder)
  subfolder <- gsub(":","-",subfolder)
  subfolder <- gsub("-","",subfolder)
  output_folder <- paste(visr.param.Output_directory,subfolder,sep = "/")
} else {
  output_folder <- visr.param.Output_directory
}
dir.create(output_folder)

# shut off exisiting device
# visr.dev <- dev.cur()
# if (exists("visr_app_pdf_dev") && !is.null(visr_app_pdf_dev)){dev.off(which = visr_app_pdf_dev)}
# pdf(file = paste(output_folder,plot_output,sep="/"))
# visr_app_pdf_dev <- dev.cur()

startReport(output_dir = output_folder)

# first page
par(oma = c(0, 0, 12, 0))
plot.new()
mtext(text = "Seurat Ouput Plots", outer = T, cex = 2,line = 1)

# Load Data ---------------------------------------------------------------
min_fraction_of_cells <- visr.param.min_fraction_cells
min_number_of_genes <- visr.param.min_nGenes
gbmData <- load_data()

if (visr.param.workflow == "integrated"){
  groups <- levels(as.factor(gbmData@meta.data$group))
}

# data summary
# plotSeuratSummary(gbmData)

# Additional checkpoint ---------------------------------------------------

if (visr.param.workflow == "integrated" && visr.param.Import_method2 == "load_one"){
  assert_that(!is.null(gbmData@meta.data$group),msg = "Need to merge datasets first")
}

if (visr.param.workflow == "integrated" && visr.param.Import_method2 == "load_two"){
  assert_that(visr.param.dataset1_name != visr.param.dataset2_name, msg = "Label of the two datasets must be different.")
}

if (length(visr.param.export_gene_list) > 0){
  valid_gene_probes <- check_gene_list(gbmData,visr.param.export_gene_list)
}

# Analysis ----------------------------------------------------------------

if (visr.param.workflow == "single"){
  # Find variable genes
  if (visr.param.Find_Variable_Genes){
    gbmData <- find_variable_genes(gbmData)
  }
  
  if (visr.param.Dim_Reduction){
    #PCA
    if (visr.param.Run_PCA){
      gbmData <- run_PCA(gbmData)
    }
    
    # Jackstraw
    if (visr.param.jackstraw){
      gbmData <- run_jackstraw(gbmData)
    }
    
    # PC heatmap
    if (visr.param.PC_heatmap){
      gbmData <- draw_PCHeatmap(gbmData)
    }
    #elbow plot
    if (visr.param.elbow){
      gbmData <- draw_elbow(gbmData)
    }
    
    if (visr.param.Run_tSNE){
      gbmData <- run_tSNE(gbmData)
    }
  }
  
  # Cluster cells
  if (visr.param.Cluster_Cells){
    gbmData <- cluster_cells(gbmData, "pca")
  }
  
  # DE analysis
  if (visr.param.Find_Marker_Genes){
   # top_genes_n <- min(visr.param.Top_gene_n,nrow(gbmData@raw.data))
   table <- diff_exp(gbmData)
   #plot_DE(gbmData, table)
  }
  
  # output
  if (visr.param.rename_id){
    gbmData <- rename_cluster(gbmData, "pca")
  }
  
  if (nchar(visr.param.export_gene_list) > 0){
    visualize_gene(gbmData)
  }
  
}else{
  if (visr.param.Dim_Reduction){
    if (visr.param.Run_CCA){
      gbmData <- run_CCA(gbmData) 
    }
    if (visr.param.bicor){
      gbmData <- plot_bicor(gbmData)
    }
    if (visr.param.align_CCA){
      gbmData <- align_subspace(gbmData)
    }
    if (visr.param.Run_tSNE2){
      gbmData <- run_tSNE2(gbmData) 
    }
  }
  
  if (visr.param.Cluster_Cells){
    gbmData <- cluster_cells(gbmData, "cca.aligned")
  }
  
  
  if (visr.param.Find_Marker_Genes){
    
    if (visr.param.choose_de == "conserved"){
      table <- diff_exp_conserved(gbmData)
    }else{
      table <- diff_exp_across(gbmData)
    }
    # plot_DE(gbmData, table)
  }
  
  # output
  if (visr.param.rename_id){
    gbmData <- rename_cluster(gbmData, "cca.aligned")
  }
  
  if (nchar(visr.param.export_gene_list) > 0){
    visualize_gene(gbmData)
  }
}

# Save Object -------------------------------------------------------------
export_var_genes(gbmData)
export_cells(gbmData)

plotSeuratSummary(gbmData)

output_parameters <- visr.getParams()
sink(file = paste(output_folder,"parameters.txt",sep = "/"))
print(output_parameters)
sink()

# close current device
# dev.off(which=visr_app_pdf_dev)
# rm(visr_app_pdf_dev)

finishReport()

browseURL(paste(output_folder,plot_output,sep = "/"))

# Reduce file size of saved object
if (!is.null(gbmData@calc.params$RunTSNE)){
  gbmData@calc.params$RunTSNE$... <- NULL
}

if (visr.param.save_obj){
  # save object
  print(paste("Saving object"))
  saveRDS(gbmData, file = paste(output_folder,"Seurat.rds",sep = "/"))
}

rm(gbmData)
