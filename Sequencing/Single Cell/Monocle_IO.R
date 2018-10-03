################################################################
visr.category(label="Input")
################################################################

INPUT_TYPE_10X = "10X single cell dataset"
INPUT_TYPE_MONOCLE_OBJECT = "Load Monocle App Object (.RDS)"
INPUT_TYPE_COUNT_MATRIX_TXT = "Count matrix (.txt)"
#TODO:
INPUT_TYPE_SPARSE_MATRIX = "Sparse matrix"
visr.param("input_type", label = "Choose Import Method",
           info = "Specify the format of your input data. Whether it is a 10X cell ranger dataset or existing monocle object",
           items = c(INPUT_TYPE_10X, INPUT_TYPE_MONOCLE_OBJECT, INPUT_TYPE_COUNT_MATRIX_TXT, INPUT_TYPE_SPARSE_MATRIX),
           debugvalue = INPUT_TYPE_MONOCLE_OBJECT)

visr.param("data_dir_10x",
           label="10X dataset directory",
           info="10X Cell Ranger dataset directory. It should contain a sub-directory named 'outs'",
           type="filename", filename.mode = "dir",
           active.condition = sprintf("visr.param.input_type == '%s'", INPUT_TYPE_10X),
           debugvalue= "~/SFU/Datasets/SingleCell/pbmc3k/") # Peripheral Blood Mononuclear Cells (PBMCs) from a healthy donor

visr.param("path_to_monocle_object", type="filename",filename.mode = "load",
           info = "Path to the existing monocle object (monocle_app_object.RDS)",
           active.condition = sprintf("visr.param.input_type == '%s'", INPUT_TYPE_MONOCLE_OBJECT),
           debugvalue= "~/SFU/Datasets/SingleCell/output_monocle/monocle_app_object.RDS")

visr.param("expression_matrix", type = "filename",
           info = "Numeric matrix of expression values, where rows are genes, and columns are cells",
           active.condition = sprintf("(visr.param.input_type == '%s') || (visr.param.input_type == '%s')", INPUT_TYPE_COUNT_MATRIX_TXT, INPUT_TYPE_SPARSE_MATRIX))

visr.param("sample_sheet", type = "filename",
           info = "Pheno data where rows are cells, and columns are cell attributes (such as cell type, culture condition, day captured, etc.). It should have same number of rows as the number of columns in expression value matrix.",
           active.condition = sprintf("(visr.param.input_type == '%s') || (visr.param.input_type == '%s')", INPUT_TYPE_COUNT_MATRIX_TXT, INPUT_TYPE_SPARSE_MATRIX))

visr.param("gene_annotation", type = "filename",
           info = "Feature data where rows are features (e.g. genes), and columns are gene attributes (such as biotype, gc content, etc.). It should have same number of rows as the number of rows in expression value matrix.",
           active.condition = sprintf("(visr.param.input_type == '%s') || (visr.param.input_type == '%s')", INPUT_TYPE_COUNT_MATRIX_TXT, INPUT_TYPE_SPARSE_MATRIX))

DATA_TYPE_UMI = "UMIs, Raw transcript counts"
DATA_TYPE_FPKM = "FPKM, TPM"
DATA_TYPE_LOG_FPKM = "log FPKM/TPMs, Ct values from SC qPCR"
visr.param("data_type",
           label = "Data type",
           info = "Specify the appropriate data type. Monocle works well with both count-based measures (e.g. UMIs) and relative expression data.",
           items=c(DATA_TYPE_UMI, DATA_TYPE_FPKM, DATA_TYPE_LOG_FPKM),
           active.condition = "(visr.param.data_dir_10x != '') || (visr.param.expression_matrix != '' && visr.param.sample_sheet != '' && visr.param.gene_annotation != '')")

################################################################
visr.category("Output",
              active.condition = "(visr.param.data_dir_10x != '') || (visr.param.path_to_monocle_object != '') || (visr.param.expression_matrix != '' && visr.param.sample_sheet != '' && visr.param.gene_annotation != '')")
################################################################

visr.param("output_dir",
           label="Output Directory to save the results",
           info="Output directory where the analysis results will be saved to",
           type="filename", filename.mode = "dir",
           debugvalue= "~/SFU/Datasets/SingleCell/output_monocle/")

visr.param("create_subdir", default = TRUE,
           label = "Create new sub-direcory",
           info = "Create a new sub directory with the name DATE_TIME (YYYYMMDD_hhmmss)")


################################################################
################################################################
################################################################

#' prepares the output directory and output report
prepare_output <- function() {
  # validate or create output subdirectory as needed
  output_dir <<- validateOutputDirectory(visr.param.output_dir, visr.param.create_subdir)

  # output current parameters into a file
  writeParameters(output_dir)

  # prepare the pdf to save the plots to
  output_plot_file <- startReport(output_dir)

  return(output_dir)
}

#' loads an RDS object
load_rds_object <- function(path){
  object <- tryCatch(readRDS(path),
                      error = function(e) {
                        visr.message(sprintf("Cannot read the input object format '%s'.\nMake sure the object is saved using the 'saveRDS()' function.", path))
                      })
  return(object)
}

#'
#' Load input data
#' @return monocle_app_object
#'
load_input <- function() {
  enable_estimations <- TRUE # estimate dispersions
  if (visr.param.input_type == INPUT_TYPE_10X) { # if data is 10x data
    visr.assert_file_exists(visr.param.data_dir_10x, '10X data directory')
    visr.assert_that(file.exists(paste0(visr.param.data_dir_10x, '/outs')), msg = paste("Cannot find 'outs' sub-directory inside the specified 10x directory:\n", visr.param.data_dir_10x))
    visr.librarySource("cellrangerRkit", "http://cf.10xgenomics.com/supp/cell-exp/rkit-install-2.0.0.R")

    visr.logProgress("Loading data from the 10x pipeline")
    gbm <- load_cellranger_matrix(visr.param.data_dir_10x)
    # class(gbm)
    # dim(Biobase::exprs(gbm))
    Biobase::exprs(gbm)[1:10, 1:4] # the first 10 genes in the first 4 cells

    # dim(pData(gbm)) # the phenotypic data
    # head(pData(gbm))

    # dim(fData(gbm)) # the feature information
    # head(fData(gbm)) # this table will be useful for matching gene IDs to symbols

    # create a CellDataSet object
    # note: Monocle expects that the gene symbol column in the feature data is called gene_short_name
    my_feat <- fData(gbm)
    names(my_feat) <- c('id', 'gene_short_name') # rename gene symbol column
    familyFunction <- negbinomial.size() # appropriate expression family for UMI data
    if (visr.param.data_type == DATA_TYPE_FPKM) {
      #TODO: http://cole-trapnell-lab.github.io/monocle-release/docs/#converting-tpm-fpkm-values-into-mrna-counts-alternative
      visr.message("DATA_TYPE_FPKM not implemented")
      familyFunction <- tobit() #TODO: parameter for tobit()
      enable_estimations <- FALSE
    } else if (visr.param.data_type == DATA_TYPE_LOG_FPKM) {
      # TODO
      visr.message("DATA_TYPE_LOG_FPKM not implemented")
      familyFunction <- gaussianff()
      enable_estimations <- FALSE
    }

    #TODO:If your data contains relative counts (e.g. FPKM or TPM values), use relative2abs() to convert these measurements into absolute counts: http://cole-trapnell-lab.github.io/monocle-release/docs/#converting-tpm-fpkm-values-into-mrna-counts-alternative
    visr.logProgress("Creating CellDateSet object")
    my_cds <- newCellDataSet(Biobase::exprs(gbm),
                              phenoData = new("AnnotatedDataFrame", data = pData(gbm)),
                              featureData = new("AnnotatedDataFrame", data = my_feat),
                              lowerDetectionLimit = 0.5,
                              expressionFamily = familyFunction)
    # my_cds
    # slotNames(my_cds)
  }
  else if (visr.param.input_type == INPUT_TYPE_MONOCLE_OBJECT) {
    if (!file.exists(visr.param.path_to_monocle_object))
      visr.message(sprintf("'%s' is not a valid path to a monocle object", visr.param.path_to_monocle_object))
    visr.logProgress("Loading data from the moncle app object ...")
    monocle_app_object <- load_rds_object(visr.param.path_to_monocle_object)
    return(monocle_app_object)
  }
  else if (visr.param.input_type == INPUT_TYPE_COUNT_MATRIX_TXT) {
    visr.assert_file_exists(visr.param.expression_matrix, 'expression matrix')
    visr.assert_file_exists(visr.param.sample_sheet, 'sample sheet')
    visr.assert_file_exists(visr.param.gene_annotation, 'gene annotation')

    visr.logProgress("Loading data from the expression matrix")

    expr_matrix <- read.table(visr.param.expression_matrix) # e.g. "fpkm_matrix.txt"
    sample_sheet <- read.delim(visr.param.sample_sheet) # e.g. "cell_sample_sheet.txt"
    gene_annotation <- read.delim(visr.param.gene_annotation) # e.g. "gene_annotations.txt"

    visr.assert_that(nrow(expr_matrix) == ncol(sample_sheet), msg = "Expression value matrix does not have the same number of columns as the sample sheet (pheno data) has rows.")
    visr.assert_that(nrow(expr_matrix) == nrow(gene_annotation), msg = "Expression value matrix does not have the same number of rows as the gene_annotation (feature data) data frame has rows.")
    visr.assert_that(all(rownames(expr_matrix) == colnames(sample_sheet)), msg = "row names of the sample sheet (pheno data) does not match the column names of the expression matrix.")
    visr.assert_that(all(rownames(expr_matrix) == rownames(gene_annotation)), msg = "row names of the gene annotation (feature data) should match row names of the expression matrix.")
    visr.assert_that("gene_short_name" %in% rownames(gene_annotation), msg = "one of the columns of the gene annotation (feature data) should be named 'gene_short_name'")

    pd <- new("AnnotatedDataFrame", data = sample_sheet)
    fd <- new("AnnotatedDataFrame", data = gene_annotation)
    my_cds <- newCellDataSet(as.matrix(expr_matrix), phenoData = pd, featureData = fd)

  } else if (visr.param.input_type == INPUT_TYPE_SPARSE_MATRIX) {
    visr.message("input_type INPUT_TYPE_SPARSE_MATRIX Not Implemented Yet")
  }

  ############################################################
  #  perform normalization and variance estimation steps
  ############################################################

  if (enable_estimations) {
    # perform normalization and variance estimation steps, which will be used in the differential expression analyses later on.
    # estimateSizeFactors() and estimateDispersions() will only work, and are only needed, if you are working with a CellDataSet
    # with a negbinomial() or negbinomial.size() expression family.
    visr.logProgress("Performing normalization and variance estimation (takes a few minutes)")
    my_cds <- estimateSizeFactors(my_cds)
    my_cds <- estimateDispersions(my_cds)
  }

  monocle_app_object <- list(cds = my_cds)
  return(monocle_app_object)
}


finalize_output <- function(output_dir, monocle_app_object) {
  finishReport()

  browseURL(getOutputPlotFile(output_dir))

  cds_dims <- data.frame(t(monocle::reducedDimA(monocle_app_object$cds)))
  colnames(cds_dims) <- c("Component1", "Component2")
  cell_table <- cbind(pData(monocle_app_object$cds), cds_dims)
  visr.writeDataTable(cell_table, paste0(output_dir, "/cells.txt"))

  visr.writeDataTable(monocle_app_object$disp_table, paste0(output_dir, "/genes_dispersion.txt"))

  monocle_app_object$params <- visr.getParams()
  saveRDS(monocle_app_object, file = paste0(output_dir, "/monocle_app_object.RDS"))
}
