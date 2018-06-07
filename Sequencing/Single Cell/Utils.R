
cleanEnvir <- function(pattern) {
  # https://stackoverflow.com/questions/4837477/remove-objects-in-globalenv-from-within-a-function/44280757
  objs <- ls(pos = ".GlobalEnv")
  rm(list = objs[grep(pattern, objs)], pos = ".GlobalEnv")
}

validateOutputDirectory <- function(output_dir, create_subdir) {
  visr.assert_that(!is.null(output_dir) && output_dir != "", msg = "Output directory not specified")
  if (create_subdir)
    output_dir <- paste0(output_dir, format(Sys.time(), "/%Y%m%d_%H%M%S"))
  dir.create(output_dir, recursive = T, showWarnings = F)
  visr.assert_that(dir.exists(output_dir), msg = sprintf("Output directory '%s' is not valid", output_dir))
  return(output_dir)
}

writeParameters <- function(output_dir) {
  write(capture.output(print(visr.getParams())), file = paste0(output_dir, "/visr_params.txt"))
}

# prepare the pdf to save the plots to
startReport <- function(output_dir) {
  assign("visr.dev", dev.cur(), envir = .GlobalEnv)

  pdf_dev <- get0("visr_app_pdf_dev", ifnotfound = NULL, envir = .GlobalEnv)
  if (!is.null(pdf_dev)) {
    dev.off(which = pdf_dev) # close the previous one
  }
  output_plot_file <- paste0(output_dir, "/plots.pdf")
  pdf(file=output_plot_file)
  assign("visr_app_pdf_dev", dev.cur(), envir = .GlobalEnv)
  return(output_plot_file)
}

switchPlotToScreen <- function() {
  visr_dev <- get0("visr.dev", ifnotfound = NULL, envir = .GlobalEnv)
  if (!is.null(visr_dev))
    dev.set(which = visr_dev)
}

switchPlotToReport <- function() {
  pdf_dev <- get0("visr_app_pdf_dev", ifnotfound = NULL, envir = .GlobalEnv)
  if (!is.null(pdf_dev))
    dev.set(which = pdf_dev)
}

finishReport <- function() {
  pdf_dev <- get0("visr_app_pdf_dev", ifnotfound = NULL, envir = .GlobalEnv)
  if (!is.null(pdf_dev)) {
    dev.off(which=pdf_dev)
    cleanEnvir("visr_app_pdf_dev")
  }
}

plotTableSummary <- function(dataTable, title) {
  plot.new()
  mtext(text = title, adj=0.5, side=3, line = 1)
  mtext(text=dataTable[,1], adj=0, col=c("gray10","gray60"), side=3, line=c(1:nrow(dataTable))*-1)
  mtext(text=dataTable[,2], adj=1, col=c("gray10","gray60"), side=3, line=c(1:nrow(dataTable))*-1)
}

plotGBMSummary <- function(gbm) {
  gbm_summary <- cbind(rownames(t(gbm@summary)), t(gbm@summary))
  colnames(gbm_summary) <- c("name", "value")
  rownames(gbm_summary) <- NULL
  plotTableSummary(gbm_summary, "Input data summary")
  return(gbm_summary)
}

