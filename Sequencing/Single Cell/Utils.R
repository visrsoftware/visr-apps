
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


#' return name of the plot file
getOutputPlotFile <- function(output_dir) {
  return(paste0(output_dir, "/plots.pdf"))
}

# prepare the pdf to save the plots to
startReport <- function(output_dir) {
  assign("visr.dev", dev.cur(), envir = .GlobalEnv)

  pdf_dev <- get0("visr_app_pdf_dev", ifnotfound = NULL, envir = .GlobalEnv)
  if (!is.null(pdf_dev)) {
    dev.off(which = pdf_dev) # close the previous one
  }
  output_plot_file <- getOutputPlotFile(output_dir)
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

#' Writes the first two columns of a table as a new plot page
#' @param dataTable input data table
#' @param title     title to be shown at the top of the page
plotTableSummary <- function(dataTable, title) {
  plot.new()
  mtext(text = title, adj=0.5, side=3, line = 1)
  mtext(text=dataTable[,1], adj=0, col=c("gray10","gray60"), side=3, line=c(1:nrow(dataTable))*-1)
  mtext(text=dataTable[,2], adj=1, col=c("gray10","gray60"), side=3, line=c(1:nrow(dataTable))*-1)
}

#' Writes a title into a new plot page
#' @param title   title of the page
#' @param size    font size
plotTitlePage <- function(title, text_size = 3, text_color = "black") {
  # plot.new()
  # mtext(text = title, padj = 0.5, adj=0.5, side=3, line = -10, cex = size)
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, title, cex = text_size, col = text_color)
}

plotGBMSummary <- function(gbm) {
  gbm_summary <- cbind(rownames(t(gbm@summary)), t(gbm@summary))
  colnames(gbm_summary) <- c("name", "value")
  rownames(gbm_summary) <- NULL
  plotTableSummary(gbm_summary, "Input data summary")
  return(gbm_summary)
}

extract_param <- function(path){
  con <- file(path, open = "r")
  lines <- readLines(con = con)
  close(con = con)

  n <- length(lines)
  i <- 1
  op <- "\\("
  cp <- ")"
  params <- list()

  read_complete_line <- function(line,i,start){
    while(!endsWith(line,cp)){
      i <- i+1
      line <- paste(line,lines[i],sep = "")
    }
    res <- list()
    res$i <- i
    start <- substr(start,1,nchar(start)-1)
    line <- gsub(line,pattern = paste(start,op,sep = ""),replacement = "list(")
    args <- eval(parse(text = line))
    res$name <- args[[1]]
    if (is.null(args$active.condition) || eval(parse(text = args$active.condition))){
      res$active <- T
    } else {
      res$active <- F
    }
    if (is.null(args$label)){
      res$label <- gsub(res$name,pattern = "_", replacement = " ")
    } else {
      res$label <- args$label
    }
    return(res)
  }
  while (T){
    line <- lines[i]
    if (startsWith(line,"visr.app.end")){
      break
    }
    for (start in c("visr.app.category(","visr.category(")){
      if (startsWith(line,start)){
        res <- read_complete_line(line,i,start)
        category <- res$name
        i <- res$i
        params[[category]] <- cbind(name = category, label = res$label, value = "",used = res$active)
      }
    }
    if (startsWith(line,"visr.param(")){
      res <- read_complete_line(line,i,"visr.param(")
      param <- res$name
      i <- res$i
      value <- eval(parse(text = paste("visr.param",param,sep = ".")))
      if (is.null(value)){
        value <- ""
      }
      if (params[[category]][1,ncol(params[[category]])] == F){res$active <- F}
      params[[category]] <- rbind(params[[category]],c(param, res$label, value,res$active))
    }
    i <- i+1
  }

  # output <- rbind(c("Name","Label","Value","Used"))
  # for (i in 1:length(params)){
  #   output <- rbind(output, params[[i]], rep("",ncol(output)))
  # }
  return(params)
}

#' get library ids from barcodes
get_lib_id <- function(barcodes){
  get_id <- function(barcode){
    id <- strsplit(x = barcode,split = "-")[[1]][2]
    if (is.na(id)) {id <- "1"}
    return(id)
  }
  lib_ids <- sapply(X = barcodes, FUN = get_id, USE.NAMES = F)
  return(lib_ids)
}

#' Plain format (no scientific) to be used for plot axis scales
#' example ggplot(...) + scale_x_log10(labels = axis_plain_format) + ...
axis_plain_format <- function(x,...)
{
  format(x, ..., scientific = FALSE, trim = TRUE)
}
