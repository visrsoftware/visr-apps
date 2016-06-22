error_message <- ""

visr.isGUI<-function(){
  if (exists("visr.var.isGUI")) {
    if (visr.var.isGUI)
      return(TRUE)
  }
  return(FALSE)
}

visr.applyParameters <- function() {
  #if (exists("visr.var.message.ignore")) rm(visr.var.message.ignore)
  dummylocalvar<-"dummyvalue"
}

visr.message<-function(text, type=c("error","warning"))
{
  #TODO: replace error_message with visr.var.message
  if (exists("error_message") && is.character(error_message) && nzchar(error_message)) {
    # There is an unhandled error message already. Concatenate this to it
    error_message <<- paste(error_message,"\n", match.arg(type), ": ", text, sep = "")
  } else {
    error_message <<- paste(match.arg(type),": ", text, sep = "")
  }

  if (!visr.isGUI() && !exists("visr.var.message.ignore")) {
    print(error_message)
    error_message <<- ""
    invisible(user_choice<-readline(prompt="(s)top / (i)gnore / ignore (a)ll ? (s/i/a)"))
    if (user_choice == "s")
      stop("Terminated", call. = FALSE, domain = NA)
    else if (user_choice == "a")
      visr.var.message.ignore <<- TRUE
  }
}

# Loads a CRAN package. If not already installed, tries to install the package from CRAN.
visr.library<-function (pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org", dependencies=TRUE)
    #update.packages(ask = TRUE)
    # adding some delay, since loading the package right after installation may not work.
    numtries=10
    while (numtries > 0 && !require(pkg, character.only = TRUE)) {
      Sys.sleep(0.1)
      numtries=numtries-1
    }

    if (!require(pkg, character.only = TRUE)) {
      visr.message(paste("Unable to load package", pkg))
    }
  }
}

visr.libraryURL<-function (pkg,url) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(url, repos = NULL, type="source", dependencies=TRUE)
    # adding some delay, since loading the package right after installation may not work.
    numtries=10
    while (numtries > 0 && !require(pkg, character.only = TRUE)) {
      Sys.sleep(0.1)
      numtries=numtries-1
    }

    if (!require(pkg, character.only = TRUE)) {
      visr.message(paste("Unable to load package", pkg))
    }
  }
}

# Loads a bioconductor package. If not already installed, tries to install the package from bioconductor.
visr.biocLite<-function (pkg) {
  if (!require(pkg, character.only = TRUE)) {
    source("http://bioconductor.org/biocLite.R")
    biocLite(pkg,
             suppressUpdates=FALSE,     # suppress automatic updating of all installed packages.
             suppressAutoUpdate=FALSE,  # whether the BiocInstaller package updates itself.
             ask=FALSE)                 # whether to prompt user before installed packages are updated
    numtries=10
    while (numtries > 0 && !require(pkg, character.only = TRUE)) {
      Sys.sleep(0.1)
      numtries=numtries-1
    }

    if (!require(pkg, character.only = TRUE)) {
      visr.message(paste("Unable to load package", pkg))
    }
  }
}

# used in tryCatch
visr.internal.handleError <- function(e)
{
  #todo, use a different variable
  error_message <<- e$message
}

# used in tryCatch
visr.internal.handleWarning <- function(w)
{
  #todo, use a different variable
  error_message <<- w$message
}

visr.rebuildPackages <- function()
{
  #biocLite()
  #Update all/some/none? [a/s/n]:
  #  a
  #Do you want to install from sources the packages which need compilation?
  #y/n: n

  pkgs = installed.packages()
  idx = pkgs[,"Built"] != "3.2.2"
  for (pkg in rownames(pkgs[idx,]))
  {
    visr.biocLite(pkg)
  }
}

# returns the user home directory
visr.getHomeDir <- function()
{
  if (.Platform$OS.type == "windows") {
    return(Sys.getenv("UserProfile"))
  }
  return(Sys.getenv("HOME"))
}

# returns the user library directory
visr.getHomeLibPath <- function()
{
  return(paste(visr.getHomeDir(),"/VisRseq/RLibs",sep=""))
}

#the function currently wraps print so that it doesn't print in VisRseq causing the SIGPIPE error.
visr.print<-function(msg) {
  if (!visr.isGUI())
    print(msg)
}


# utility function to convert column names to correct format by replacing invalid characters with _
visr.internal.validname <- function(columnnames) {
  return (make.names(gsub("[^a-zA-Z0-9_]", "_", columnnames)))
}

# utility function to open data tables with corrected column names (useful for debugging within R studio)
visr.readDataTable <-function(file) {
  t <- read.csv(file, sep = "\t", check.names = F)
  colnames(t) <- visr.internal.validname(colnames(t))
  return (t)
}

visr.setLogDir <- function(logDir) {
  if (FALSE) {
    visr.var.logDir <<- logDir;
    if (nchar(logDir) > 0) {
      sinkFile <- file(paste(logDir,"/all.txt",sep=""), open = "wt")
      sink(sinkFile)
      sink(sinkFile, type = "message")
      #print(date())
    } else {
      ## back to the console
      if (sink.number() > 0) {
        #print(date())
        sink(type = "message")
        sink()
      }
    }
  }
}


###############################################
#                  visr.app
#
# functions to create an app json file in R
###############################################

visr.var.appJSON <- ""
visr.var.definedCategory <- FALSE
visr.var.definedParam <- FALSE

# indents all lines of a string where lines are separated by \n
visr.internal.indent <- function(txt, indents=2) {
  spaces <- paste(rep(" ", indents), collapse="")
  paste(spaces, gsub("\n", paste("\n", spaces, sep=""), txt), sep="")
}

# appends @param(txt) to the current json output
visr.internal.appendJSON <- function(txt) {
  visr.var.appJSON <<- paste(visr.var.appJSON, txt, sep="")
}

#' Starts definition of parameters for an R app
#' @param  name app name
#' @param  info  app info shown as tooltip
#' @param  debugdata  debug dataframe to be used when debuggin the app in R / RStudio
visr.app.start <- function(name, info = "", debugdata = NULL) {
  if (visr.isGUI())
    return()

  visr.var.appJSON <<- paste('{\n  "label": "', name, '",\n  "info": "', info, '",\n  "categories":[', sep='')
  visr.var.definedCategory <<- FALSE
  visr.var.definedParam <<- FALSE

  visr.input  <<- debugdata
  input_table <<- debugdata
}

#' finishes the current apps parameter definition
#' @param printjson     whether to print the generated json file to console
#' @param writefile     whether to write the generated json to a file
#' @param filename      path to the filename to write the json to. If not specified and
#'                      writefile is TRUE, a json file is generated from the caller
#'                      source file path by replacing .R with .json
visr.app.end <- function(printjson = FALSE, writefile = FALSE, filename = NULL) { # preview=FALSE
  if (visr.isGUI())
    return()

  if (visr.var.definedCategory)
    visr.internal.appendJSON('\n    }\n  }')

  visr.internal.appendJSON(']\n}')

  if (printjson) {
    cat(visr.var.appJSON)
  }

  if (writefile) {
    if (is.null(filename)) {
      # auto generate the json file name from the R source filename
      srcfilename <- parent.frame(3)$ofile
      if (!is.null(srcfilename))
        filename <- paste(dirname(srcfilename), "/", gsub("\\.R", ".json", basename(srcfilename)), sep="")
    }

    if (!is.null(filename)) {
      print(paste("Writing app parameter description to", filename))
      write(visr.var.appJSON, file=filename)
    }
  }
}

#' Starts a new category of parameters for the app
#' @param label   category label to be shown in VisRseq
#' @param info    additional information about category shown as tooltip
visr.category <- function(label, info = "") {
  if (visr.isGUI())
    return()

  if (visr.var.definedCategory)
    visr.internal.appendJSON('\n    }\n  },\n')

  visr.internal.appendJSON(paste('  {\n    "label": "', label, '",\n    "info": "', info, '",\n    "variables": {\n', sep=""))
  visr.var.definedCategory <<- TRUE
  visr.var.definedParam <<- FALSE
}

#' Adds a new parameter to the app
#'
#' @param name    parameter name. will be appended to "visr.param." to create the full variable name in R
#' @param label   the label for the parameter's GUI control. will use name if NULL
#' @param info    additional information about parameter shown as tooltip in VisRseq
#' @param default default value for the parameter
visr.param <- function(name, label = NULL, info = NULL,
                       default = NULL, min = NULL, max = NULL,
                       items = NULL, item.labels = NULL,
                       filename.mode = c("load", "save", "dir"),
                       type = c("string", "character", "int", "integer", "double", "boolean", "logical", "multi-string",
                               "column", "multi-column", "column-numerical", "multi-column-numerical",
                               "color", "multi-color", "filename",
                               "output-column", "output-multi-column", "output-table"),
                       debugvalue = NULL
                      ) {
  if (visr.isGUI())
    return()

  paramname = paste("visr.param", name, sep=".") #full parameter name

  if (missing(type) && !is.null(default)) {
    #guess type from the default value
    if (is.numeric(default) && (default %% 1 == 0))
      type <- "int"
    else if (is.numeric(default))
      type <- "double"
    else if (is.logical(default))
      type <- "boolean"
  }
  type <- match.arg(type)

  type <- if (type == "character") {"string"} else {type}
  type <- if (type == "logical") {"boolean"} else {type}
  type <- if (type == "integer") {"int"} else {type}

  if (!is.null(min) && !is.numeric(min))
    stop("argument min should be numeric")

  if (!is.null(max) && !is.numeric(max))
    stop("argument max should be numeric")

  if (type == "filename") {
    filename.mode <- match.arg(filename.mode)
  } else {
    if (!missing(filename.mode)) {
      warning("filename.mode is ignored when type != 'filename'")
    } else {
      filename.mode <- NULL
    }
  }

  # check that type matches default
  if (!is.null(default)) {
    if (((type=="int" || type=="double") && !is.numeric(default)) ||
        ( type=="boolean" && !is.logical(default)) ||
        ( type=="string" && !is.character(default)))
      stop ("default value does not match the type")

    if (type == "color") {
      default <- paste('#', paste(as.hexmode(col2rgb(default)), collapse=""), '', sep='')
      #if (length(default) > 1) default <- apply(as.character(as.hexmode(col2rgb(default))), 2, paste, collapse="")
    }
  }

  # try to guess a debug value from the other properties
  if (is.null(debugvalue)) {
    if (!is.null(default)) {
      debugvalue <- default
    } else if (!is.null(items)) {
      debugvalue <- items[1]
    } else if (type == "multi-column" && !is.null(visr.input)) {
      debugvalue <- colnames(visr.input)
    } else if (type == "multi-column-numerical" && !is.null(visr.input)) {
      debugvalue <- colnames(visr.input)[which(sapply(visr.input, is.numeric))]
    }
  }
  assign(paramname, debugvalue, envir = .GlobalEnv)

  default.ischar <- is.character(default)
  if (type == "boolean" && !is.null(default))
    default <- tolower(default)

  properties <- c(
    if (!is.null(label))    {paste('"label": ',   label,    '', sep='"')} else {NULL},
    if (!is.null(info))     {paste('"info": ',    info,     '', sep='"')} else {NULL},
    if (!is.null(type))     {paste('"type": ',    type,     '', sep='"')} else {NULL},
    if (!is.null(default))  {paste('"default": ', default , '', sep = ifelse(default.ischar, '"', ''))} else {NULL},
    if (!is.null(min))      {paste('"min": ',     min,          sep='')} else {NULL},
    if (!is.null(max))      {paste('"max": ',     max,          sep='')} else {NULL},
    if (!is.null(items))       {paste('"items": ',       paste('[', paste('"', items,       '"', sep="", collapse=",") ,"]"), sep='')} else {NULL},
    if (!is.null(item.labels)) {paste('"item-labels": ', paste('[', paste('"', item.labels, '"', sep="", collapse=",") ,"]"), sep='')} else {NULL},
    if (!is.null(filename.mode)) {paste('"filename.mode": ', filename.mode, '', sep='"')} else {NULL}
  )

  properties <- properties[which(!is.na(properties))]

  jsonstr <- paste(
    paste('"', paramname,'": {\n', sep=''),
    visr.internal.indent(paste(properties, collapse = ",\n")),
    '\n}', sep=''
  )

  if (!visr.var.definedCategory)
    visr.category(label="")

  if (visr.var.definedParam)
    visr.internal.appendJSON(",\n")
  visr.internal.appendJSON(visr.internal.indent(jsonstr, 6))
  visr.var.definedParam <<- TRUE
}

#unit test
visr.internal.test.param <- function() {
  visr.app.start("test-app", info="A test app", debugdata = iris)
  visr.param("test-minimal")
  visr.param("test-auto-int", default = 2)
  visr.param("test-auto-double", default = 0.5)
  visr.param("test-auto-bool", default = FALSE)
  visr.category("group2", "info for group2")
  #visr.param("test-mismatch", type="char", default=2)
  visr.param("test-color", label="foreground", info="foreground color", type="color", default = "yellow")
  visr.param("test-min-max", default=3, min=1, max=10)
  visr.param("test-filename", type="filename", filename.mode = "load")
  visr.param("test-items", items = c("i1","i2","i3"))
  visr.app.end(printjson=TRUE, filename="~/testapp.json")
  #cat(visr.var.appJSON)
}

#.libPaths(c(visr.getHomeLibPath(), .libPaths()))