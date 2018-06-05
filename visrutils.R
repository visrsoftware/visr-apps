assign("error_message", "", .GlobalEnv)

#' Check if running in VisR GUI
#'
#' Checks if running from within the VisR
#' @return  Will return \code{TRUE} if the code is running from within VisR.
#' @export
visr.isGUI<-function(){
  return (get0("visr.var.isGUI", ifnotfound = FALSE, envir = .GlobalEnv))
}

if (visr.isGUI()) {
  # https://stat.ethz.ch/R-manual/R-devel/library/base/html/options.html
  options(install.packages.check.source = "no") # "yes" "no"
  options(install.packages.compile.from.source = "never") # "interactive", "never", "always"
  options(pkgType = "both") # "both", "source", "binary", "mac.binary.el-capitan",  "win.binary"
  # interactive() ?
}

visr.getSelectedRows<-function() {
  selected = get0("visr.input.selectedRows", ifnotfound = NULL, envir = .GlobalEnv)
  if (!is.null(selected)) {
    selected <- selected + 1 # VisR is 0-based
  }
  return(selected)
}

visr.setSelectedRows<-function(x) {
  assign("visr.input.selectedRows", x, .GlobalEnv)
}

visr.get <- function(paramName) {
  return(get0(paste0("visr.param.",paramName), ifnotfound = NULL, envir = .GlobalEnv))
}

#' Apply GUI parameters
#'
#' Applies (imports) parameters spcified in the VisR GUI to the R environment.
#' @export
visr.applyParameters <- function() {
  #if (exists("visr.var.message.ignore")) rm(visr.var.message.ignore)
  dummylocalvar<-"dummyvalue"
  if (!visr.isGUI()) {
    print(visr.params)
  }
}

visr.getParams <- function() {
  return (get0("visr.params", ifnotfound = NULL, envir = .GlobalEnv))
}

readline.orig <- base::readline
visr.readline <- function(prompt = "", ...) {
  cat("visr.readline")
  if (!visr.isGUI()) {
    return (readline.orig(prompt, ...))
  } else {
    visr.message(prompt, "prompt")
    return("\n")
  }
}

if (visr.isGUI()) {
    # disable readline when running inside visr

    #unlockBinding("readline", .GlobalEnv)
    #readline <- visr.readline
    assign("readline", visr.readline, .GlobalEnv)
    #assignInNamespace("readline", visr.readline, ns="base")
}


#' Show message dialog
#'
#' Shows a message dialog to the user in VisR.
#' @param msg   Message text to be shown.
#' @param type   Message type. (\code{"error"} or \code{"warning"} or \code{"prompt"})
#' @examples
#' if (any(is.na(visr.input)))
#'     visr.message("There are NA values in the input", type="error")
#' @export
visr.message<-function(msg, type=c("error","warning", "prompt"))
{
  #TODO: replace error_message with visr.var.message
  if (exists("error_message") && is.character(error_message) && nzchar(error_message)) {
    # There is an unhandled error message already. Concatenate this to it
    assign("error_message",
           paste(error_message,"\n", match.arg(type), ": ", msg, sep = ""),
           .GlobalEnv)
  } else {
    assign("error_message",
           paste(match.arg(type),": ", msg, sep = ""),
           .GlobalEnv)
  }

  if (visr.isGUI() && type == "error") {
      stop(msg)
  }

  if (!visr.isGUI() && !exists("visr.var.message.ignore")) {
    print(error_message)
    assign("error_message", "", .GlobalEnv)
    invisible(user_choice<-readline.orig(prompt="(s)top / (i)gnore / ignore (a)ll ? (s/i/a)"))
    if (user_choice == "s") {
      stop("Terminated", call. = FALSE, domain = NA)
    } else if (user_choice == "a") {
      assign("visr.var.message.ignore", TRUE, .GlobalEnv)
    }
  }
}

visr.require <- function(pkg) {
  return (suppressWarnings(require(pkg, character.only = TRUE, quietly = TRUE)))
}

# Loads a CRAN package. If not already installed, tries to install the package from CRAN.
visr.library<-function (pkg, quiet = F, repos = NULL) {
  if (!quiet)
    visr.logProgress(paste("Loading package:", pkg))
  if (!visr.require(pkg)) {
    if (!quiet)
      visr.logProgress(paste("Installing cran package:", pkg))
    if (is.null(repos))
      repos = "http://cran.us.r-project.org"
    install.packages(pkg, repos = repos, dependencies=TRUE)
    #update.packages(ask = TRUE)
    # adding some delay, since loading the package right after installation may not work.
    numtries=10
    while (numtries > 0 && !visr.require(pkg)) {
      Sys.sleep(0.1)
      numtries=numtries-1
    }

    if (!visr.require(pkg)) {
      visr.message(paste("Unable to load package", pkg))
      if (!quiet)
        visr.logProgress(paste("Failed loading package:", pkg))
    } else {
      if (!quiet) {
        visr.logProgress(paste(pkg, 'version', packageVersion(pkg)))
        visr.logProgress("");
      }
    }
  } else {
    if (!quiet) {
      visr.logProgress(paste(pkg, 'version', packageVersion(pkg)))
      visr.logProgress("");
    }
  }
}

visr.libraryURL<-function (pkg,url) {
  visr.logProgress(paste("Loading package:", pkg))
  if (!visr.require(pkg)) {
    visr.logProgress(paste("Installing package", pkg,"from", url))
    install.packages(url, repos = NULL, type="source", dependencies=TRUE)
    # adding some delay, since loading the package right after installation may not work.
    numtries=10
    while (numtries > 0 && !visr.require(pkg)) {
      Sys.sleep(0.1)
      numtries=numtries-1
    }

    if (!visr.require(pkg)) {
      visr.message(paste("Unable to load package", pkg))
      visr.logProgress(paste("Failed loading package:", pkg))
    } else {
      visr.logProgress(paste(pkg, 'version', packageVersion(pkg)))
      visr.logProgress("");
    }
  } else {
    visr.logProgress(paste(pkg, 'version', packageVersion(pkg)))
    visr.logProgress("");
  }

}

visr.librarySource<-function (pkg, url) {
  visr.logProgress(paste("Loading package:", pkg))
  if (!visr.require(pkg)) {
    visr.logProgress(paste("Installing package", pkg,"from", url))
    source(url)
    # adding some delay, since loading the package right after installation may not work.
    numtries=10
    while (numtries > 0 && !visr.require(pkg)) {
      Sys.sleep(0.1)
      numtries=numtries-1
    }

    if (!visr.require(pkg)) {
      visr.message(paste("Unable to load package", pkg))
      visr.logProgress(paste("Failed loading package:", pkg))
    } else {
      visr.logProgress(paste(pkg, 'version', packageVersion(pkg)))
      visr.logProgress("");
    }
  } else {
    visr.logProgress(paste(pkg, 'version', packageVersion(pkg)))
    visr.logProgress("");
  }
}

visr.libraryGithub<-function (pkg, repo) {
  visr.logProgress(paste("Loading package:", pkg))
  if (!visr.require(pkg)) {
    visr.logProgress(paste("Installing package", pkg,"from github repo: ", repo))
    visr.library("devtools")
    install_github(repo) # https://github.com/cole-trapnell-lab/monocle-release/issues/118

    # adding some delay, since loading the package right after installation may not work.
    numtries=10
    while (numtries > 0 && !visr.require(pkg)) {
      Sys.sleep(0.1)
      numtries=numtries-1
    }

    if (!visr.require(pkg)) {
      visr.message(paste("Unable to load package", pkg))
      visr.logProgress(paste("Failed loading package:", pkg))
    } else {
      visr.logProgress(paste(pkg, 'version', packageVersion(pkg)))
      visr.logProgress("");
    }
  } else {
    visr.logProgress(paste(pkg, 'version', packageVersion(pkg)))
    visr.logProgress("");
  }
}

#' Load/install app dependent package
#'
#' Loads a bioconductor package. If not already installed, tries to install the package from bioconductor.
#' @param pkg package name
visr.biocLite<-function (pkg) {
  visr.logProgress(paste("Loading bioconductor package:", pkg))
  if (!visr.require(pkg)) {
    visr.logProgress(paste("Installing bioconductor package:", pkg))
    source("http://bioconductor.org/biocLite.R")
    biocLite(pkg,
             suppressUpdates=FALSE,     # suppress automatic updating of all installed packages.
             suppressAutoUpdate=FALSE,  # whether the BiocInstaller package updates itself.
             ask=FALSE)                 # whether to prompt user before installed packages are updated
    numtries=10
    while (numtries > 0 && !visr.require(pkg)) {
      Sys.sleep(0.1)
      numtries=numtries-1
    }

    if (!visr.require(pkg)) {
      visr.message(paste("Unable to load package", pkg))
      visr.logProgress(paste("Failed loading package:", pkg))
    } else {
      visr.logProgress(paste(pkg, 'version', packageVersion(pkg)))
      visr.logProgress("");
    }
  } else {
    visr.logProgress(paste(pkg, 'version', packageVersion(pkg)))
    visr.logProgress("");
  }
}

# used in tryCatch
visr.internal.handleError <- function(e)
{
  #todo, use a different variable
  assign("error_message", e$message, .GlobalEnv)
}

# used in tryCatch
visr.internal.handleWarning <- function(w)
{
  #todo, use a different variable
  assign("error_message", w$message, .GlobalEnv)
}

visr.rebuildPackages <- function()
{
  #biocLite()
  #Update all/some/none? [a/s/n]:
  #  a
  #Do you want to install from sources the packages which need compilation?
  #y/n: n

  pkgs = installed.packages()
  idx = which(pkgs[,"Built"] != getRversion())
  if (length(idx)) {
    for (pkg in rownames(pkgs[idx,]))
    {
      visr.biocLite(pkg)
    }
  }
}

#' Returns the user home directory.
#' @export
visr.getHomeDir <- function()
{
  if (.Platform$OS.type == "windows") {
    return(Sys.getenv("UserProfile"))
  }
  return(Sys.getenv("HOME"))
}

#' Returns the user library directory
#'
#' This is typically the ~/VisR/RLibs (*nix) or /Users/[username]/VisR/RLibs (windows).
#' Packages are default installed in this directory to ensure write access permission.
#' @export
visr.getHomeLibPath <- function()
{
  return(paste(visr.getHomeDir(),"/VisR/RLibs",sep=""))
}

# Wraps R's print() function (legacy, as previously print() in VisR was causing a SIGPIPE error)
visr.print<-function(msg) {
  #if (!visr.isGUI())
  print(msg)
}


#' Convert names to syntactically valid names.
#'
#' Utility function to convert column names to syntactically valid names by replacing invalid characters with _
#' @param columnnames vector of character column names
#' @export
visr.makeValidNames <- function(columnnames) {
  return (make.names(gsub("[^a-zA-Z0-9_]", "_", columnnames)))
}

#' Read tables with corrected column names
#'
#' Utility function to read tables with corrected column names, as it is sent from VisR to the apps
#' Mainly useful for debugging apps within R studio
#' @param file  path to the tab delimited file
#' @export
visr.readDataTable <-function(file, delim = "\t") {
  t <- read.csv(file, sep = delim, check.names = F)
  colnames(t) <- visr.makeValidNames(colnames(t))
  return (t)
}

#' Write tables to tab-delimited output files
#'
#' Utility function to write tables to tab-delimited output files that correctly open in VisR
#' @param x
#' data table to write
#' @param filename
#' path output filename
#' @param row.names
#' either logical to indicate outputting row names, or a character vector of row names to be written.
visr.writeDataTable <- function(x, filename, row.names = FALSE) {
  if (isTRUE(row.names) || is.character(row.names)) {
    cn <- colnames(x)
    if (is.character(row.names))
      x <- cbind(row.names, x)
    else
      x <- cbind(rownames(x), x)
    colnames(x) <- c("name", cn)
  }
  write.table(x, file = filename, quote = FALSE, sep="\t", row.names = FALSE)
}


#' Sets the log directory.
#'
#' Sets the current log directory where the R output should be diverted to.
#' @param current log directory. "" to go back to logging to console
#' Warning: Don't call this within an app as it is called internatlly by the VisR
visr.internal.setLogDir <- function(logDir) {
  if (TRUE) {
    assign("visr.var.logDir", logDir, .GlobalEnv)
    assign("visr.var.progressFile", paste(logDir, "/progress.txt", sep=""), .GlobalEnv)
    if (!is.null(logDir) && nchar(logDir) > 0) {
      sinkFile <- file(paste(logDir,"/all.txt",sep=""), open = "wt")
      sink(sinkFile)
      sink(sinkFile, type = "message")
      assign("visr.var.sinkFile", sinkFile, .GlobalEnv)
      #print(date())
    } else {
      ## back to the console
      if (sink.number() > 0) {
        #print(date())
        sink(type = "message")
        sink()
        if (exists("visr.var.sinkFile")) {
          close(visr.var.sinkFile)
          rm(visr.var.sinkFile)
        }
      }
    }
  }
}

#' Logs the current progress of the app execution.
#' @param  message      message to be shown (e.g. the current action being performed)
#' @param  percentage   the current progress percentage of the app (if known)
#' @examples
#' visr.logProgress("Normalizing the input", 0.8)
#' @export
visr.logProgress <- function(message, percentage = -1)
{
  if (exists("visr.var.progressFile")) {
    cat(paste(percentage, message, sep="\t"), file=visr.var.progressFile, append=TRUE, sep = "\n")
  } else {
    if (message != "")
      print(message)
  }
}

#' safely shuts down a graphcis device specified by a character name
#' @param devName device name (character)
#' @return Logical, true if and only the device was found and shut down
visr.internal.dev_off_safe <- function(devName) {
  if (!is.character(devName))
    return(FALSE)
  devValue <- get0(devName, envir = .GlobalEnv)
  if (!is.null(devValue) && !identical(devValue, 1)) {
    dev.off(devValue)
    rm(list=devName, envir = .GlobalEnv)
    return(TRUE)
  }
  return(FALSE)
}

#' safely sets the current graphcis device to the device specified by a character name
#' @param devName device name (character)
#' @return Logical, true if and only the device was found
visr.internal.dev_set_safe <- function(devName) {
  if (!is.character(devName))
    return(FALSE)

  devValue <- get0(devName, envir = .GlobalEnv)
  if (!is.null(devValue) && !identical(devValue, 1) && devValue %in% dev.list()) {
    dev.set(devValue)
    return(TRUE)
  }
  return(FALSE)
}


#' Prints the session info and installed packages (useful for diagnosis)
#' @export
visr.printSessionInfo <- function() {
  cat("\nsessionInfo()\n")
  print(sessionInfo())

  cat("\n.libPaths()\n")
  print(.libPaths())

  cat("\ninstalled.packages()\n")
  ip <- as.data.frame(installed.packages()[,c(1,3:4)])
  rownames(ip) <- NULL
  ip <- ip[is.na(ip$Priority),1:2,drop=FALSE]
  print(ip, row.names=FALSE)

  cat("\ncapabilities()\n")
  print(capabilities())
}

###############################################
#                  visr.app
#
# functions to create an app json file in R
###############################################

assign("visr.var.appJSON", "", .GlobalEnv)
assign("visr.var.definedCategory", FALSE, .GlobalEnv)
assign("visr.var.definedParam", FALSE, .GlobalEnv)

# indents all lines of a string where lines are separated by \n
visr.internal.indent <- function(txt, indents=2) {
  spaces <- paste(rep(" ", indents), collapse="")
  paste(spaces, gsub("\n", paste("\n", spaces, sep=""), txt), sep="")
}

# appends @param(txt) to the current json output
visr.internal.appendJSON <- function(txt) {
  assign("visr.var.appJSON", paste(visr.var.appJSON, txt, sep=""), .GlobalEnv)
}

#' Start app definition
#'
#' Starts definition of parameters for an R app.
#' @param  name app name
#' @param  info  app info shown as tooltip
#' @param  input.type input type of the app. "dragAndDrop" if requires a data node to be drag and dropped. "none" otherwise.
#' @param  debugdata  debug dataframe to be used when debuggin the app in R / RStudio
#' @examples
#' visr.app.start("kmeans", info="kmeans clustering")
#' visr.app.start("kmeans", debugdata = iris) # will assign visr.input <- iris in debug mode
#' @export
visr.app.start <- function(name, info = "", input.type=c("none", "dragAndDrop"), instructions = "", debugdata = NULL) {
  if (visr.isGUI())
    return()

  input.type <- match.arg(input.type)

  assign("visr.var.appJSON",
         paste('{\n  "label": "', name,
               '",\n  "info": "', gsub('"', '\\\\"', gsub('\n', '\\\\n', info)),
               '",\n  "instructions": "', gsub('"', '\\\\"', gsub('\n', '\\\\n', instructions)),
               '",\n  "input": "', input.type,
               '",\n  "categories":[', sep=''), .GlobalEnv)
  assign("visr.var.definedCategory", FALSE, .GlobalEnv)
  assign("visr.var.definedParam", FALSE, .GlobalEnv)

  assign("visr.input", debugdata, .GlobalEnv)
  assign("input_table", debugdata, .GlobalEnv)
  assign("visr.params", list(), .GlobalEnv)
}

#' End app definition
#'
#' Finishes the current apps parameter definition.
#' @param printjson     whether to print the generated json file to console
#' @param writefile     whether to write the generated json to a file
#' @param filename      path to the filename to write the json to. If not specified and
#'                      writefile is TRUE, a json file is generated from the caller
#'                      source file path by replacing .R with .json
#' @export
visr.app.end <- function(printjson = TRUE, writefile = TRUE, filename = NULL) { # preview=FALSE
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
    } else {
      warning("Could not write the json file as filename parameter is not specified and could not be inferred from source R file.")
    }
  }
}

#' Start category
#'
#' Starts a new category of parameters for the app
#' @param label   category label to be shown in VisR
#' @param info    additional information about category shown as tooltip
#' @param collapsed whether the category should be initially shown in collapsed (minimized) mode
#' @examples
#' visr.category("clustering parameters")
#' @export
visr.category <- function(label, info = "", collapsed = FALSE, active.condition = NULL) {
  if (visr.isGUI())
    return()

  if (visr.var.definedCategory)
    visr.internal.appendJSON('\n    }\n  },\n')

  if (!is.null(active.condition))
    eval(parse(text=active.condition)) # to make sure it is valid

  visr.internal.appendJSON(paste('  {\n',
                                 '    "label": "', gsub('"', '\\\\"', label), '",\n',
                                 '    "info": "', gsub('"', '\\\\"', gsub('\n', '\\\\n', info)), '",\n',
                                 '    "collapsed": ', ifelse(collapsed, "true", "false"), ',\n',
                                 if (is.null(active.condition)) {''} else {paste(
                                 '    "active-condition": "', active.condition, '",\n', sep='')},
                                 '    "variables": {\n', sep=""))
  assign("visr.var.definedCategory", TRUE, .GlobalEnv)
  assign("visr.var.definedParam", FALSE, .GlobalEnv)
}
visr.app.category <- visr.category

visr.app.paramNamePrefix = "visr.param."

#' Set parameter name prefix
#'
#' Sets the prefix that is added to parameter name when using visr.param()
#'
#' @param prefix    the name prefix
visr.app.setParamNamePrefix <- function(prefix = "visr.param.") {
  assign("visr.app.paramNamePrefix", prefix, .GlobalEnv)
}


#' Add app parameter
#'
#' Adds a new parameter to the curent app.
#'
#' @param name    Parameter name. will be appended to "visr.param." to create
#'                the full variable name in R.
#' @param label   The label for the parameter's GUI control.
#'                Will use variable name if not specified (NULL).
#' @param info    Additional information about parameter.
#'                Will be shown as tooltip in VisR.
#' @param type    Parameter type
#' @param default Initial default value for the parameter.
#' @param min     Minimum value for numerical type parameters
#'                (\code{"int"}, \code{"range-int"}, \code{"double"}, \code{"range-double"})
#' @param max     Maximum value for numerical type parameters
#'                (\code{"int"}, \code{"range-int"}, \code{"double"}, \code{"range-double"})
#' @param items   Specify a vector of items for a \code{"string"} or \code{"multi-string"} variable to select from.
#' @param item.labels Specify a vector of items to be used as the labels for \code{items} argument.
#' @param filename.mode  Specify the file dialog mode
#'                (file load, file save or directory) for a
#'                \code{"filename"} type parameter
#' @param active.condition  A logical expression to control when this parameter will be visible.
#'                    e.g. visr.param("title_font", active.condition = "visr.param.title != ''")
#' @param options     various options to control the specific behaviors of the parameters.
#'                    Should be specified in the form of "option1=value, option2=value, ..."
#'                    Current defined options:
#'                    "importRowNames" for type="output-table": if importRowNames=false, row names won't be imported. default is true.
#' @param debugvalue  The value to be assigned to the R variable in debug mode.
#'                    Useful for unit testing.
#'
#'
#' @details
#' For parameter type of "multi-color" the default should be in the form of "PaletteName" or "PaletteName  ColorCount".
#' PaletteName is a color brewer palette name which can be found by display.brewer.all() of "RColorBrewer" package.
#' (http://www.sthda.com/english/wiki/colors-in-r#using-rcolorbrewer-palettes)
#'
#'
#' @examples
#' visr.param("k", type = "integer") # specify type
#' visr.param("k", default = 3L) # will infer type from default value
#' visr.param("k", label = "number of clusters")  # explicitly specify label
#' visr.param("title") # no type or default value: treated as a "string" type
#' visr.param("algorithm", items = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"), active.condition="visr.param.k > 1")
#' visr.param("columns", type = "multi-column-numerical")
#' visr.param("output.clusterid", type = "output-column") # column appended to input table
#' @export
visr.param <- function(name, label, info,
                       type = c("string", "character",
                                "int", "integer", "double",
                                "range-int", "range-double",
                                "boolean", "logical", "multi-string",
                                "column", "multi-column",
                                "column-numerical", "multi-column-numerical",
                                "color", "multi-color", "filename",
                                "output-column", "output-multi-column", "output-table"),
                       default, min, max,
                       items, item.labels,
                       filename.mode = c("load", "save", "dir"),
                       active.condition,
                       options,
                       debugvalue) {
  if (visr.isGUI()) # don't generate parameters when running within VisR
    return()

  paramname = paste0(visr.app.paramNamePrefix, name) #full parameter name

  if (missing(default) && !missing(min) && !missing(max) && !missing(type) &&
      (type == 'range-int' || type == 'range-double'))
    default <- c(min, max)

  if (missing(default) && !missing(min)) # guess default from min
    default <- min

  if (missing(default) && !missing(max)) # guess default from max
    default <- max

  if (missing(type) && !missing(default)) { # try to guess type from default
    #guess type from the default value
    if (is.numeric(default) && is.integer(default)) {
      if (length(default) == 2)
        type <- "range-int"
      else
        type <- "int"
    } else if (is.numeric(default)) {
      if (length(default) == 2)
        type <- "range-double"
      else
        type <- "double"
    } else if (is.logical(default)) {
      type <- "boolean"
    }
  }
  type <- match.arg(type)

  type <- if (type == "character") {"string"} else {type}
  type <- if (type == "logical") {"boolean"} else {type}
  type <- if (type == "integer") {"int"} else {type}

  if (!missing(min) && !is.numeric(min))
    stop("argument min should be numeric")

  if (!missing(max) && !is.numeric(max))
    stop("argument max should be numeric")

  if (type == "filename") {
    filename.mode <- match.arg(filename.mode)
  } else {
    if (!missing(filename.mode)) {
      warning("filename.mode is ignored when type != 'filename'")
    }
  }

  # check that type matches default
  if (!missing(default)) {
    defaultNotInItems <- !missing(items) && !(default %in% items)
    if (((type=="int" || type=="range-int" || type=="double" || type=="range-double") && !is.numeric(default) && defaultNotInItems) ||
          ( type=="boolean" && !is.logical(default)) ||
          ( type=="string" && !is.character(default)))
      stop ("default value does not match the type")

    if (type == "color") {
      default <- paste('#', paste(as.hexmode(col2rgb(default)), collapse=""), '', sep='')
      #if (length(default) > 1) default <- apply(as.character(as.hexmode(col2rgb(default))), 2, paste, collapse="")
    }
  }

  # try to guess a debug value from the other properties
  if (missing(debugvalue)) {
    if (!missing(default)) {
      debugvalue <- default
    } else if (!missing(items)) {
      debugvalue <- items[1]
    } else if (type == "column" && !is.null(visr.input)) {
      debugvalue <- colnames(visr.input)[1]
    } else if (type == "column-numerical" && !is.null(visr.input)) {
      debugvalue <- colnames(visr.input)[which(sapply(visr.input, is.numeric))][1]
    } else if (type == "multi-column" && !is.null(visr.input)) {
      debugvalue <- colnames(visr.input)
    } else if (type == "multi-column-numerical" && !is.null(visr.input)) {
      debugvalue <- colnames(visr.input)[which(sapply(visr.input, is.numeric))]
    } else if (type == "multi-color") {
      debugvalue <- "#d73027,#fc8d59,#fee090,#ffffbf,#e0f3f8,#91bfdb,#4575b4"
    } else if (type == "int") {
      debugvalue = if (!missing(min)) min else 0L
    } else if (type == "double") {
      debugvalue = if (!missing(min)) min else 0
    } else if (type == "range-int") {
      debugvalue = if (!missing(min) && !missing(max)) c(min, max) else c(0L, 10L)
    } else if (type == "range-double") {
      debugvalue = if (!missing(min) && !missing(max)) c(min, max) else c(0, 1)
    } else if (type == "boolean") {
      debugvalue = F
    } else {
      debugvalue = ''
    }
  }

  assign(paramname, debugvalue, envir = .GlobalEnv)
  if (!is.null(debugvalue))
    visr.params[[paramname]] <<- debugvalue
  else
    visr.params[paramname] <<- list(NULL)

  default.ischar <- !missing(default) && is.character(default)

  if (type == "boolean" && !missing(default))
    default <- tolower(default) # TRUE -> "true"

  if (!missing(active.condition))
    eval(parse(text=active.condition)) # evaluate to make sure it is valid

  properties <- c(
    if (!missing(label))    {paste('"label": ',   gsub('"', '\\\\"', label),    '', sep='"')} else {NULL},
    if (!missing(info))     {paste('"info": ',    gsub('"', '\\\\"', gsub('\n', '\\\\n', info)),     '', sep='"')} else {NULL},
    if (!missing(type))     {paste('"type": ',    type,     '', sep='"')} else {NULL},
    if (!missing(default))  {
      if (length(default) == 1)
        paste('"default": ', default , '', sep = ifelse(default.ischar, '"', ''))
      else
        paste0('"default": ', paste0('[', paste('', default,       '', sep = ifelse(default.ischar, '"', ''), collapse=", ") ,"]")) # range-int, range-double or multi-string
    } else {NULL},
    if (!missing(min))      {paste('"min": ',     min,          sep='')} else {NULL},
    if (!missing(max))      {paste('"max": ',     max,          sep='')} else {NULL},
    if (!missing(items))       {paste0('"items": ',       paste0('[', paste('"', items,       '"', sep="", collapse=", ") ,"]"))} else {NULL},
    if (!missing(item.labels)) {paste0('"item-labels": ', paste0('[', paste('"', item.labels, '"', sep="", collapse=", ") ,"]"))} else {NULL},
    if (!missing(filename.mode)) {paste('"filename.mode": ', filename.mode, '', sep='"')} else {NULL},
    if (!missing(active.condition))     {paste('"active-condition": ',    active.condition,     '', sep='"')} else {NULL},
    if (!missing(options))     {paste('"options": ',    options,     '', sep='"')} else {NULL}
  )

  #properties <- properties[which(!is.na(properties))]

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
  assign("visr.var.definedParam", TRUE, .GlobalEnv)
}

visr.app.param <- visr.param


#' Unit test
visr.internal.test.param <- function() {
  visr.app.start("test-app", info="A test app", debugdata = iris)
  visr.app.param("test-minimal")
  visr.app.param("test-auto-int", default = 2L)
  visr.app.param("test-auto-double", default = 0.5)
  visr.app.param("test-auto-bool", default = FALSE)
  visr.app.category("group2", "info for group2")
  #visr.param("test-mismatch", type="char", default=2)
  visr.app.param("test-color", label="foreground", info="foreground color", type="color", default = "yellow")
  visr.app.param("test-min-max", default=3, min=1, max=10)
  visr.app.param("test-filename", type="filename", filename.mode = "load")
  visr.app.param("test-items", items = c("i1","i2","i3"))
  visr.app.param("test-item-labels", items = c("i1","i2","i3"),
             item.labels = c("item 1","item 2","item 3"))
  visr.app.end(printjson=TRUE, filename="~/testapp.json")
  #cat(visr.var.appJSON)
}


#' Creates a list of strings from a comma or space separated string
#'
#' Creates a list of strings from a comma or space separated string
visr.util.split_comma_string <- function(comma_separated_string) {
  return (strsplit(comma_separated_string, split = "[ ]*[,| ]+[ ]*")[[1]])
}

#' Creates a gradient color scale
#'
#' Creates a gradient color scale to be used by ggplots from a string of color values
#'
#' @param colorsString    String of comma separted color values. e.g. "#d73027,#fc8d59,#fee090,#ffffbf,#e0f3f8,#91bfdb,#4575b4"
#'
#' @param label           The name of the scale. Used as axis or legend title.
#'                        If NULL, the default, the name of the scale is taken from the first mapping used for that aesthetic.
visr.util.scale_color_gradient <- function(colorsString, label = waiver()) {
  colorsVec <- visr.util.split_comma_string(colorsString)
  return (scale_colour_gradientn(colours = colorsVec, name = label))
}

#' Creates a fill color scale
#'
#' Creates a fill color scale to be used by ggplots from a string of color values
#'
#' @param colorsString    String of comma separted color values. e.g. "#d73027,#fc8d59,#fee090,#ffffbf,#e0f3f8,#91bfdb,#4575b4"
#'
#' @param label           The name of the scale. Used as axis or legend title.
#'                        If NULL, the default, the name of the scale is taken from the first mapping used for that aesthetic.
visr.util.scale_fill_gradientn <- function(colorsString, label = waiver()) {
  colorsVec <- visr.util.split_comma_string(colorsString)
  return (scale_fill_gradientn(colours = colorsVec, name = label))
}

# just a palette of 25 colors (https://stackoverflow.com/a/9568659/8371761)
visr.var.customPalette25 <- c(
  "black", "dodgerblue2","#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "gold1", "skyblue2","#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon","orchid1","deeppink1","blue1","steelblue4",
  "darkturquoise","green1","yellow4","yellow3",
  "darkorange4","brown"
)


todo_things_to_try <- function() {
    gc()
    gcinfo(TRUE)
    .dynLibs()
    library.dynam()
    getLoadedDLLs()

    warnings()
    tail(warnings(), 2)
    last.warning()

    #package: R.utils
    setTimeLimit()
    withTimeout(expr)
}

visr.library("assertthat", quiet = T)

#' Utility function
#' TODO
#'
visr.assert_that <- function(..., env = parent.frame(), msg = NULL) {
  if (visr.isGUI()) {
    res <- assertthat::see_if(..., env = env, msg = msg)
    if (res)
      return(TRUE)
    visr.message(attr(res, "msg"))
  } else {
    return (assertthat::assert_that(..., env = env, msg = msg))
  }
}

#' Utility function
#' TODO
#'
visr.assert_file_exists <- function(filename, parameterName) {
  visr.assert_that(!is.null(filename), msg = paste0("Parameter '", parameterName, "' is not specified"))
  visr.assert_that(filename != "", msg = paste0("Parameter '", parameterName, "' is not specified"))
  visr.assert_that(file.exists(filename), msg = paste0("Unable to find the path '", filename, "'"))
  return(TRUE)
}

