error_message <- ""

visr.isGUI<-function(){
  return (exists("visr.var.isGUI") && visr.var.isGUI)
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

# utility function to open data tables with corrected column names for debugging within R studio.
visr.readInputTable <-function(file) {
  input.table <<- read.csv(file, sep = "\t", check.names = F)
  colnames(input.table) <- make.names(gsub("[^a-zA-Z0-9_]", "_", colnames(input.table)))
  input_table <<- input.table
  return (input.table)
}

visr.setLogDir <- function(logDir) {
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

#.libPaths(c(visr.getHomeLibPath(), .libPaths()))
