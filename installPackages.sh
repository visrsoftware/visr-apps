#!/bin/bash
CURR_DIR=`pwd`
INSTALL_SCRIPT=installPackages.R
rm -rf $INSTALL_SCRIPT
echo 'sessionInfo()' >> $INSTALL_SCRIPT
echo 'while (is.null(try(remove.packages("BiocInstaller"), silent=T))) {}' >> $INSTALL_SCRIPT
echo 'source("visrutils.R")' >> $INSTALL_SCRIPT
if [ "$1" != "-su" ]; then
  echo 'dir.create(visr.getHomeLibPath(), recursive=T)' >> $INSTALL_SCRIPT 
  echo '.libPaths(c(visr.getHomeLibPath(), .libPaths()))' >> $INSTALL_SCRIPT 
  echo 'while (is.null(try(remove.packages("BiocInstaller"), silent=T))) {}' >> $INSTALL_SCRIPT
fi
echo 'source("http://bioconductor.org/biocLite.R")' >> $INSTALL_SCRIPT
echo 'try(biocValid(fix=T))' >> $INSTALL_SCRIPT

find `pwd` -name "*.R" -print0 | xargs -0 egrep -h "(\bvisr.biocLite\b\(\"|\bvisr.library\b\(\"|\bvisr.libraryURL\b\(\")" | awk '!seen[$0] {print "try(" $0 ")"} {++seen[$0]}' >> $INSTALL_SCRIPT

echo 'try(biocValid(fix=T))' >> $INSTALL_SCRIPT

Rscript --version
Rscript --vanilla --verbose $INSTALL_SCRIPT
