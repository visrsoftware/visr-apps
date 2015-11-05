#!/bin/bash
CURR_DIR=`pwd`
INSTALL_SCRIPT=installPackages.R
rm -rf $INSTALL_SCRIPT
echo 'source("visrutils.R")' > $INSTALL_SCRIPT
if [ "$1" != "-su" ]; then
  echo '.libPaths(c(visr.getHomeLibPath(), .libPaths()))' >> $INSTALL_SCRIPT 
fi

find `pwd` -name "*.R" -print0 | xargs -0 egrep -h "(\bvisr.biocLite\b\(\"|\bvisr.library\b\(\"|\bvisr.libraryURL\b\(\")" | awk '!seen[$0] {print "try(" $0 ")"} {++seen[$0]}' >> $INSTALL_SCRIPT

Rscript $INSTALL_SCRIPT
