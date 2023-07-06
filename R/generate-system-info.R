
## ===================================================================
## The following script can be used to get a list of all the packages
## used by scripts in this directory.

## #!/bin/bash

## # Loop through all files with the .R extension in the current directory
## for file in *.R; do
##   # Grep lines containing the library() call and extract only the alphanumeric part within the parentheses
##   grep -oP 'library\(\K[[:alnum:]]+(?=\))' "$file"

##   # Grep lines containing alphanumeric strings (including underscores) ending with :: and extract only the alphanumeric part before the ::
##   grep -oP '[[:alnum:]_]+(?=::)' "$file"
## done | sort | uniq
## ===================================================================

#' Note
#' ====
#'
#' You may need to install beastio directly from GitHub
#'
#' > devtools::install_github("laduplessis/beastio")
#'

library(ape)
library(argparse)
library(bayestestR)
library(beastio)
library(coda)
library(cowplot)
library(devtools)
library(dplyr)
library(ggplot2)
library(magrittr)
library(phangorn)
library(purrr)
library(RColorBrewer)
library(reshape2)
library(stringr)
library(tools)
library(xml2)
library(xtable)

## Sink the following output to a file "session-info.txt"
sink("session-info.txt")
print(sessionInfo())
sink()
