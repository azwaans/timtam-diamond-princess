#'
#' Description
#' ===========
#'
#' Check some diagnostics on the MCMC samples to ensure that nothing
#' has gone obviously wrong.
#'
#' Usage
#' =====
#'
#' $ Rscript ./R/postprocessing-part-5.R
#'
library(ggplot2)
library(xml2)
library(dplyr)
library(cowplot)
library(RColorBrewer)
library(coda)
palette <- brewer.pal(5, "Dark2")[1:3]
palette_green <- "#1B9E77"
palette_orange <- "#D95F02"
palette_purple <- "#7570B3"

## Configuration
## =============

param_log_file <- "./out/log-files/timtam-timeseries-stage-2-with-removal-prevalence-HKY-weekly-histories_fasta.1.log" # nolint
output_file <- "./out/manuscript/diagnostics.txt"

## Check that the input files exist
if (!file.exists(param_log_file)) {
  stop("The parameter log file does not exist.")
}

## Helper functions
## ================

#' Read a BEAST2 log file into a data frame.
#'
#' @param filename is the path to the log file.
#' @param burn is the number to remove from the start.
#' @param take_last is the number to take from the end.
#'
#' @return data frame containing the samples.
#'
read_beast2_log <- function(filename, burn = 0, take_last = NA) {
  y <- read.csv(filename, sep = "\t", comment.char = "#")
  if (is.na(take_last) && burn >= 0) {
    return(tail(y, nrow(y) - burn))
  } else if (!is.na(take_last) && burn == 0) {
    return(tail(y, take_last))
  } else {
    stop("Unsupported arguments given to read_beast2_log.")
  }
}

## ====================
## End helper functions

## Main
## ====
## -TTPropPsi.1,-TTPropPsi.2,-TTPropPsi.3,
log_mcmc <- read_beast2_log(param_log_file) |>
  select(-Sample, -matches("TTPropPsi"),  -matches("TTPropTS")) |>
  as.mcmc()

## Log the minimum effective sample file to the output file using the sink command
sink(output_file, append = FALSE, split = FALSE)
print("Minimum effective sample sizes:")
print(min(effectiveSize(log_mcmc)))
sink()
