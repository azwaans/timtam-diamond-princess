#!/usr/bin/env Rscript
#'
#' Combine LTT and History estimates
#' =================================
#'
#' This script generates a summary of the posterior samples for the prevalence
#' of infection through time taking into account the number of lineages in the
#' reconstructed tree at each point in time.
#'
#' Usage
#' =====
#'
#' Here is an example of how to use this script.
#'
#' $ ./R/combine-ltt-and-histories.R -v \\
#' -x xml/timtam-timeseries-stage-2-prevalence.xml \\
#' -t out/timeseries/log_files/timtam-timeseries-stage-2-with-removal-prevalence-HKY-weekly-histories_diamond.1 \\
#' -p out/timeseries/log_files/timtam-timeseries-stage-2-with-removal-prevalence-HKY-weekly-histories.1.log \\
#' -o out/timeseries/foo.csv
#'
#' Note
#' ====
#'
#' You may need to install beastio directly from GitHub
#'
#' > devtools::install_github("laduplessis/beastio")
#'


library(argparse)
library(ape)
library(beastio)
library(ggplot2)
library(xml2)

parser <- ArgumentParser()

parser$add_argument(
  "-v",
  "--verbose",
  action = "store_true",
  default = FALSE,
  help = "Verbose output (default FALSE)"
)
parser$add_argument(
  "-x",
  "--xml-file",
  type = "character",
  help = "Filepath for BEAST XML."
)
parser$add_argument(
  "-t",
  "--trees-file",
  type = "character",
  help = "Filepath for log file containing trees."
)
parser$add_argument(
  "-p",
  "--parameter-file",
  type = "character",
  help = "Filepath for log file containing the parameters."
)
parser$add_argument(
  "-o",
  "--output-file",
  type = "character",
  help = "Filepath for the CSV to write result to."
)

#' Read a single BEAST log file
#'
#' Read a single BEAST log file and return as a coda "mcmc" object.
#'
#' @param filename The name of the log file to read.
#' @param burnin Discard this proportion of samples at the start of the chain
#'        (if `burninAsSamples` == FALSE). Otherwise discard this many samples
#'        at the start of the chain.
#' @param maxsamples If > 0 stop after reading in this many lines
#'        (this option is only for testing and should generally not be used).
#' @param as.mcmc If FALSE then return an object of class "data.frame", else
#'        return an "mcmc" object
#' @param burninAsSamples if TRUE burnin is given as the number of samples,
#'        if FALSE burnin is a proportion (0 <= burnin < 1) of samples.
#'        (default = FALSE).
#'
#' @return An "mcmc" object (\code{\link[coda]{mcmc}}) or data frame
#'         (\code{\link{data.frame}}) object containing all of the parameters
#'         in the MCMC chain, with the burn-in discarded.
#'
#' @seealso \code{\link{readLog}}
#'
#' @examples
readSingleLog <- function(filename, burnin = 0.1, maxsamples = -1, as.mcmc = TRUE, burninAsSamples = FALSE) {
  if (!burninAsSamples && burnin > 1) {
    stop("Error: Burnin must be a proportion of samples between 0 and 1.", call. = FALSE)
  }

  if (burninAsSamples && burnin != round(burnin)) {
    stop("Error: Burnin must be an integer number of states.", call. = FALSE)
  }

  logfile <- read.table(filename, sep = "\t", header = TRUE, nrows = maxsamples)
  n <- nrow(logfile)

  if (!burninAsSamples) {
    burnSamples <- floor(burnin * n)
  } else {
    burnSamples <- burnin + 1
  }

  if (burnSamples >= n) {
    stop("Error: Discarding all samples in the log file.", call. = FALSE)
  }

  logfile <- logfile[burnSamples:n, ]

  if (as.mcmc == TRUE) {
    if (is.null(logfile$Sample)) {
      if (is.null(logfile$state)) {
        logfile$Sample <- as.numeric(rownames(logfile))
      } else {
        logfile$Sample <- logfile$state
        logfile$state <- NULL
      }
    }

    # Sometimes BEAST2 log files have an empty column at the end (lines end in tabs)
    # This needs to be removed before converting to mcmc object
    for (i in names(logfile)) {
      if (all(is.na(logfile[[i]]))) logfile[[i]] <- NULL
    }

    start <- logfile$Sample[1]
    thin <- logfile$Sample[2] - logfile$Sample[1]
    rownames(logfile) <- logfile$Sample
    logfile$Sample <- NULL

    return(coda::mcmc(logfile, start = start, thin = thin))
  } else {
    return(logfile)
  }
}

#' Returns the LTT value at the input time
#'
#' @param time a numeric time.
#' @param ltt_coords the coordinates from \code{ltt.plot.coords}
#'
find_ltt_value <- function(time, ltt_coords) {
  if (time < min(ltt_coords[, "time"]) || time > 0) {
    # the history time is outside the time range of the tree, return 0
    return(0)
  } else {
    # the history time is within the time range of the tree, the closest corresponding LTT value
    if (ltt_coords[[which(abs(ltt_coords[, "time"] - time) == min(abs(ltt_coords[, "time"] - time)))[1], "time"]] < time) {
      return(ltt_coords[[which(abs(ltt_coords[, "time"] - time) == min(abs(ltt_coords[, "time"] - time)))[1], "N"]])
    } else {
      return(ltt_coords[[which(abs(ltt_coords[, "time"] - time) == min(abs(ltt_coords[, "time"] - time)))[1] - 1, "N"]])
    }
  }
}


main <- function(args) {
  beast_xml <- args$xml_file
  parameter_trace_file <- args$parameter_file
  trees_file <- args$trees_file
  output_csv <- args$output_file

  ## It is easy to forget to unzip the tree file so there is a check here to
  ## provide some guidance if you end doing this.
  if (tools::file_ext(trees_file) != "trees") {
    warning("The posterior tree log file has a strange file extension.")
    if (tools::file_ext(trees_file) == "zip") {
      stop("Cannot read posterior trees from a zip file.")
    }
  }

  if (args$verbose) {
    if (!file.exists(beast_xml)) {
      stop("\n\n\tCould not find BEAST XML: ", beast_xml, "\n\n")
    }
    if (!file.exists(parameter_trace_file)) {
      stop("\n\n\tCould not find trace file: ", parameter_trace_file, "\n\n")
    }
    if (!file.exists(trees_file)) {
      stop("\n\n\tCould not find tree file: ", trees_file, "\n\n")
    }
  }

  beast_node <- read_xml(beast_xml)
  times_extract <- beast_node |>
    xml_find_first(xpath = "//parameter[@name=\"historyTimes\"]") |>
    xml_text() |>
    strsplit(" ") |>
    unlist() |>
    as.numeric()

  trees_trace <- readTreeLog(trees_file, burnin = 0.0)
  parameter_trace <- readSingleLog(parameter_trace_file, burnin = 0.0)


  ## We only want the trees that have a corresponding sample in the parameters
  ## log file. The following snippet will check the parameter samples are a
  ## superset of the tree samples and extract only the parameter samples. If
  ## both log files where sampled at the same iterations of the MCMC this should
  ## do nothing because they are already matched.
  trees_trace_sample_ixs <-
    trees_trace |>
    names() |>
    stringr::str_extract(pattern = "[0-9]+")
  parameter_trace_sample_ixs <- dimnames(parameter_trace)[[1]]
  ## This has been added to handle the case where the input vectors
  ## are not of suitable lengths to be compared with equality. If they
  ## are not of equal length then they cannot all be equal!
  if (length(trees_trace_sample_ixs) == length(parameter_trace_sample_ixs)) {
    ixs_identical_p <- all(trees_trace_sample_ixs == parameter_trace_sample_ixs)
  } else {
    ixs_identical_p <- FALSE
  }
  trees_for_each_sample_p <-
    all(
      is.element(
        el = trees_trace_sample_ixs,
        set = parameter_trace_sample_ixs
      )
    )
  geq_samples_trees_p <- nrow(parameter_trace) >= length(trees_trace)
  if (ixs_identical_p) {
    message("Sample ixs identical for log and tree files.")
  } else if (trees_for_each_sample_p && geq_samples_trees_p) {
    message("Subsampling the log to match the trees.")
    tmp <- parameter_trace |> as.data.frame()
    tmp_mask <- is.element(el = parameter_trace_sample_ixs,
                           set = trees_trace_sample_ixs)
    parameter_trace <- coda::as.mcmc(tmp[tmp_mask, ])
    rm(tmp)
    rm(tmp_mask)
  } else {
    stop("The parameter log file and the tree log file appear to be incompatible in that there are trees samples which do not have a matching parameter sample.") # nolint
  }

  if (nrow(parameter_trace) != length(trees_trace) ) {
    stop("Error: HistorySizes and Trees traces were not sampled at the same time during MCMC. They cannot be combined", call. = FALSE)
  }

  all_ltt <- c()
  for (i in seq_along(trees_trace)) {
    coords <- ltt.plot.coords(trees_trace[[i]])
    extracted_ltt_values <-
      sapply(times_extract, \(t) find_ltt_value(t, coords))
    all_ltt <- rbind(all_ltt, extracted_ltt_values)
  }


  for (i in seq.int(ncol(all_ltt))) {
    all_ltt[, i] <-
      all_ltt[, i] +
      as.vector(
        parameter_trace[seq_along(all_ltt[, i]),
                        paste0("TTHistorySizes.", i)]
      )
  }

  all_ltt <- data.frame(all_ltt)
  n <- nrow(all_ltt)
  burnin <- 0.1
  all_ltt <- all_ltt[floor(burnin*n):n,]

  if (nrow(all_ltt) < 20) {
    stop("There are very few posterior samples, the bayestestR::hdi function may return NAs. Try re-running this with a longer MCMC chain.") # nolint
  }

  hdi_intervals <-
    rbind(
      as.vector(sapply(all_ltt, mean)),
      sapply(all_ltt, \(x) as.numeric(bayestestR::hdi(x)[c(2, 3)]))
    )
  rownames(hdi_intervals) <- c("mean", "HDIlow", "HDIup")
  hdi_intervals <- data.frame(
    history_times = times_extract,
    t(hdi_intervals)
  )

  write.table(
    x = hdi_intervals,
    file = output_csv,
    sep = ",",
    row.names = FALSE
  )
}

if (!interactive()) {
  args <- parser$parse_args()
  main(args)
} else {
  message("Testing testing testing: Copying the given example!")
  test_args_1 <- list(
    output_file = "out/timeseries/foo.csv",
    parameter_file = "out/timeseries/log_files/timtam-timeseries-stage-2-prevalence.log", # nolint
    trees_file = "out/timeseries/log_files/timtam-timeseries-stage-2-prevalence-diamond.trees", # nolint
    verbose = TRUE,
    xml_file = "xml/timtam-timeseries-stage-2-prevalence.xml"
  )
  test_args_2 <- list(
    output_file = "out/prevalence-estimate-HKY-weekly-histories-TESTING.csv",
    parameter_file = "out/log-files/timtam-timeseries-stage-2-with-removal-prevalence-HKY-weekly-histories.1.log", # nolint
    trees_file = "out/log-files/timtam-timeseries-stage-2-with-removal-prevalence-HKY-weekly-histories-diamond.1.trees", # nolint
    verbose = TRUE,
    xml_file = "xml/timtam-timeseries-stage-2-with-removal-prevalence-HKY-weekly-histories.xml" # nolint
  )
  main(test_args_2)
}
