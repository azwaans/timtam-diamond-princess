#'
#' Description
#' ===========
#'
#' Generate a plot of the R0 and prevalence estimates based on the data read
#' from the files specified below. Note that there are also a couple of dates
#' hard coded there to make some things a bit simpler.
#'
#' Usage
#' =====
#'
#' $ ./R/postprocessing-part-4.R -v \\
#' -x xml/timtam-timeseries-stage-2-with-removal-prevalence-HKY-weekly-histories.xml \\
#' -t out/log-files/timtam-timeseries-stage-2-with-removal-prevalence-HKY-weekly-histories-diamond.1.trees \\
#' -p out/log-files/timtam-timeseries-stage-2-with-removal-prevalence-HKY-weekly-histories.1.log \\
#' -r out/manuscript/r0-estimates.png \\
#' -n out/manuscript/prevalence-estimates.png \\
#'
#'
library(ggplot2)
library(xml2)
library(dplyr)
library(cowplot)
library(gridExtra)
library(RColorBrewer)
library(argparse)

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
  "-c",
  "--hdi-file",
  type = "character",
  help = "Filepath for the csv file containing the combined prevalence HDIs."
)
parser$add_argument(
  "-r",
  "--output-file-r",
  type = "character",
  help = "Filepath for the png to write Re plot to."
)
parser$add_argument(
  "-n",
  "--output-file-n",
  type = "character",
  help = "Filepath for the png to write prevalence plot to."
)



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

#' Make a data frame for a ribbon step plot
#'
#' @param xs x values
#' @param ymins y minimum values
#' @param ymaxs y maximum values
#'
#' @examples
#' plot_df <- make_ribbonstep_df(xs = c(1, 2, 3),
#'                              ymins = c(1, 2),
#'                             ymaxs = c(3, 5))
#' print(plot_df)
#' ggplot() +
#'  geom_ribbon(data = plot_df,
#'              mapping = aes(x = x, ymin = ymin, ymax = ymax),
#'              color = "red",
#'              fill = "red",
#'              alpha = 0.2) +
#'  theme_minimal()
make_ribbonstep_df <- function(xs, ymins, ymaxs) {
  stopifnot(length(xs) == length(ymins) + 1)
  stopifnot(length(ymins) == length(ymaxs))
  middle_xs <- xs[-c(1, length(xs))]
  data.frame(x = c(xs[1], rep(middle_xs, each = 2), xs[length(xs)]),
             ymin = rep(ymins, each = 2),
             ymax = rep(ymaxs, each = 2))
}

## Using a function ensures that the style of the figures is
## consistent across each set of estimates and makes it easier to
## construct the figures.
make_r0_plot <- function(rbn_df, est_df, colour, label, hjust, y_text_p = TRUE) {
  ggplot() +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_ribbon(data = rbn_df,
                mapping = aes(x = x, ymin = ymin, ymax = ymax),
                fill = colour, alpha = 0.3,
                linetype = 0) +
    geom_step(data = est_df,
              mapping = aes(x = xs, y = ys),
              color = colour,
              linewidth = 0.5) +
    ## make sure that the text is aligned to the right of the specified point
    geom_label(data = data.frame(x = as.Date("2020-01-21"), y = 8.5),
               mapping = aes(x = x, y = y, label = label),
               size = 4,
               hjust = 0,
               fill = "white",
               color = colour) +
    scale_y_continuous(limits = c(0, 9)) +
    theme_bw() +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 7),
          axis.text.x = element_text(angle = -40)) +
          ## axis.ticks.x = element_blank()) +
    {if (!y_text_p)
       theme(axis.text.y = element_blank(),
             axis.ticks.y = element_blank())
    }
}

## ====================
## End helper functions

main <- function(args) {
  beast_xml <- args$xml_file
  param_log_file <- args$parameter_file
  hist_hdi_file <- args$hdi_file
  trees_file <- args$trees_file
  output_png_r0 <- args$output_file_r
  output_png_prev <- args$output_file_n
  
  # set colors
  palette <- brewer.pal(5, "Dark2")[1:3]
  palette_green <- "#1B9E77"
  palette_orange <- "#D95F02"
  palette_purple <- "#7570B3"
  
  ## Check that the input files exist
  if (!file.exists(beast_xml)) {
    stop("The beast xml file does not exist.")
  }
  if (!file.exists(hist_hdi_file)) {
    stop("The history HDI file does not exist.")
  }
  if (!file.exists(param_log_file)) {
    stop("The parameter log file does not exist.")
  }
  
  date_of_last_seq <- "2020-02-17"
  date_range <- as.Date(c("2020-01-20", "2020-02-27"))
  
  hist_hdi_df <-
    read.csv(hist_hdi_file) |>
    mutate(history_dates = as.Date(-history_times, origin = date_of_last_seq),
           est_source = "timtam")
  
  ## Andréoletti
  andreoletti2022estimates_prev_df <-
    data.frame(
      history_dates = hist_hdi_df$history_dates,
      HDIup = c(3.002684,  7.049955, 10.100358, 18.493290, 53.674918,  4.000000),
      HDIlow = c(1.002684,  1.049955,  3.100358,  8.493290, 42.674918,  1.000000),
      mean = c(2.002684,  3.049955,  7.100358, 13.493290, 48.674918,  2.000000),
      est_source = "andreoletti2022estimates"
    )
  
  prev_plot_df <-
    bind_rows(hist_hdi_df, andreoletti2022estimates_prev_df)
  
  plot_prevalence <-
    ggplot() +
    geom_label(
      data = data.frame(x = rep(as.Date("2020-01-21"), 2),
                        y = c(80,95),
                        label = c("\u25b2 Andréoletti et al. (2022)",
                                  "\u25cf Timtam")),
      mapping = aes(x = x, y = y, label = label),
      size = 3.4,
      hjust = 0, # so text is left aligned
      fill = "white",
      color = c(palette_purple, palette_green)
    ) +
    geom_linerange(
      prev_plot_df,
      mapping = aes(x = history_dates, ymin = HDIlow, ymax = HDIup,
                    colour = est_source)
    ) +
    geom_point(
      prev_plot_df,
      mapping = aes(x = history_dates, y = mean,
                    colour = est_source, shape = est_source),
      size = 3,
    ) +
    scale_colour_manual(
      values = c(palette_purple, palette_green),
      labels = c("Andréoletti et al. (2022)", "Timtam")
    ) +
    scale_shape_manual(
      values = c(17, 16),
      labels = c("Andréoletti et al. (2022)", "Timtam")
    ) +
    scale_y_continuous(name = "Prevalence") +
    scale_x_date(limits = date_range) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text = element_text(size = 7),
      axis.text.x = element_text(angle = -40, hjust = 0)
    )
  
  #' Because we want to be able to compare things on calendar times we need to do
  #' a little bit of messy work to get the change times in terms of absolute
  #' dates.
  r0_change_dates <-
    read_xml(beast_xml) |>
    xml_find_first(xpath = "//parameter[@name=\'r0ChangeTimes\']") |>
    xml_text() |>
    strsplit(split = " ") |>
    unlist() |>
    as.numeric() |>
    (\(x) as.Date(-x, origin = date_of_last_seq))() |>
    c(date_range) |>
    sort()
  
  post_r0_df <- read_beast2_log(param_log_file) |> select(starts_with("TTR0"))
  plot_df <- bayestestR::hdi(post_r0_df) |> as.data.frame()
  if (any(is.na(plot_df))) {
    stop("Stopping because plot_df has NA values. Double check the log file used.")
  }
  plot_df$point_est <- colMeans(post_r0_df)
  
  #' Because we don't want the step function to end too early, we duplicate the
  #' last value so it continues on until the end of the study period.
  plot_df <- rbind(plot_df, tail(plot_df, 1))
  plot_df$xs <- r0_change_dates
  plot_df$ys <- plot_df$point_est
  plot_rbn_df <-
    make_ribbonstep_df(xs = plot_df$xs,
                       ymins = head(plot_df$CI_low, -1),
                       ymaxs = head(plot_df$CI_high, -1))
  
  vaughan2020estimates_df <- data.frame(
    xs = as.Date(c("2020-01-20", "2020-02-03", "2020-02-27")),
    ys = c(4.5171, 1.8713, 1.8713)
  )
  #' vValues obtained from Vaughan et al.'s log files, HDI on combined x 5 seeds. 
  vaughan2020estimates_rbn_df <-
    make_ribbonstep_df(
      xs = as.Date(c("2020-01-20", "2020-02-03", "2020-02-27")),
      ymaxs = c(6.2491, 2.2068),
      ymins = c(2.8986, 1.5495)
    )
  
  andreoletti2022estimates_df <- data.frame(
    xs = r0_change_dates,
    ys = c(4.01, 1.22, 0.996, 0.996)
  )
  andreoletti2022estimates_rbn_df <- make_ribbonstep_df(
    xs = r0_change_dates,
    ymins = c(1.02, 0.97, 0.891),
    ymaxs = c(7.54, 1.48, 1.11)
  )
  
  plot_r0_vaughan2020estimates <-
    make_r0_plot(vaughan2020estimates_rbn_df,
                 vaughan2020estimates_df,
                 palette_orange,
                 "Vaughan et al. (2024)",
                 hjust = 1, y_text_p = FALSE)
  
  plot_r0_andreoletti2022estimates <-
    make_r0_plot(andreoletti2022estimates_rbn_df,
                 andreoletti2022estimates_df,
                 palette_purple,
                 "Andréoletti et al. (2022)",
                 hjust = 0.9, y_text_p = FALSE)
  
  plot_r0_timtam <-
    make_r0_plot(plot_rbn_df,
                 plot_df,
                 palette_green,
                 "Timtam",
                 hjust = 2.85) +
    scale_y_continuous(limits = c(0, 9),
                       ## breaks = 1:8,
                       name = "Reproduction number") +
    theme(axis.title.y = element_text(size = 10))
  
  ## We use gridExtra::grid.arrange() because it is more flexible than
  ## cowplot::plot_grid(). The fudge factor here came from viewing
  ## measuring the panel widths on my screen to ensure that they are
  ## equal. The difference stems from the presence of the y-axis labels
  ## in the first panel.
  tmp_fudge_factor <- 1.14
  plot_r0 <-
    grid.arrange(
      plot_r0_timtam,
      plot_r0_andreoletti2022estimates,
      plot_r0_vaughan2020estimates,
      ncol = 3, padding = unit(0, "mm"),
      widths = c(tmp_fudge_factor, 1.0, 1.0)
    )
  
  saveRDS(list(
      plot_r0_timtam,
      plot_r0_andreoletti2022estimates,
    plot_r0_vaughan2020estimates),
    file = gsub(".png", "-components.rds", output_png_r0)
    )
  ggsave(filename = output_png_r0,
         plot = plot_r0,
         height = 7.4, width = 20,
         units = "cm")
  
  saveRDS(
    list(gg = plot_prevalence,
         data_df = prev_plot_df,
         purple = palette_purple,
         green = palette_green),
    file = gsub(".png", "-components.rds", output_png_prev))
  ggsave(filename = output_png_prev,
         plot = plot_prevalence,
         height = 7.4, width = 10.5, # A7
         units = "cm")
 
  
}

args <- parser$parse_args()
main(args)