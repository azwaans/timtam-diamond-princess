#'
#' Description
#' ===========
#'
#' Generate a plot of the R0 and prevalence estimates based on the data read
#' from the files specified below. Note that there are also a couple of dates
#' hard coded there to make some things a bit simpler.
#'
#' There are three sections:
#' - Configuration :: where some files and dates are hard coded
#' - Helper functions :: where helper functions are defined.
#' - Main figure generation :: where the final result is made.
#'
#' Usage
#' =====
#'
#' $ Rscript ./R/postprocessing-part-4.R
#'
library(ggplot2)
library(xml2)
library(dplyr)
library(cowplot)
library(gridExtra)
library(RColorBrewer)
palette <- brewer.pal(5, "Dark2")[1:3]
palette_green <- "#1B9E77"
palette_orange <- "#D95F02"
palette_purple <- "#7570B3"

## Configuration
## =============

beast_xml <- "./xml/timtam-timeseries-stage-2-with-removal-prevalence-HKY-weekly-histories.xml" # nolint
hist_hdi_file <- "./out/prevalence-estimate-HKY-weekly-histories.csv"
param_log_file <- "./out/log-files/timtam-timeseries-stage-2-with-removal-prevalence-HKY-weekly-histories.1.log" # nolint
date_of_last_seq <- "2020-02-17"
date_range <- as.Date(c("2020-01-20", "2020-02-27"))
output_file <- "./out/manuscript/timeline.png"

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

## ====================
## End helper functions

## Main figure generation
## ======================

hist_hdi_df <-
  read.csv(hist_hdi_file) |>
  mutate(history_dates = as.Date(-history_times, origin = date_of_last_seq))

#' These estimates where obtained from figures in the manuscript after
#' extracting the values with WebPlotDigitizer.
#'
#' TODO I'm not sure if these values for the prevalence account for
#' the LTT value. If they did not it would go some way to explaining
#' the difference in the values.

andreoletti2022estimates_prev_df <-
  data.frame(
    history_dates = hist_hdi_df$history_dates,
    HDIup = c(3.002684,  7.049955, 10.100358, 18.493290, 53.674918,  4.000000),
    HDIlow = c(1.002684,  1.049955,  3.100358,  8.493290, 42.674918,  1.000000),
    mean = c(2.002684,  3.049955,  7.100358, 13.493290, 48.674918,  2.000000)
  )

plot_prevalence <-
  ggplot() +
  geom_linerange(hist_hdi_df,
                 mapping = aes(x = history_dates, ymin = HDIlow, ymax = HDIup),
                 colour = palette_green) +
  geom_point(hist_hdi_df,
             mapping = aes(x = history_dates, y = mean),
             colour = palette_green) +
  geom_linerange(andreoletti2022estimates_prev_df,
                 mapping = aes(x = history_dates, ymin = HDIlow, ymax = HDIup),
                 colour = palette_purple) +
  geom_point(andreoletti2022estimates_prev_df,
             mapping = aes(x = history_dates, y = mean),
             colour = palette_purple) +
  scale_y_continuous(name = "Prevalence") +
  scale_x_date(limits = date_range) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0)
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
plot_df$point_est <- colMeans(post_r0_df)

#' Because we don't want the step function to end too early, we duplicate the
#' last value so it continues on until the end of the study period.
plot_df <- rbind(plot_df, tail(plot_df, 1))
plot_df$date <- r0_change_dates
plot_ribbonstep_df <-
  make_ribbonstep_df(xs = plot_df$date,
                     ymins = head(plot_df$CI_low, -1),
                     ymaxs = head(plot_df$CI_high, -1))


## ===================================================================
## Create some bogus data to test the visualisation with the first
## data frame provides a point estimate and the second a range.
## ===================================================================
plot_df <-
  data.frame(
    xs = r0_change_dates,
    ys = 0.5 * (c(1.02, 0.97, 0.891, 0.891) +
                c(7.54, 1.48, 1.11, 1.11))
  )
plot_rbn_df <-
  make_ribbonstep_df(
  xs = r0_change_dates,
  ymins = c(1.02, 0.97, 0.891) + rexp(n = 3, rate = 10),
  ymaxs = c(7.54, 1.48, 1.11) - rexp(n = 3, rate = 10)
  )
## ===================================================================

## True values of estimates have been hard-coded from the manuscripts.
vaughan2020estimates_df <- data.frame(
  xs = as.Date(c("2020-01-20", "2020-02-03", "2020-02-27")),
  ys = c(5.361, 1.84, 1.84)
)

vaughan2020estimates_rbn_df <-
  make_ribbonstep_df(
    xs = as.Date(c("2020-01-20", "2020-02-03", "2020-02-27")),
    ymins = c(7.6752, 3.2938),
    ymaxs = c(2.1708, 1.5002)
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
              size = 0.5) +
    ## make sure that the text is aligned to the right of the specified point
    geom_text(data = data.frame(x = as.Date("2020-02-20"), y = 8.5),
              mapping = aes(x = x, y = y, label = label),
              size = 4,
              hjust = hjust,
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

plot_r0_vaughan2020estimates <-
  make_r0_plot(vaughan2020estimates_rbn_df,
               vaughan2020estimates_df,
               palette_orange,
               "Vaughan et al. (2020)",
               hjust = 1, y_text_p = FALSE)

plot_r0_andreoletti2022estimates <-
  make_r0_plot(andreoletti2022estimates_rbn_df,
               andreoletti2022estimates_df,
               palette_purple,
               "Andreoletti et al. (2022)",
               hjust = 0.9, y_text_p = FALSE)

plot_r0_timtam <-
  make_r0_plot(plot_rbn_df,
               plot_df,
               palette_green,
               "Timtam",
               hjust = 2.8) +
  scale_y_continuous(limits = c(0, 9),
                     breaks = 1:8,
                     name = "Reproduction number") +
  theme(axis.title.y = element_text(size = 10))

## We use gridExtra::grid.arrange() because it is more flexible than
## cowplot::plot_grid().
plot_r0 <-
  grid.arrange(
    plot_r0_timtam,
    plot_r0_andreoletti2022estimates,
    plot_r0_vaughan2020estimates,
    ncol = 3, padding = unit(0, "mm")
  )

ggsave(filename = "tweaked-r0-plot.png",
       plot = plot_r0,
       height = 0.5 * 14.8,
       width = 21.0,
       units = "cm")


## combined_plot <- plot_grid(plot_r0, plot_prevalence, ncol = 1, align = "v", axis = "b")

## ggsave(filename = output_file,
##        plot = combined_plot,
##        height = 14.8,
##        width = 10.5,
##        units = "cm")
