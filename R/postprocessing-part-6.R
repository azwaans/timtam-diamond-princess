library(dplyr)
library(ggplot2)
library(ggpattern)
library(reshape2)
library(cowplot)
library(gridExtra)

TEXT_SIZE <- 12
TITLE_SIZE <- 15
LINE_SIZE <- 1.0
POINT_SIZE <- 4
ANNOTATION_SIZE <- 5

FIGURE_WIDTH <- 160 # mm

disaster_data_csv <- "out/disaster-data.csv"
sample_times_csv <- "out/diamond-times.csv"
r0_plot_components_rds <- "out/manuscript/r0-estimates-components.rds"
prev_plot_components_rds <- "out/manuscript/prevalence-estimates-components.rds"

plot_panels_ab <- "out/manuscript/diamond-figure-panels-ab.svg"
plot_panel_c <- "out/manuscript/diamond-figure-panel-c.svg"

if (!file.exists(disaster_data_csv)) {
  stop("The file ", disaster_data_csv, " does not exist.")
}
if (!file.exists(sample_times_csv)) {
  stop("The file ", sample_times_csv, " does not exist.")
}
if (!file.exists(r0_plot_components_rds)) {
  stop("The file ", r0_plot_components_rds, " does not exist.")
}

disaster_df <-
  disaster_data_csv |>
  read.csv() |>
  mutate(
    type = "Cases",
    time = floor(time)
  )
sample_df <-
  sample_times_csv |>
  read.csv() |>
  select(-num) |>
  mutate(
    count = 1,
    time = floor(time)
  ) |>
  group_by(time) |>
  summarise(count = sum(count)) |>
  mutate(
    type = "Sequences",
    time = 28 - time
  ) |>
  as.data.frame()
data_df <-
  bind_rows(disaster_df, sample_df) |>
  mutate(date = as.Date("2020-02-17") - time)


interval_df <- data.frame(
  start_date = as.Date(c("2020-01-20", "2020-02-04", "2020-02-11",
                         "2020-02-03", "2020-02-04", "2020-02-15")),
  end_date = as.Date(c("2020-02-03", "2020-02-11", "2020-02-27",
                       "2020-02-27", "2020-02-19", "2020-02-17")),
  interval = c("Cruise", "Sympt. testing", "Increased testing",
               "Harbour quarantine", "Cabin quarantine", "Sequencing"),
  y = 15 * c(1, 2, 3, 4, 5, 2) + 80
)

date_range <- as.Date(c("2020-01-20", "2020-02-27"))

data_gg <-
  ggplot() +
  geom_segment(
    data = interval_df,
    mapping = aes(x = start_date, xend = end_date, y = y, yend = y)
  ) +
  geom_point(
    data = melt(select(interval_df, -interval), id.vars = "y"),
    mapping = aes(x = value, y = y)
  ) +
  geom_text(
    data = interval_df,
    mapping = aes(x = start_date, y = y, label = interval),
    vjust = -1, hjust = 0
  ) +
  geom_col_pattern(
    data = data_df,
    mapping = aes(x = date, y = count,
                  pattern = type),
    fill = "white",
    pattern_spacing = 0.015, pattern_angle = 45
  ) +
  scale_pattern_manual(
    values = c("stripe", "crosshatch")
  ) +
  scale_y_continuous(
    breaks = seq(from = 0, to = 120, by = 40),
    minor_breaks = seq(from = 0, to = 120, by = 20),
    limits = c(0, max(interval_df$y) + 10),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_x_date(limits = date_range) +
  labs(x = NULL, y = "Daily count", pattern = NULL, pattern_angle = NULL) +
  theme_bw() +
  theme(
    legend.position = c(0.2, 0.3),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = TITLE_SIZE),
    axis.text = element_text(size = TEXT_SIZE)
  )

prev_plot_components <- readRDS(prev_plot_components_rds)
palette_purple <- prev_plot_components$purple
palette_green <- prev_plot_components$green
prev_plot_df <- prev_plot_components$data_df

prev_plot_gg <-
  ggplot() +
    geom_label(
      data = data.frame(x = rep(as.Date("2020-01-21"), 2),
                        y = c(80,95),
                        label = c("\u25b2 Andréoletti et al. (2022)",
                                  "\u25cf  Timtam")),
      mapping = aes(x = x, y = y, label = label),
      size = ANNOTATION_SIZE,
      hjust = 0, # so text is left aligned
      fill = "white",
      color = c(palette_purple, palette_green)
    ) +
    geom_linerange(
      prev_plot_df,
      mapping = aes(x = history_dates, ymin = HDIlow, ymax = HDIup,
                    colour = est_source),
      linewidth = LINE_SIZE
    ) +
    geom_point(
      prev_plot_df,
      mapping = aes(x = history_dates, y = mean,
                    colour = est_source, shape = est_source),
      size = POINT_SIZE,
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
      axis.title.y = element_text(size = TITLE_SIZE),
      axis.text = element_text(size = TEXT_SIZE)
    )


timed_gg <- plot_grid(data_gg, prev_plot_gg, align = "vh", ncol = 1, rel_heights = c(1, 0.7))
ggsave(plot_panels_ab, timed_gg,
       width = FIGURE_WIDTH, height = 160,
       units = "mm",
       dpi = 300)




r0_plot_components <- readRDS(r0_plot_components_rds)
r0_comp_a <- r0_plot_components[[1]] +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = TITLE_SIZE),
    axis.text = element_text(size = TEXT_SIZE),
    axis.text.x = element_text(size = TEXT_SIZE, angle = 90)
  )
r0_comp_b <- r0_plot_components[[2]] +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(size = TEXT_SIZE),
    axis.text.x = element_text(size = TEXT_SIZE, angle = 90)
  )
r0_comp_c <- r0_plot_components[[3]] +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(size = TEXT_SIZE),
    axis.text.x = element_text(size = TEXT_SIZE, angle = 90)
  )

r0_gg <- plot_grid(
  r0_comp_a,
  r0_comp_b,
  r0_comp_c,
  align = "t", ncol = 3, rel_heights = c(1, 1, 1)
)

ggsave(plot_panel_c, r0_gg,
       width = FIGURE_WIDTH, height = 65,
       units = "mm",
       dpi = 300)
