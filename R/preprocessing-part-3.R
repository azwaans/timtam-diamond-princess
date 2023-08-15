library(dplyr)
library(ggplot2)
library(ggpattern)
library(reshape2)

disaster_data_csv <- "out/disaster-data.csv"
sample_times_csv <- "out/diamond-times.csv"

if (!file.exists(disaster_data_csv)) {
  stop("The file ", disaster_data_csv, " does not exist.")
}
if (!file.exists(sample_times_csv)) {
  stop("The file ", sample_times_csv, " does not exist.")
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
  labs(x = NULL, y = NULL, pattern = NULL, pattern_angle = NULL) +
  theme_bw() +
  theme(
    legend.position = c(0.2, 0.3)
  )

ggsave(
  filename = "out/manuscript/data-plot.png",
  plot = data_gg,
  height = 0.75 * 14.8, width = 21.0,
  units = "cm"
)
