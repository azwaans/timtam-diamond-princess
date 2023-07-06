library(ape)
library(phangorn)
library(stringr)
library(ggplot2)
set.seed(1)

main <- function() {
  diamond_fasta <- "data/diamond.fasta"
  diamond_csv <- "out/diamond-times.csv"

  if (!file.exists(diamond_fasta)) {
    stop("The required FASTA file ", diamond_fasta, " is missing.")
  }
  seq_times <-
    diamond_fasta |>
    read.FASTA() |>
    names() |>
    str_extract("[.0-9]+$") |>
    as.numeric() |>
    sort() |>
    (\(x) data.frame(time = x, num = seq_along(x)))()

  write.table(
    x = seq_times,
    file = diamond_csv,
    sep = ",",
    row.names = FALSE
  )

  seq_gg <-
    ggplot() +
    geom_point(
      data = seq_times,
      mapping = aes(time, num)
    ) +
    geom_vline(xintercept = 25:30) +
    labs(x = "Day number", y = NULL, title = "Sequences") +
    theme_bw()

  ggsave(
    filename = "out/sequence-times.png",
    plot = seq_gg,
    height = 10.5, width = 14.8,
    ## A5 height = 14.8, width = 21.0,
    ## A6 height = 10.5, width = 14.8,
    ## A7 height = 7.4, width = 10.5,
    units = "cm"
  )
}

main()
