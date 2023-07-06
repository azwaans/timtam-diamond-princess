library(stringr)

times <- seq(from = 13.4655172, to = -9.5344828, by = -1)
counts <- c(10, 10, 41, 3, 6, 65,
  20, 19, 44, 41, 26, 50, 77,
  60, 79, 13, 0, 0, 57, 0,
  0, 14, 0, 0
)

write.table(x = data.frame(time = times, count = counts),
            file = "out/disaster-data.csv",
            sep = ",",
            row.names = FALSE)
