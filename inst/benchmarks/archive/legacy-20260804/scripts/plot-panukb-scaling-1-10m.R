#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1) args[[1]] else "outputs/benchmarks/panukb-scaling-1-10m.csv"
output <- if (length(args) >= 2) args[[2]] else "outputs/benchmarks/panukb-scaling-1-10m.png"

dat <- read.csv(input, check.names = FALSE)
dat$rows_m <- dat$rows / 1e6
dat$format <- factor(
  dat$format,
  levels = c("TSV.gz", "Parquet q9", "CompreSSoR Pcodec")
)
min_rows_m <- min(dat$rows_m)
max_rows_m <- max(dat$rows_m)
fine_scale <- max_rows_m / min_rows_m >= 50
row_label <- function(x) {
  ifelse(x < 1, paste0(format(round(x * 1000), trim = TRUE), "k"), paste0(x, "m"))
}

storage <- data.frame(
  rows_m = dat$rows_m,
  format = dat$format,
  metric = "Storage (MB)",
  value = dat$storage_bytes / 1e6
)
ratio <- data.frame(
  rows_m = dat$rows_m,
  format = dat$format,
  metric = "Compression ratio vs TSV.gz (x)",
  value = dat$relative_to_tsvgz
)
read_time <- data.frame(
  rows_m = dat$rows_m,
  format = dat$format,
  metric = "Whole-file read time (seconds)",
  value = dat$read_median_seconds
)
plot_dat <- rbind(storage, ratio, read_time)
plot_dat$metric <- factor(
  plot_dat$metric,
  levels = c(
    "Storage (MB)",
    "Compression ratio vs TSV.gz (x)",
    "Whole-file read time (seconds)"
  )
)

x_scale <- if (fine_scale) {
  scale_x_log10(
    breaks = c(0.01, 0.1, 1),
    labels = c("10k", "100k", "1m")
  )
} else {
  scale_x_continuous(breaks = seq(1, max_rows_m, by = 1))
}

p <- ggplot(plot_dat, aes(x = rows_m, y = value, colour = format, group = format)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.1) +
  facet_grid(metric ~ ., scales = "free_y") +
  x_scale +
  scale_colour_manual(
    values = c(
      "TSV.gz" = "#3B6FB6",
      "Parquet q9" = "#D27A2C",
      "CompreSSoR Pcodec" = "#7A4FA3"
    ),
    drop = FALSE
  ) +
  labs(
    title = "Storage and read-speed scaling on Pan-UKB-derived GWAS rows",
    subtitle = paste0(
      row_label(min_rows_m), "–", row_label(max_rows_m),
      " valid autosomal biallelic SNP rows; three complete reads per point; BP compute node"
    ),
    x = if (fine_scale) "Rows (millions; log scale)" else "Rows (millions)",
    y = NULL,
    colour = "Format"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text.y = element_text(angle = 0, hjust = 0),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "#555555")
  )

ggsave(output, p, width = 11, height = 10, dpi = 180, bg = "white")
message(output)
