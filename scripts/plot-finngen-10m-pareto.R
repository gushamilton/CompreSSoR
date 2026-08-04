#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(dplyr)
})

root <- if (interactive()) "." else {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args)) args[[1]] else "."
}

summary_path <- file.path(
  root, "inst", "benchmarks", "finngen-10m-bp-20260804",
  "format-screen", "whole-file-10m-summary.csv"
)
output_path <- file.path(
  root, "inst", "benchmarks", "finngen-10m-bp-20260804",
  "format-screen", "pareto-10m-bp.png"
)

dat <- read_csv(summary_path, show_col_types = FALSE) |>
  filter(status == "OK") |>
  mutate(
    storage_ratio = as.numeric(source_bytes) / as.numeric(storage_bytes),
    read_seconds = as.numeric(read_median_seconds),
    family = factor(
      family,
      levels = c("Text", "VCF + Tabix", "Binary", "Framed binary",
                 "Parquet", "CompreSSoR")
    )
  )

stopifnot(
  nrow(dat) == 34L,
  all(dat$runs == 5L),
  all(dat$rows == 10000000L),
  all(is.finite(dat$storage_ratio)),
  all(is.finite(dat$read_seconds)),
  max(abs(dat$storage_ratio - dat$compression_ratio)) < 1e-10
)

# A point is Pareto-efficient when no other format is at least as compact and
# at least as fast, with one strict improvement.
dat <- dat |>
  rowwise() |>
  mutate(
    pareto = !any(
      dat$read_seconds <= read_seconds &
        dat$storage_ratio >= storage_ratio &
        (dat$read_seconds < read_seconds | dat$storage_ratio > storage_ratio)
    )
  ) |>
  ungroup() |>
  mutate(
    label = case_when(
      format_id == "pcodec_native" ~ "CompreSSoR Pcodec\n(native C++, 4 threads)",
      format_id == "parquet_q9" ~ "Parquet q9",
      format_id == "tsv_gzip" ~ "TSV.gz",
      TRUE ~ NA_character_
    )
  )

family_colours <- c(
  "Text" = "#4D4D4D",
  "VCF + Tabix" = "#D55E00",
  "Binary" = "#009E73",
  "Framed binary" = "#56B4E9",
  "Parquet" = "#0072B2",
  "CompreSSoR" = "#CC79A7"
)

frontier_line <- dat |> filter(pareto) |> arrange(read_seconds)

p <- ggplot(dat, aes(read_seconds, storage_ratio, colour = family, shape = family)) +
  geom_hline(yintercept = 1, colour = "grey55", linetype = "dashed", linewidth = 0.35) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(
    data = dat |> filter(!is.na(label)),
    aes(label = label),
    seed = 42, box.padding = 0.55, point.padding = 0.35,
    min.segment.length = 0, size = 3.6, show.legend = FALSE,
    colour = "grey15"
  ) +
  scale_x_log10(
    breaks = c(0.4, 0.5, 0.75, 1, 2, 5, 10, 20),
    labels = scales::label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.08, 0.14))
  ) +
  scale_y_continuous(
    breaks = scales::breaks_extended(n = 8),
    labels = function(x) paste0(format(x, trim = TRUE, digits = 3), "×"),
    expand = expansion(mult = c(0.06, 0.12))
  ) +
  scale_colour_manual(values = family_colours, drop = FALSE) +
  scale_shape_manual(
    values = c("Text" = 16, "VCF + Tabix" = 4, "Binary" = 15,
               "Framed binary" = 0, "Parquet" = 17, "CompreSSoR" = 18),
    drop = FALSE
  ) +
  labs(
    title = "10m FinnGen SNPs: keyed storage Pareto frontier",
    subtitle = "Five-run BluePebble medians; all formats include self-contained variant identity",
    x = "Whole-file read time (seconds; log scale)",
    y = "Storage compression relative to TSV.gz",
    colour = "Format family",
    shape = "Format family",
    caption = "Source: finngen-10m-bp-20260804. Shared reference excluded; TSV.gz baseline = 1×."
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey35"),
    plot.caption = element_text(colour = "grey40", hjust = 0),
    plot.margin = margin(12, 18, 10, 12)
  )

if (nrow(frontier_line) > 1L) {
  p <- p + geom_line(
    data = frontier_line,
    linewidth = 0.55, colour = "grey30", inherit.aes = FALSE,
    aes(read_seconds, storage_ratio)
  )
}

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
ggsave(output_path, p, width = 11, height = 7.5, units = "in", dpi = 220, bg = "white")

message("wrote ", output_path)
message("formats=", nrow(dat), "; frontier=", sum(dat$pareto),
        "; pcodec_ratio=", sprintf("%.3f", dat$storage_ratio[dat$format_id == "pcodec_native"]),
        "; pcodec_read_s=", sprintf("%.3f", dat$read_seconds[dat$format_id == "pcodec_native"]))
