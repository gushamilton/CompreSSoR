#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

summary_path <- Sys.getenv(
  "COMPRESSOR_PARETO_SUMMARY",
  unset = "inst/benchmarks/pareto-chr1-summary.csv"
)
output_base <- Sys.getenv(
  "COMPRESSOR_PARETO_PLOT",
  unset = "inst/figures/compressor-pareto"
)
if (!file.exists(summary_path)) stop("missing benchmark summary: ", summary_path)

dat <- fread(summary_path)
key_formats <- c(
  "TSV gzip", "TSV uncompressed", "Parquet self-contained exact",
  "VCF bgzip + Tabix", "CompreSSoR Pcodec self-contained"
)
dat <- dat[format %in% key_formats]
if (!nrow(dat)) stop("benchmark summary contains no variant-key formats")
dat[, family := fifelse(grepl("CompreSSoR", format), "CompreSSoR",
                 fifelse(grepl("Parquet", format), "Parquet",
                 fifelse(grepl("VCF", format), "VCF + Tabix",
                 fifelse(grepl("TSV", format), "Text", "Binary"))))]
dat[, frontier := vapply(seq_len(.N), function(i) {
  !any(dat$median_seconds <= dat$median_seconds[i] &
       dat$compression_ratio >= dat$compression_ratio[i] &
       (dat$median_seconds < dat$median_seconds[i] |
        dat$compression_ratio > dat$compression_ratio[i]))
}, logical(1))]
dat[, label := fifelse(frontier | format == "TSV gzip", format, NA_character_)]
dat[format == "CompreSSoR Pcodec self-contained", label :=
      "CompreSSoR Pcodec self-contained\n(identity included)"]
labels <- dat[!is.na(label)]
frontier_data <- dat[frontier == TRUE][order(median_seconds)]

family_levels <- c("CompreSSoR", "Parquet", "Text", "VCF + Tabix")
dat[, family := factor(family, levels = family_levels)]
labels[, family := factor(family, levels = family_levels)]

p <- ggplot(dat, aes(median_seconds, compression_ratio, colour = family,
                      shape = family)) +
  geom_hline(yintercept = 1, colour = "grey60", linetype = "dashed",
             linewidth = 0.45) +
  geom_path(data = frontier_data,
            aes(median_seconds, compression_ratio),
            inherit.aes = FALSE,
            colour = "grey25", linewidth = 0.75) +
  geom_point(aes(size = ifelse(frontier, 3.4, 2.4)), alpha = 0.9) +
  geom_text_repel(data = labels, aes(label = label), seed = 20260804,
                  size = 3.3, lineheight = 0.9, box.padding = 0.5,
                  point.padding = 0.25, force = 1.7, force_pull = 0.3,
                  min.segment.length = 0, max.overlaps = Inf,
                  show.legend = FALSE) +
  scale_size_identity() +
  scale_x_log10(breaks = c(0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2),
                labels = label_number(accuracy = 0.01),
                expand = expansion(mult = c(0.04, 0.25))) +
  scale_y_continuous(labels = function(x) paste0(number(x, accuracy = 0.1), "×"),
                     expand = expansion(mult = c(0.04, 0.08))) +
  scale_colour_manual(values = c(CompreSSoR = "#8e44ad", Parquet = "#1769aa",
                                  Text = "#6b7280",
                                  `VCF + Tabix` = "#dc2626"), drop = FALSE) +
  scale_shape_manual(values = c(CompreSSoR = 17, Parquet = 16,
                                 Text = 18, `VCF + Tabix` = 4),
                     drop = FALSE) +
  labs(
    title = "FinnGen chr1 variant-key compression/access Pareto frontier",
    subtitle = paste0(format(dat$rows[1], big.mark = ","),
                      " SNVs; five-run access medians; same input for every point"),
    x = "Whole-file access time (seconds; log scale)",
    y = "Storage relative to source TSV.gz",
    colour = NULL, shape = NULL,
    caption = paste0("Only formats that store a variant identity key are shown. ",
                     "Lower time and higher ratio are better.")
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey35"),
        plot.caption = element_text(colour = "grey35", hjust = 0))

ggsave(paste0(output_base, ".png"), p, width = 12, height = 8, dpi = 180, bg = "white")
fwrite(dat, "inst/benchmarks/pareto-chr1-frontier.csv")
message("Wrote ", output_base, ".png")
