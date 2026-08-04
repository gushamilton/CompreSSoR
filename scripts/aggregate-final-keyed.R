#!/usr/bin/env Rscript

# Aggregate the five independent BP final keyed-format runs and draw the
# headline Pareto plot.  All formats in the measured input must include the
# exact position + directed REF->ALT identity key.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

input_root <- Sys.getenv("COMPRESSOR_FINAL_ROOT", unset = "")
output_root <- Sys.getenv("COMPRESSOR_FINAL_OUTPUT", unset = "")
if (!nzchar(input_root) || !dir.exists(input_root)) stop("missing COMPRESSOR_FINAL_ROOT")
if (!nzchar(output_root)) stop("missing COMPRESSOR_FINAL_OUTPUT")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

files <- sort(list.files(input_root, pattern = "^final-keyed-run\\.csv$",
                         recursive = TRUE, full.names = TRUE))
if (length(files) != 5L) stop("expected exactly five run CSVs, found ", length(files))
runs <- rbindlist(lapply(files, fread), fill = TRUE)
if (uniqueN(runs$run_id) != 5L) stop("run IDs are not five independent runs")
if (any(runs$identity_included != TRUE)) stop("non-keyed candidate present")
if (any(runs$status != "OK")) stop("failed candidate present: ", paste(unique(runs[status != "OK", format]), collapse = ", "))

counts <- runs[, .(n = .N, runs = uniqueN(run_id)), by = format_id]
if (any(counts$n != 5L | counts$runs != 5L)) stop("not every format has five runs")

summary <- runs[, .(
  format = format[1L], family = family[1L], identity_included = all(identity_included),
  rows = rows[1L], source_bytes = source_bytes[1L],
  storage_bytes = as.numeric(median(storage_bytes)),
  compression_ratio = as.numeric(median(compression_ratio)),
  write_median_seconds = as.numeric(median(write_seconds)),
  write_min_seconds = as.numeric(min(write_seconds)),
  write_max_seconds = as.numeric(max(write_seconds)),
  full_read_median_seconds = as.numeric(median(full_read_seconds)),
  full_read_min_seconds = as.numeric(min(full_read_seconds)),
  full_read_max_seconds = as.numeric(max(full_read_seconds)),
  max_abs_z_error = max(max_abs_z_error),
  max_abs_se_error = max(max_abs_se_error),
  max_abs_eaf_error = max(max_abs_eaf_error),
  checksum_consistent = uniqueN(checksum) == 1L
), by = format_id][order(full_read_median_seconds, -compression_ratio)]

summary[, pareto_frontier := vapply(seq_len(.N), function(i) {
  !any(compression_ratio >= compression_ratio[i] &
       full_read_median_seconds <= full_read_median_seconds[i] &
       (compression_ratio > compression_ratio[i] |
        full_read_median_seconds < full_read_median_seconds[i]))
}, logical(1))]

registry <- unique(runs[, .(format_id, format, family, identity_included)])
registry[, `:=`(status = "MEASURED", reason = "Five-run keyed benchmark passed")]
unavailable <- data.table(
  format_id = c("zarr_pcodec", "zarr_blosc2", "blosc2_frame", "vortex", "python_pcodec"),
  format = c("Zarr + Pcodec", "Zarr + Blosc2", "Blosc2 frame", "Vortex", "Python Pcodec"),
  family = c("Zarr", "Zarr", "Blosc2", "Vortex", "Python"),
  identity_included = TRUE,
  status = "UNAVAILABLE",
  reason = c(
    "Python zarr/pcodec bindings were not installed in the BP R environment",
    "Python zarr/numcodecs/blosc2 bindings were not installed in the BP R environment",
    "blosc2 Python binding was not installed in the BP R environment",
    "Vortex runtime was not installed in the BP environment",
    "Python pcodec binding was not installed in the BP environment"
  )
)

fwrite(runs, file.path(output_root, "final-keyed-runs.csv"))
fwrite(summary, file.path(output_root, "final-keyed-summary.csv"))
fwrite(summary[pareto_frontier == TRUE], file.path(output_root, "final-keyed-frontier.csv"))
fwrite(rbindlist(list(registry, unavailable), fill = TRUE),
       file.path(output_root, "final-keyed-registry.csv"))

plot_data <- copy(summary)
plot_data[, family := factor(family, levels = c("CompreSSoR", "Binary", "Framed binary", "Parquet", "VCF + Tabix", "Text"))]
plot_data[, label := ifelse(pareto_frontier | format_id == "tsv_gzip", format, NA_character_)]
plot_data[format_id == "pcodec_native", label := "CompreSSoR native Pcodec\n(default)"]
plot_data[, point_size := ifelse(pareto_frontier, 3.1, 2.0)]

plot <- ggplot(plot_data, aes(full_read_median_seconds, compression_ratio,
                              colour = family, shape = family)) +
  geom_point(aes(size = point_size), alpha = 0.88) +
  scale_size_identity() +
  ggrepel::geom_text_repel(
    data = plot_data[!is.na(label)], aes(label = label),
    seed = 42, box.padding = 0.55, point.padding = 0.35,
    min.segment.length = 0, max.overlaps = Inf, direction = "both",
    show.legend = FALSE, size = 3.3
  ) +
  scale_x_log10() +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.08))) +
  labs(
    title = "Final keyed-format screen on FinnGen chr1",
    subtitle = paste0(nrow(summary), " measured formats; five independent BP runs; ",
                      "exact identity key included in every file"),
    x = "Whole-file access time (seconds; log scale)",
    y = "Storage compression relative to TSV.gz (higher is smaller)",
    colour = "Format family", shape = "Format family",
    caption = "Labels show only Pareto-frontier formats plus the TSV.gz baseline.\nAll measured candidates passed exact identity and numeric round-trip checks."
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(output_root, "compressor-pareto-final-keyed.png"), plot,
       width = 13, height = 8, dpi = 180, bg = "white")
ggsave(file.path(output_root, "compressor-pareto-final-keyed.pdf"), plot,
       width = 13, height = 8, device = cairo_pdf, bg = "white")
print(summary[, .(format, storage_bytes, compression_ratio,
                  full_read_median_seconds, pareto_frontier)])
