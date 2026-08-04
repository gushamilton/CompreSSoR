#!/usr/bin/env Rscript

# Build five scoped Pareto plots from the five BP optimization screens, then
# replace the historical Pcodec point with the locked five-run rerun.  The
# factor plots deliberately keep their access scopes in the title/caption;
# they are evidence plots, not interchangeable whole-file benchmarks.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

bench_root <- Sys.getenv("COMPRESSOR_FACTOR_BENCH_ROOT",
                        unset = "inst/benchmarks/final-factors-20260804")
fig_root <- Sys.getenv("COMPRESSOR_FACTOR_FIG_ROOT",
                       unset = "inst/figures")
dir.create(bench_root, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_root, recursive = TRUE, showWarnings = FALSE)

old_summary_path <- Sys.getenv(
  "COMPRESSOR_OLD_FINAL_SUMMARY",
  unset = "inst/benchmarks/final-keyed-20260804/final-keyed-summary.csv"
)
old_summary <- fread(old_summary_path)
source_bytes <- old_summary$source_bytes[1L]

selected_runs_path <- file.path(bench_root, "final-selected-runs.csv")
if (!file.exists(selected_runs_path)) {
  stop("missing selected rerun: ", selected_runs_path)
}
selected_runs <- fread(selected_runs_path)
if (any(selected_runs$status != "OK") || uniqueN(selected_runs$run_id) != 5L) {
  stop("selected rerun is not five successful independent runs")
}
selected_summary <- selected_runs[, .(
  format = format[1L], family = family[1L], identity_included = all(identity_included),
  rows = rows[1L], source_bytes = source_bytes[1L],
  storage_bytes = median(storage_bytes),
  compression_ratio = median(compression_ratio),
  write_median_seconds = median(write_seconds),
  full_read_median_seconds = median(full_read_seconds),
  max_abs_z_error = max(max_abs_z_error),
  max_abs_se_error = max(max_abs_se_error),
  max_abs_eaf_error = max(max_abs_eaf_error),
  checksum_consistent = uniqueN(checksum) == 1L
), by = format_id]
fwrite(selected_summary, file.path(bench_root, "final-selected-summary.csv"))

make_frontier <- function(x) {
  x[, pareto_frontier := vapply(seq_len(.N), function(i) {
    !any(storage_ratio >= storage_ratio[i] & access_seconds <= access_seconds[i] &
           (storage_ratio > storage_ratio[i] | access_seconds < access_seconds[i]))
  }, logical(1))]
  x[]
}

save_factor_plot <- function(x, title, subtitle, caption, filename,
                             xlab = "Access time (seconds; log scale)") {
  x <- copy(x)
  x <- x[is.finite(access_seconds) & is.finite(storage_ratio) & access_seconds > 0]
  x <- make_frontier(x)
  x[, label := ifelse(pareto_frontier, label, NA_character_)]
  frontier <- x[pareto_frontier == TRUE][order(access_seconds)]
  p <- ggplot(x, aes(access_seconds, storage_ratio, colour = label_group,
                     shape = label_group))
  if (nrow(frontier) > 1L) {
    p <- p + geom_path(data = frontier, aes(access_seconds, storage_ratio),
                       inherit.aes = FALSE, colour = "grey30", linewidth = 0.65)
  }
  p <- p + geom_point(size = 3, alpha = 0.9) +
    ggrepel::geom_text_repel(
      data = x[!is.na(label)], aes(label = label), seed = 20260804,
      box.padding = 0.5, point.padding = 0.25, min.segment.length = 0,
      max.overlaps = Inf, show.legend = FALSE, size = 3.15
    ) +
    scale_x_log10() +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.10))) +
    labs(title = title, subtitle = subtitle, x = xlab,
         y = "Storage relative to source TSV.gz", colour = NULL, shape = NULL,
         caption = caption) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          plot.caption = element_text(colour = "grey35", hjust = 0))
  ggsave(file.path(fig_root, filename), p, width = 11, height = 7, dpi = 180,
         bg = "white")
  fwrite(x, file.path(bench_root, sub("\\.png$", ".csv", filename)))
  invisible(x)
}

selected_storage <- selected_summary$storage_bytes[1L]
selected_ratio <- source_bytes / selected_storage

# 1. Public access policy: four threads for full keyed access, one for
# selective access.  This plot is the directly comparable whole-file read.
threads <- fread(file.path(bench_root, "public-thread-probe.csv"))
thread_dat <- threads[,
  .(access_seconds = median_seconds,
    storage_ratio = selected_ratio,
    label = paste0("Public read; ", threads, " thread",
                   ifelse(threads == 1, "", "s")),
    label_group = "Thread policy")]
save_factor_plot(thread_dat,
  "Thread policy Pareto screen",
  "FinnGen chr1 final Pcodec store; full keyed read; three repeats per setting",
  paste0("All points use the same bytes. The package default is four threads for a full store;\n",
         "region and sparse-key reads default to one thread."),
  "pareto-factor-threads.png")

# 2. Rust CPU flags: the benchmark's full_numeric workload isolates decoding
# without mixing in key materialization.  It is therefore labelled as a
# numeric-only build screen.
cpu_files <- c(
  portable = "opt2-cpu-portable.csv",
  native = "opt2-cpu-native.csv",
  avx2_bmi = "opt2-cpu-avx2.csv",
  skylake_avx512 = "opt2-cpu-avx512.csv"
)
cpu <- rbindlist(lapply(names(cpu_files), function(nm) {
  x <- fread(file.path(bench_root, cpu_files[[nm]]))
  x[, build := nm]
  x
}), fill = TRUE)
cpu_dat <- cpu[store == "finngen_chr1_final" & workload == "full_numeric",
  .(access_seconds = median_seconds,
    storage_ratio = selected_ratio,
    label = paste0("Pcodec ", flag, " build"),
    label_group = "Rust build")]
save_factor_plot(cpu_dat,
  "Rust build flags Pareto screen",
  "FinnGen chr1 final Pcodec store; numeric-only decode; three repeats per setting",
  paste0("CPU-specific target flags did not improve the measured store. The portable build is\n",
         "locked as the distribution default; COMPRESSOR_RUSTFLAGS remains an override."),
  "pareto-factor-cpu.png")

# 3. Destination-buffer experiment.  Sum the real stored streams within each
# repetition so every point represents a complete numeric-store decode.
buffer <- fread(file.path(bench_root, "opt2-buffer.csv"))
buffer_dat <- buffer[store == "finngen_chr1_final" & source_kind == "actual" &
                       status == "success",
  .(access_seconds = sum(seconds), storage_bytes = sum(compressed_bytes)),
  by = .(mode, run)]
buffer_dat <- buffer_dat[, .(access_seconds = median(access_seconds),
                             storage_ratio = source_bytes / median(storage_bytes),
                             label = fifelse(mode == "alloc_then_copy_current",
                                             "Current allocation + copy",
                                             fifelse(mode == "direct_destination_fresh",
                                                     "Direct destination; fresh",
                                                     "Direct destination; reused")),
                             label_group = "Buffer path"), by = mode]
save_factor_plot(buffer_dat,
  "Direct-buffer Pareto screen",
  "FinnGen chr1 final store; actual Pcodec streams; full numeric decode",
  paste0("This is a native buffer microbenchmark, not the public keyed R read. On the actual\n",
         "stored streams, direct destination buffers were effectively neutral, so the current\n",
         "allocation path remains locked."),
  "pareto-factor-buffer.png")

# 4. Fused native-reader screen.  The fused prototype produced a valid numeric
# checksum but its keyed checksum was NaN; do not present that unvalidated path
# as the production implementation.
fused_line <- paste(readLines(file.path(bench_root, "opt2-fused.txt")), collapse = " ")
fused_seconds <- as.numeric(sub(".*seconds=([0-9.]+).*", "\\1", fused_line))
current_seconds <- cpu[store == "finngen_chr1_final" & workload == "full_numeric" &
                         flag == "portable" & threads == 1, median_seconds]
fused_dat <- data.table(
  access_seconds = c(current_seconds, fused_seconds),
  storage_ratio = selected_ratio,
  label = c("Current native numeric path", "Fused numeric prototype"),
  label_group = "Native reader"
)
save_factor_plot(fused_dat,
  "Fused-reader Pareto screen",
  "Numeric-only prototype; keyed validation result was NaN and is excluded",
  paste0("The fused path is retained as an experiment only. It is not locked into the package\n",
         "because the complete keyed round-trip did not validate."),
  "pareto-factor-fused.png")

# 5. I/O geometry and connection policy.  Keep only successful full reads at
# the production stream-major/shared-connection/one-thread access point.
io <- fread(file.path(bench_root, "opt2-io-factorial.csv"))
io_dat <- io[status == "success" & workload == "full" & access == "stream_major" &
              connection == "shared_connection" & threads == 1,
  .(access_seconds = median(median_seconds), storage_bytes = median(bytes_read)),
  by = geometry]
io_dat[, `:=`(storage_ratio = source_bytes / storage_bytes,
             label = paste0("Pcodec ", geometry, "-row key blocks"),
             label_group = "I/O geometry")]
save_factor_plot(io_dat,
  "I/O geometry Pareto screen",
  "FinnGen chr1 store; full stream-major read; shared connection; one thread",
  paste0("The existing 131,072-row key geometry is retained: file size and access differences\n",
         "were sub-1%, while it is already validated by the production reader."),
  "pareto-factor-io.png")

# Summary: replace only the historical Pcodec point with the locked rerun and
# retain the complete 34-format keyed comparison for context.
summary <- old_summary[format_id != "pcodec_native"]
replacement <- old_summary[format_id == "pcodec_native"]
replacement[, `:=`(
  storage_bytes = selected_storage,
  compression_ratio = selected_ratio,
  write_median_seconds = selected_summary$write_median_seconds[1L],
  full_read_median_seconds = selected_summary$full_read_median_seconds[1L],
  max_abs_z_error = selected_summary$max_abs_z_error[1L],
  max_abs_se_error = selected_summary$max_abs_se_error[1L],
  max_abs_eaf_error = selected_summary$max_abs_eaf_error[1L],
  checksum_consistent = selected_summary$checksum_consistent[1L]
)]
summary <- rbindlist(list(summary, replacement), fill = TRUE)
summary[, pareto_frontier := vapply(seq_len(.N), function(i) {
  !any(compression_ratio >= compression_ratio[i] &
         full_read_median_seconds <= full_read_median_seconds[i] &
         (compression_ratio > compression_ratio[i] |
            full_read_median_seconds < full_read_median_seconds[i]))
}, logical(1))]
setorder(summary, full_read_median_seconds, -compression_ratio)
fwrite(summary, file.path(bench_root, "final-summary-updated.csv"))
fwrite(summary[pareto_frontier == TRUE], file.path(bench_root, "final-frontier-updated.csv"))

summary_plot <- copy(summary)
summary_plot[, label := ifelse(pareto_frontier | format_id == "tsv_gzip", format, NA_character_)]
summary_plot[format_id == "pcodec_native", label := "CompreSSoR native Pcodec\n(locked default)"]
summary_plot[, label_group := family]
frontier <- summary_plot[pareto_frontier == TRUE][order(full_read_median_seconds)]
summary_plot[, family := factor(family, levels = c("CompreSSoR", "Binary",
                                                    "Framed binary", "Parquet",
                                                    "VCF + Tabix", "Text"))]
frontier[, family := factor(family, levels = levels(summary_plot$family))]
p <- ggplot(summary_plot, aes(full_read_median_seconds, compression_ratio,
                              colour = family, shape = family)) +
  geom_path(data = frontier, aes(full_read_median_seconds, compression_ratio),
            inherit.aes = FALSE, colour = "grey30", linewidth = 0.7) +
  geom_point(size = 2.7, alpha = 0.9) +
  ggrepel::geom_text_repel(data = summary_plot[!is.na(label)], aes(label = label),
                           seed = 20260804, box.padding = 0.55,
                           point.padding = 0.3, min.segment.length = 0,
                           max.overlaps = Inf, show.legend = FALSE, size = 3.2) +
  scale_x_log10() +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.09))) +
  labs(title = "Locked CompreSSoR Pcodec against keyed alternatives",
       subtitle = paste0(nrow(summary), " keyed formats; historical five-run BP screen;\n",
                         "Pcodec point replaced by the locked five-run rerun"),
       x = "Whole-file access time (seconds; log scale)",
       y = "Storage compression relative to TSV.gz",
       colour = "Format family", shape = "Format family",
       caption = paste0("Labels show only Pareto-frontier formats plus the TSV.gz baseline.\n",
                        "Every plotted file includes the self-contained variant identity key.")) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.caption = element_text(colour = "grey35", hjust = 0))
ggsave(file.path(fig_root, "compressor-pareto-summary.png"), p,
       width = 13, height = 8, dpi = 180, bg = "white")
ggsave(file.path(fig_root, "compressor-pareto-summary.pdf"), p,
       width = 13, height = 8, device = cairo_pdf, bg = "white")

message("Wrote five factor plots and summary plot to ", fig_root)
