#!/usr/bin/env Rscript

# Plot the two downstream MR access paths from one or more full-file runs.
# Inputs and outputs are compact; raw stores remain outside the repository.

input_paths <- strsplit(
  Sys.getenv("COMPRESSOR_PVALUE_PLOT_INPUTS"), ",", fixed = TRUE
)[[1L]]
input_paths <- trimws(input_paths[nzchar(trimws(input_paths))])
if (!length(input_paths)) stop("COMPRESSOR_PVALUE_PLOT_INPUTS is required")
labels <- strsplit(
  Sys.getenv("COMPRESSOR_PVALUE_PLOT_LABELS", unset = ""), ",", fixed = TRUE
)[[1L]]
labels <- trimws(labels[nzchar(trimws(labels))])
if (length(labels) != length(input_paths)) {
  labels <- basename(dirname(input_paths))
}
result_root <- Sys.getenv(
  "COMPRESSOR_PVALUE_PLOT_RESULTS",
  unset = file.path(getwd(), "inst/benchmarks/pvalue-sidechannel-macmini-20260806")
)
dir.create(result_root, recursive = TRUE, showWarnings = FALSE)

tables <- lapply(seq_along(input_paths), function(i) {
  if (!file.exists(input_paths[[i]])) stop("summary CSV not found: ", input_paths[[i]])
  value <- read.csv(input_paths[[i]], check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("threshold", "source_hit_rows", "p_filter_fetch_seconds",
                "flag_fetch_seconds")
  missing <- setdiff(required, names(value))
  if (length(missing)) stop("summary is missing: ", paste(missing, collapse = ", "))
  value$dataset <- labels[[i]]
  value
})
combined <- do.call(rbind, tables)
combined <- combined[order(combined$dataset, -log10(combined$threshold)), , drop = FALSE]
combined$relative_savings_percent <- 100 *
  (1 - combined$flag_fetch_seconds / combined$p_filter_fetch_seconds)
write.csv(combined,
          file.path(result_root, "pvalue-sidechannel-mr-comparison.csv"),
          row.names = FALSE)

series <- c("p_filter_fetch_seconds", "flag_fetch_seconds")
series_labels <- c("p filtering + MR field fetch",
                   "flag inspection + MR field fetch")
series_colours <- c("#2C7FB8", "#D95F0E")

finite_values <- c(combined[[series[[1L]]]], combined[[series[[2L]]]])
finite_values <- finite_values[is.finite(finite_values)]
y_limit <- range(finite_values, c(0, 0))
y_limit[[2L]] <- max(y_limit[[2L]] * 1.15, 0.01)

draw_plot <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  par(mfrow = c(1L, length(input_paths)), mar = c(6.5, 5, 3, 1),
      oma = c(0, 0, 2, 0))
  for (i in seq_along(input_paths)) {
    value <- combined[combined$dataset == labels[[i]], , drop = FALSE]
    value <- value[order(-log10(value$threshold)), , drop = FALSE]
    x <- seq_len(nrow(value))
    axis_labels <- paste0(
      "p<=1e-", format(round(-log10(value$threshold)), trim = TRUE),
      "\n(", format(round(value$source_hit_rows), big.mark = ","), " hits)"
    )
    plot(x, value[[series[[1L]]]], type = "o", pch = 16, lwd = 2,
         col = series_colours[[1L]], ylim = y_limit, xaxt = "n",
         xlab = "selection threshold", ylab = "seconds",
         main = labels[[i]])
    axis(1L, at = x, labels = axis_labels, las = 1L, cex.axis = 0.8)
    lines(x, value[[series[[2L]]]], type = "o", pch = 17, lwd = 2,
          col = series_colours[[2L]])
    legend("topleft", legend = series_labels, col = series_colours,
           pch = c(16, 17), lwd = 2, bty = "n", cex = 0.8)
  }
  mtext("MR-relevant selected-field access: p filtering versus aligned flag",
        outer = TRUE, cex = 1.1, font = 2)
}

draw_relative_plot <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  par(mar = c(6.5, 5, 3, 1))
  datasets <- labels
  value_by_dataset <- lapply(datasets, function(dataset) {
    value <- combined[combined$dataset == dataset, , drop = FALSE]
    value[order(-log10(value$threshold)), , drop = FALSE]
  })
  x <- seq_len(nrow(value_by_dataset[[1L]]))
  y_values <- unlist(lapply(value_by_dataset, function(value) value$relative_savings_percent))
  y_limit <- c(min(0, y_values), max(10, max(y_values) * 1.15))
  plot(x, value_by_dataset[[1L]]$relative_savings_percent, type = "o",
       pch = 16, lwd = 2, col = series_colours[[1L]], ylim = y_limit,
       xaxt = "n", xlab = "selection threshold", ylab = "time saved (%)",
       main = "Relative benefit of flag inspection")
  axis(1L, at = x,
       labels = paste0("p<=1e-", format(round(-log10(value_by_dataset[[1L]]$threshold)), trim = TRUE)),
       las = 1L)
  for (i in seq_along(value_by_dataset)) {
    value <- value_by_dataset[[i]]
    lines(x, value$relative_savings_percent, type = "o", pch = 15 + i,
          lwd = 2, col = series_colours[[i]])
  }
  abline(h = 0, col = "grey60", lty = 2)
  legend("topleft", legend = datasets, col = series_colours,
         pch = c(16, 17), lwd = 2, bty = "n")
}

png(file.path(result_root, "pvalue-sidechannel-mr-comparison.png"),
    width = max(1200, 900 * length(input_paths)), height = 900, res = 140)
draw_plot()
dev.off()

png(file.path(result_root, "pvalue-sidechannel-mr-relative.png"),
    width = 1200, height = 850, res = 140)
draw_relative_plot()
dev.off()

cat("Wrote pvalue-sidechannel-mr-comparison.csv and absolute/relative .png in ",
    result_root, "\n", sep = "")
