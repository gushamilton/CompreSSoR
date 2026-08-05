# Archived with the pre-0.5 harmonisation workflow. This file is not installed.

#' Lift GWAS summary statistics to GRCh38
#'
#' @param input A data frame or summary-statistics file accepted by
#'   `import_sumstats()`.
#' @param input_build Source genome build, normally `"GRCh37"`.
#' @param chain A local GRCh37-to-GRCh38 chain file.
#' @param drop_unmapped Drop rows that do not map exactly once.
#' @return A canonical data.frame with liftover diagnostics.
liftover_sumstats <- function(input, input_build = "GRCh37", chain = NULL,
                              drop_unmapped = TRUE) {
  if (length(drop_unmapped) != 1L || !is.logical(drop_unmapped) || is.na(drop_unmapped)) {
    stop("drop_unmapped must be TRUE or FALSE", call. = FALSE)
  }
  data <- import_sumstats(input)
  lifted <- lift_table_to_grch38(data, input_build = input_build, chain = chain)
  status <- lifted$.compressor_liftover_status %||% rep("not_needed", nrow(lifted))
  stats <- as.list(table(status, useNA = "ifany"))
  if (isTRUE(drop_unmapped)) {
    keep <- status %in% c("mapped", "not_needed")
    lifted <- lifted[keep, , drop = FALSE]
    status <- status[keep]
  }
  lifted$liftover_status <- status
  if (".compressor_source_variant_id" %in% names(lifted)) {
    lifted$source_variant_id <- lifted$.compressor_source_variant_id
    lifted$.compressor_source_variant_id <- NULL
  }
  lifted$.compressor_liftover_status <- NULL
  lifted$.compressor_liftover_row <- NULL
  row.names(lifted) <- NULL
  attr(lifted, "liftover_stats") <- stats
  attr(lifted, "genome_build") <- "GRCh38"
  lifted
}
