#!/usr/bin/env Rscript

# Bounded regression runner for the archived reference-backed harmoniser.
# The live package deliberately does not load archive/harmonisation/.

suppressPackageStartupMessages({
  library(CompreSSoR)
  library(testthat)
})

source_root <- normalizePath(".", mustWork = TRUE)
source_env <- new.env(parent = asNamespace("CompreSSoR"))
archive_root <- file.path(source_root, "archive", "harmonisation")
source_files <- list.files(file.path(archive_root, "R"), pattern = "[.]R$",
                           full.names = TRUE)
for (source_file in sort(source_files)) sys.source(source_file, source_env)

test_file <- file.path(archive_root, "tests", "test-observability-legacy.R")
testthat::test_file(test_file, reporter = "summary", env = source_env)
