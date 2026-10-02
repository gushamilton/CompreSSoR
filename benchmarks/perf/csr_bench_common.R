# Shared helpers for CompreSSoR read/ingest performance benchmarks.
# Self-contained: depends only on CompreSSoR + data.table.

suppressPackageStartupMessages({
  library(CompreSSoR)
  library(data.table)
})

csr_grch38_lengths <- c(
  248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
  145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
  101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983,
  50818468)

# Realistic-ish synthetic GWAS: sorted unique positions across chr1-22,
# random distinct REF/ALT, Beta-distributed EAF, SE from N and EAF, null Z
# plus `n_signals` causal loci with LD-like decaying Z peaks (so p<=5e-8 hits
# cluster in regions, like real data). Deterministic given `seed`.
csr_synthetic_gwas <- function(n, seed = 1L, n_signals = 60L, N = 50000, trait = 1L) {
  set.seed(seed)
  w <- csr_grch38_lengths / sum(csr_grch38_lengths)
  per_chr <- as.integer(round(n * w))
  per_chr[1] <- per_chr[1] + (n - sum(per_chr))
  chr <- rep(as.character(1:22), per_chr)
  pos <- unlist(lapply(seq_along(per_chr), function(i) {
    sort(sample.int(csr_grch38_lengths[i] - 1L, per_chr[i]))
  }), use.names = FALSE)
  bases <- c("A", "C", "G", "T")
  ref_i <- sample.int(4L, n, replace = TRUE)
  alt_i <- ((ref_i - 1L + sample.int(3L, n, replace = TRUE)) %% 4L) + 1L
  # Shared variant panel across traits (same seed); trait-specific statistics.
  set.seed(seed * 7919L + trait)
  eaf <- pmin(0.995, pmax(0.005, rbeta(n, 0.6, 0.6)))
  se <- 1 / sqrt(2 * N * eaf * (1 - eaf)) * exp(rnorm(n, 0, 0.03))
  z <- rnorm(n)
  # signals: peak Z ~ 6..25, decaying over ~200 neighbouring rows
  centres <- sample.int(n, n_signals)
  for (c0 in centres) {
    peak <- runif(1, 6, 25) * sample(c(-1, 1), 1)
    idx <- max(1L, c0 - 400L):min(n, c0 + 400L)
    z[idx] <- z[idx] + peak * exp(-abs(idx - c0) / 120) * runif(length(idx), 0.6, 1)
  }
  data.table(
    chromosome = chr, base_pair_location = as.integer(pos),
    reference_allele = bases[ref_i], alternate_allele = bases[alt_i],
    effect_allele = bases[alt_i], other_allele = bases[ref_i],
    beta = z * se, standard_error = se, effect_allele_frequency = eaf,
    p_value = 2 * pnorm(-abs(z)), z = z
  )
}

csr_dir_bytes <- function(path) {
  sum(file.info(list.files(path, full.names = TRUE, recursive = TRUE))$size,
      na.rm = TRUE)
}

# Order-independent numeric checksum of a data.frame (or list of them).
csr_checksum <- function(x) {
  if (is.null(x)) return(NA_real_)
  if (is.list(x) && !is.data.frame(x)) {
    return(sum(vapply(x, csr_checksum, numeric(1))))
  }
  if (is.data.frame(x)) {
    return(sum(vapply(x, function(v) {
      if (is.character(v)) sum(nchar(v)) + length(unique(v))
      else sum(as.numeric(v), na.rm = TRUE)
    }, numeric(1))) + nrow(x))
  }
  sum(as.numeric(x), na.rm = TRUE) + length(x)
}

csr_env_row <- function() {
  cpu <- tryCatch({
    if (file.exists("/proc/cpuinfo")) {
      l <- grep("^model name", readLines("/proc/cpuinfo"), value = TRUE)[1]
      trimws(sub("^model name\\s*:", "", l))
    } else system("sysctl -n machdep.cpu.brand_string", intern = TRUE)
  }, error = function(e) NA_character_)
  sha <- Sys.getenv("CSR_GIT_SHA", "")
  if (!nzchar(sha)) {
    sha <- tryCatch(system2("git", c("-C", Sys.getenv("CSR_REPO", "."),
                                     "rev-parse", "HEAD"), stdout = TRUE,
                            stderr = FALSE)[1], error = function(e) NA_character_)
  }
  data.table(
    hostname = Sys.info()[["nodename"]], cpu_model = cpu,
    nproc = parallel::detectCores(), r_version = as.character(getRversion()),
    compressor_version = as.character(packageVersion("CompreSSoR")),
    git_sha = sha, label = Sys.getenv("CSR_LABEL", "unlabelled"),
    slurm_job = Sys.getenv("SLURM_JOB_ID", ""),
    slurm_task = Sys.getenv("SLURM_ARRAY_TASK_ID", "")
  )
}
