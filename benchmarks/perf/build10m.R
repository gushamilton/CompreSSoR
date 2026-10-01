source("csr_bench_common.R")
d <- csr_synthetic_gwas(1e7, seed = 1)
cat("hits 5e-8:", sum(d$p_value <= 5e-8), " 0.01:", sum(d$p_value<=0.01), "\n")
t0 <- proc.time()[[3]]
s <- compress_sumstats(as.data.frame(d), "stores/s10m.cpr", qc = "none", overwrite = TRUE, pvalue_order = TRUE, threads = 4)
cat("compress 10M qc=none t=4:", proc.time()[[3]]-t0, "\n")
str(s$manifest$preparation$phase_timings %||% s$manifest$timings)
print(names(s$manifest))
cat(csr_dir_bytes("stores/s10m.cpr"), "\n")
