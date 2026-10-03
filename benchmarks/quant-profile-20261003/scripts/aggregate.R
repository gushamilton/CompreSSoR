#!/usr/bin/env Rscript
Q <- "/user/work/fh6520/cs-quant"
.libPaths(c(file.path(Q, "lib"), "/user/work/fh6520/r_packages",
            "/software/local/languages/miniforge3/envs/r-4.5.1/lib/R/library"))
suppressPackageStartupMessages(library(data.table))
R <- file.path(Q, "results"); g <- file.path(R, "grid"); m <- file.path(R, "mr")
rd <- function(dir, pat) rbindlist(lapply(list.files(dir, pat, full.names = TRUE), fread), fill = TRUE)
meta <- rd(g, "_meta.csv$"); errs <- rd(g, "_errors.csv$"); flips <- rd(g, "_flips.csv$")
base <- meta[profile == "z9/eaf8/se6", .(dataset, base_bytes = bytes)]
meta <- merge(meta, base, by = "dataset")[, size_vs_z9se6_pct := 100 * (bytes / base_bytes - 1)]
fwrite(meta[order(dataset, profile)], file.path(R, "grid_size_time.csv"))
fwrite(errs[order(dataset, field, profile)], file.path(R, "grid_errors.csv"))
fwrite(flips[order(dataset, reference, threshold, profile)], file.path(R, "grid_flips.csv"))
print(meta[order(dataset, profile), .(dataset, profile, mb = round(mb, 2), bpv = round(bytes_per_variant, 3),
  size_pct = round(size_vs_z9se6_pct, 1), z_mb = round(z_bytes / 1e6, 2), se_mb = round(se_bytes / 1e6, 2),
  write_s = round(write_s_med, 1), read_s = round(read_s_med, 2))])
print(dcast(errs[field %in% c("se", "z", "beta_in_se_units", "nlp_vs_exact_z")], dataset + profile ~ field,
            value.var = "max_rel")[])
print(errs[, .(dataset, profile, field, max_abs = signif(max_abs, 3), p999_abs = signif(p999_abs, 3),
               max_rel = signif(max_rel, 3), p999_rel = signif(p999_rel, 3))])
print(flips[reference == "exact_z_p", .(dataset, profile, threshold, n_ref, lost, gained, flips, flip_pct = round(flip_pct, 3))])
# MR
mr <- rd(m, "_mr.csv$"); key <- function(d) d[, k := paste(id.exposure, id.outcome, method, sep = "|")]
ex <- key(mr[arm == "exact"]); ar <- key(mr[arm != "exact"])
j <- merge(ar, ex[, .(k, b.ex = b, se.ex = se, nsnp.ex = nsnp)], by = "k")
j[, `:=`(db_se = abs(b - b.ex) / se.ex, dse_rel = abs(se - se.ex) / se.ex)]
est <- j[, .(pairs = .N, same_nsnp = mean(nsnp == nsnp.ex), med_db_se = median(db_se), p90_db_se = quantile(db_se, 0.9),
             max_db_se = max(db_se), med_dse_rel = median(dse_rel), max_dse_rel = max(dse_rel)), by = .(arm, method)]
fwrite(est[order(method, arm)], file.path(R, "mr_deviation.csv"))
print(est[order(method, arm)][, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 3) else x)])
conv <- rd(m, "_conversion.csv$")
if (nrow(conv)) { conv[, arm := NA_character_] }
cf <- list.files(m, "_conversion.csv$", full.names = TRUE)
conv <- rbindlist(lapply(cf, function(f) fread(f)[, arm := sub("_conversion.csv$", "", basename(f))]))
cs <- conv[, .(traits = .N, total_mb = sum(bytes) / 1e6, mean_encode_s = mean(encode_s)), by = arm]
fwrite(cs, file.path(R, "mr_store_sizes.csv")); print(cs)
