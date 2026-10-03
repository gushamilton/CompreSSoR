Q <- "/user/work/fh6520/cs-quant"; SHOW <- "/user/work/fh6520/showcase"
.libPaths(c(file.path(Q, "lib"), "/user/work/fh6520/r_packages"))
suppressPackageStartupMessages({library(data.table); library(CompreSSoR)})
e <- fread(file.path(Q, "results/mr/exact_mr.csv")); 
for (p in c("z9_eaf8_se8", "z10_eaf8_se8", "z10_eaf8_se6")) {
  a <- fread(file.path(Q, "results/mr", paste0(p, "_mr.csv")))
  j <- merge(a, e, by = c("id.exposure", "id.outcome", "method"))[method == "Inverse variance weighted"]
  j[, d := (b.x - b.y) / se.y]; assign(paste0("d_", p), j[order(id.exposure, id.outcome)]$d)
}
cat("cor z9se8 vs z10se8 signed IVW dev:", cor(d_z9_eaf8_se8, d_z10_eaf8_se8), "\n")
print(summary(abs(d_z10_eaf8_se8))); print(summary(abs(d_z9_eaf8_se8)))
inst <- readRDS(file.path(SHOW, "e2e/results/rep1/10x10_C8/result.rds"))$instruments
# per-SNP outcome errors for exp01 -> out01..03, z10se8 store vs TSV
for (o in c("out01", "out02")) {
  k <- inst$exp01
  t <- fread(file.path(SHOW, "prep/traits/tsv", paste0(o, ".tsv.gz")))
  t[, key := compressor_variant_key(chromosome, base_pair_location, other_allele, effect_allele)]
  t <- t[match(k, key)]
  s <- as.data.table(read_sumstats(file.path(Q, "work/mr/z10_eaf8_se8", paste0(o, ".cpr")), variants = k,
        columns = c("chromosome","base_pair_location","effect_allele","other_allele","beta","standard_error","z")))
  s[, key := compressor_variant_key(chromosome, base_pair_location, other_allele, effect_allele)]
  s <- s[match(k, key)]
  zt <- t$beta / t$standard_error
  cat(o, "max |dz|", max(abs(s$z - zt)), " max rel dse", max(abs(s$standard_error / t$standard_error - 1)),
      " max |dbeta|/se", max(abs(s$beta - t$beta) / t$standard_error), " n |z|>=3.5:", sum(abs(zt) >= 3.5), "\n")
}
