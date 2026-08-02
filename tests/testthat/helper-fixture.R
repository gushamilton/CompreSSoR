make_fixture <- function(n = 1200L) {
  data.frame(
    chromosome = rep("1", n),
    base_pair_location = seq.int(100001L, length.out = n),
    effect_allele = rep(c("A", "C", "G", "T"), length.out = n),
    other_allele = rep(c("C", "G", "T", "A"), length.out = n),
    beta = sin(seq_len(n) / 13) / 5,
    standard_error = 0.02 + (seq_len(n) %% 17) / 1000,
    effect_allele_frequency = 0.1 + (seq_len(n) %% 20) / 100,
    p_value = 2 * pnorm(-abs(sin(seq_len(n) / 13) / 5 / (0.02 + (seq_len(n) %% 17) / 1000))),
    variant_id = paste0("1_", seq.int(100001L, length.out = n)),
    rsid = paste0("rs", seq_len(n)),
    annotation = rep(c("a", "b", "c"), length.out = n),
    stringsAsFactors = FALSE
  )
}
