#!/usr/bin/env bash
set -euo pipefail

source_sumstats=$1
if [ "$#" -ge 2 ]; then
  output_dir=$2
else
  output_dir=vcf-tabix-benchmark
fi

mkdir -p "$output_dir/sorttmp"

gzip -dc "$source_sumstats" |
  awk 'BEGIN {
    OFS = sprintf("%c", 9)
    print "##fileformat=VCFv4.3"
    for (i = 1; i <= 22; i++) print "##contig=<ID=" i ">"
    print "##contig=<ID=X>"
    print "##contig=<ID=Y>"
    print "##contig=<ID=MT>"
    print "##INFO=<ID=BETA,Number=1,Type=Float,Description=Effect estimate>"
    print "##INFO=<ID=SE,Number=1,Type=Float,Description=Standard error>"
    print "##INFO=<ID=P,Number=1,Type=Float,Description=P value>"
    print "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"
  }
  NR > 1 {
    print $2, $3, $1, $5, $4, ".", "PASS",
      "BETA=" $6 ";SE=" $7 ";P=" $8
  }' |
  bcftools sort -T "$output_dir/sorttmp" -m 4G -Oz +    -o "$output_dir/sumstats.vcf.gz"

tabix -f -p vcf "$output_dir/sumstats.vcf.gz"
echo "Wrote $output_dir/sumstats.vcf.gz and .tbi"
