#!/bin/bash -l
#SBATCH -A sscm013902
#SBATCH -p test
#SBATCH -c 2
#SBATCH --mem=12G
#SBATCH -t 0:30:00
#SBATCH -o /user/work/fh6520/cs-quant/logs/compat-%j.log
module load languages/R/4.5.1
Q=/user/work/fh6520/cs-quant
S="$Q/work/legacy-z9se6-a27b32d.cpr /user/work/fh6520/showcase/prep/traits/cpr/exp01.cpr /user/work/fh6520/showcase/prep/traits/cpr/out01.cpr"
Rscript $Q/scripts/compat.R $Q/lib-base $Q/results/compat_old.rds $S
Rscript $Q/scripts/compat.R $Q/lib3 $Q/results/compat_new.rds $S
Rscript -e "o <- readRDS(\"$Q/results/compat_old.rds\"); n <- readRDS(\"$Q/results/compat_new.rds\"); cat(o\$version, \"vs\", n\$version, \"\n\"); for (s in names(o)[-(1:2)]) cat(s, \"identical reads:\", sum(o[[s]] == n[[s]]), \"of\", length(o[[s]]), \"\n\"); print(o[[3]])"
