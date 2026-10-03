#!/bin/bash -l
#SBATCH -A sscm013902
#SBATCH -p test
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 0:59:00
#SBATCH --array=0-6
#SBATCH -o /user/work/fh6520/cs-quant/logs/mr-%A_%a.log
module load languages/R/4.5.1
AR=(exact z9/eaf8/se6 z9/eaf8/se7 z9/eaf8/se8 z10/eaf8/se6 z10/eaf8/se7 z10/eaf8/se8)
Rscript /user/work/fh6520/cs-quant/scripts/mr.R ${AR[$SLURM_ARRAY_TASK_ID]}
