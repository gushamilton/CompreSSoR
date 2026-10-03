#!/bin/bash -l
#SBATCH -A sscm013902
#SBATCH -p test
#SBATCH -c 16
#SBATCH --mem=64G
#SBATCH -t 0:59:00
#SBATCH --array=0-5
#SBATCH -o /user/work/fh6520/cs-quant/logs/mr2-%A_%a.log
module load languages/R/4.5.1
export CSQ_LIB=lib2
AR=(z9/eaf8/se6+xse z9/eaf8/se7+xse z9/eaf8/se8+xse z10/eaf8/se6+xse z10/eaf8/se7+xse z10/eaf8/se8+xse)
Rscript /user/work/fh6520/cs-quant/scripts/mr.R ${AR[$SLURM_ARRAY_TASK_ID]}
