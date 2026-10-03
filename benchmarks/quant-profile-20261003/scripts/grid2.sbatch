#!/bin/bash -l
#SBATCH -A sscm013902
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=12G
#SBATCH -t 0:25:00
#SBATCH --array=0-5
#SBATCH --constraint=Cascadelake
#SBATCH -o /user/work/fh6520/cs-quant/logs/grid2-%A_%a.log
module load languages/R/4.5.1
export CSQ_LIB=lib2
DS=(finngen exp01 out01)
PR=(z9/eaf8/se8+xse z10/eaf8/se8+xse)
i=$SLURM_ARRAY_TASK_ID
Rscript /user/work/fh6520/cs-quant/scripts/grid.R ${DS[$((i / 2))]} ${PR[$((i % 2))]}
