#!/bin/bash -l
#SBATCH -A sscm013902
#SBATCH -p short
#SBATCH -c 8
#SBATCH --mem=12G
#SBATCH -t 0:25:00
#SBATCH --array=10-17
#SBATCH --constraint=Cascadelake
#SBATCH -o /user/work/fh6520/cs-quant/logs/grid-%A_%a.log
module load languages/R/4.5.1
DS=(finngen exp01 out01)
PR=(z9/eaf8/se6 z9/eaf8/se7 z9/eaf8/se8 z10/eaf8/se6 z10/eaf8/se7 z10/eaf8/se8)
i=$SLURM_ARRAY_TASK_ID
Rscript /user/work/fh6520/cs-quant/scripts/grid.R ${DS[$((i / 6))]} ${PR[$((i % 6))]}
