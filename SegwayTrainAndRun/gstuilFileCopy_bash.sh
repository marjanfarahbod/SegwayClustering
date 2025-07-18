#!/bin/bash
#SBATCH --job-name=rep1_mlttrk
#SBATCH --account=def-maxwl
#SBATCH --cpus-per-task=5
#SBATCH --mem=2G
#SBATCH --time=00-20:10
#SBATCH --output=rep1_batch.out

# setting up environment
module load python/3.8
source /home/mfarahbo/gsutilENV/bin/activate

python gsutilFileCopy.py
