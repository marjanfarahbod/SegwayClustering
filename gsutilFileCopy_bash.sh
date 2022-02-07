#!/bin/bash
#SBATCH --job-name=rep1_mlttrk
#SBATCH --account=rrg-maxwl
#SBATCH --cpus-per-task=5
#SBATCH --mem=5G
#SBATCH --time=1-00:00
#SBATCH --output=rep1_batch.out

# setting up environment
module load python/3.8
source /home/mfarahbo/gsutilENV/bin/activate

python gsutilFileCopy.py
