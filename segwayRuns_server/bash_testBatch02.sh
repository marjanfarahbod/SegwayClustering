#!/bin/bash
#SBATCH --job-name=rep1_mlttrk
#SBATCH --account=def-maxwl
#SBATCH --cpus-per-task=3
#SBATCH --mem=3G
#SBATCH --time=0-01:10
#SBATCH --output=rep1_batch.out

# Your experiment setup logic here
source ~/miniconda3/etc/profile.d/conda.sh
conda activate segway_env

# Just run python
python batch_run01.py
