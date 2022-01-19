#!/bin/bash
#SBATCH --job-name=rep1_mlttrk
#SBATCH --account=def-maxwl
#SBATCH --cpus-per-task=3
#SBATCH --mem=3G
#SBATCH --time=1-00:10
#SBATCH --output=%j.out
#SBATCH --array=0-20

# Your experiment setup logic here
source ~/miniconda3/etc/profile.d/conda.sh
conda activate segway_env

# Note the actual command is run through
PYTHONDONTWRITEBYTECODE= python batch_run02.py
