#!/bin/bash
#SBATCH --job-name=rep1_mlttrk
#SBATCH --account=def-maxwl
#SBATCH --cpus-per-task=5
#SBATCH --mem=20G
#SBATCH --time=0-05:10
#SBATCH --output=%j.out
#SBATCH --array=0-10

# Your experiment setup logic here
source ~/miniconda3/etc/profile.d/conda.sh
conda activate segway_env

# Note the actual command is run through
PYTHONDONTWRITEBYTECODE= python batch_run05.py
