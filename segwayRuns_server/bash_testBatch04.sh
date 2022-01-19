#!/bin/bash
#SBATCH --job-name=rep1_mlttrk
#SBATCH --account=def-maxwl
#SBATCH --cpus-per-task=5
#SBATCH --mem=5G
#SBATCH --time=0-01:10
#SBATCH --output=bashtest04.out

# Your experiment setup logic here
source ~/miniconda3/etc/profile.d/conda.sh
conda activate segway_env

# Just run python
PYTHONDONTWRITEBYTECODE= python batch_run04.py
