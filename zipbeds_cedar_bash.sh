#!/bin/bash
#SBATCH --job-name=rep1_mlttrk
#SBATCH --account=rrg-maxwl
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-2:00
#SBATCH --output=%j_zipbeds_38.out

module load python/3.8.10

python zipbeds_cedar.py
