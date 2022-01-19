#!/bin/bash
#SBATCH --job-name=rep1_mlttrk
#SBATCH --account=def-maxwl
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0-00:10
#SBATCH --output=rep1_batch.out

# Your experiment setup logic here
source ~/miniconda3/etc/profile.d/conda.sh
conda activate segway_env

# Note the actual command is run through srun
PYTHONDONTWRITEBYTECODE= segway train --include-coords=data/encodePilotRegions.hg19.bed --num-labels=15 --resolution=100 --minibatch-fraction=0.01 --num-instances=10 --prior-strength=1.0 --segtransition-weight-scale=1.0 --ruler-scale=100 --track-weight=0.002 --max-train-rounds=25 data/DND-41.genomedata DND41train
segway posterior --include-coords=data/encodePilotRegions.hg19.bed data/DND-41.genomedata DND41train DND41posterior


