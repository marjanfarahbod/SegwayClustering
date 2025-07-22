#!/bin/bash
#SBATCH --job-name=rep1_mlttrk
#SBATCH --account=rrg-maxwl
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --time=1-0:0
#SBATCH --output=download.out

# setting up environment
module load python/3.8.10
source /home/mfarahbo/segwayENV/bin/activate

python SegwayDataDownloadFromENCODE.py
