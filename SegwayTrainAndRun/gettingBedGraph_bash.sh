#!/bin/bash
#SBATCH --job-name=rep1_mlttrk
#SBATCH --account=rrg-maxwl
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=0-03:00
#SBATCH --output=%j.out
#SBATCH --array=0-111

source /home/mfarahbo/segwayENV/activate

python wrapperForGettingBedGraph.py
