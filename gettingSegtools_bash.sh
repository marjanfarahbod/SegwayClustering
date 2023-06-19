#!/bin/bash
#SBATCH --job-name=rep1_mlttrk
#SBATCH --account=rrg-maxwl
#SBATCH --cpus-per-task=2
#SBATCH --mem=30G
#SBATCH --time=0-23:00
#SBATCH --output=%j_segtoolsAllButGMTK_1.out
#SBATCH --array=1

source /home/mfarahbo/scratch/segtools-env2/bin/activate
module load python/2.7.18
module load r/4.2.1
module load hdf5/1.10.6

python wrapperForSegtools_allButGMTK.py
