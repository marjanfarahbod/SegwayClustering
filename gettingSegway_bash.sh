#!/bin/bash
#SBATCH --job-name=rep1_mlttrk
#SBATCH --account=rrg-maxwl
#SBATCH --cpus-per-task=2
#SBATCH --mem=30G
#SBATCH --time=0-23:00
#SBATCH --output=%j_segway_3_15.out
#SBATCH --array=3:15

module load python/3.8.10
module load hdf5
module load StdEnv/2020
module load gcc/9.3.0
module load gmtk
source /home/mfarahbo/scratch/segwayENV/bin/activate

python wrapperForSegway.py