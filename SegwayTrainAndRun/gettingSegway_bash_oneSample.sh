#!/bin/bash
#SBATCH --job-name=rep1_mlttrk
#SBATCH --account=rrg-maxwl
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --time=0-24:00
#SBATCH --output=%j.out


module load python/3.8.10
module load hdf5
module load StdEnv/2020
module load gcc/9.3.0
module load gmtk
source segwayEnv/bin/activate

python wrapperForSegway_oneSample.py
