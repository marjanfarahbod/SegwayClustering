#!/bin/bash
#SBATCH --job-name=rep1_mlttrk
#SBATCH --account=rrg-maxwl
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --time=0-24:00
#SBATCH --output=%j.out

source segtools-env2/bin/activate
module load python/2.7.18
module load r/4.2.1
module load hdf5/1.10.6

python wrapperForSegtools_oneSample.py
