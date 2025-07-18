#!/bin/bash
#SBATCH --job-name=rep1_mlttrk
#SBATCH --account=rrg-maxwl
#SBATCH --cpus-per-task=2
#SBATCH --mem=180G
#SBATCH --time=0-08:00
#SBATCH --output=%j.out
#SBATCH --array=0-111

module unload
module load python/3.7.9
module load hdf5/1.10.6
source gdarchive-env/bin/activate

python wrapperForGenomeData.py
