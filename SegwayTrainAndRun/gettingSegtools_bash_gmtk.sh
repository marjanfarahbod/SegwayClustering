#!/bin/bash
#SBATCH --job-name=rep1_mlttrk
#SBATCH --account=rrg-maxwl
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=0-00:30
#SBATCH --output=%j_segtools_gmtk_5_37.out
#SBATCH --array=5-37

module load python/3.8.10
module load hdf5
module load StdEnv/2020
module load gcc/9.3.0 r/4.0.2
module load gmtk
source /home/mfarahbo/scratch/segwayENV/bin/activate

python wrapperForSegtools_gmtk.py
