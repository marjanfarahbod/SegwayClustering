#!/bin/bash
#SBATCH --job-name=rep1_mlttrk
#SBATCH --account=rrg-maxwl
#SBATCH --cpus-per-task=5
#SBATCH --mem=5G
#SBATCH --time=0-03:10
#SBATCH --output=rep1_batch.out

# setting up environment
module load python/3.8.10
module load hdf5
module load StdEnv/2020
module load gcc/9.3.0
module load gmtk
source /home/mfarahbo/segwayENV/bin/activate

# running segway
segway train --include-coords=data/encodePilotRegions.hg19.bed --num-labels=15 --resolution=100 --minibatch-fraction=0.01 --num-instances=10 --prior-strength=1.0 --segtransition-weight-scale=1.0 --ruler-scale=100 --track-weight=0.002 --max-train-rounds=25 data/adrenal_gland.genomedata ag_train
segway posterior --include-coords=data/encodePilotRegions.hg19.bed data/adrenal_gland.genomedata ag_train ag_posterior
gunzip ag_posterior/segway.bed.gz 
segtools-length-distribution ag_posterior/segway.bed
segtools-gmtk-parameters ag_train/params/params.params
