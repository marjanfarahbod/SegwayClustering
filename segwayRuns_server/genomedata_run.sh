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

# running genomedata
genomedata-load -s hg38.chrom.sizes --sizes -t CTCF=ENCFF181ESK.bedgraph -t H3K9me3=ENCFF517GCE.bedgraph -t H3K36me3=ENCFF217BUR.bedgraph --verbose output.genomedata
