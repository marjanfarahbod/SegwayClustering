# wrapper for Segtools, that is use the command, one line at a time from a commnad line.
# make the genomedata command, run it with array. Each job for getting the genomedata command needs 20G. Once it is done, the log will print that it is done.

import math
import os
import shutil
import pickle
import glob

#parameterIndex = int(os.environ['SLURM_ARRAY_TASK_ID'])
#parameterIndex = 12
#print(parameterIndex)

### cedar add
dataFolder = '/home/mfarahbo/scratch/segway112Samples/'

## commenting the below items for now since they don't work in python2
## load the address file
#inputFile = dataFolder + 'accessionList.pkl'
#with open(inputFile, 'rb') as f:
#    accessionList = pickle.load(f)

# Get the Segway accession with the parameter index from the list
#segwayAccession = accessionList[int(parameterIndex)]

segwayAccession = 'ENCSR020PJR'

# sampleFolder
segdir = dataFolder + segwayAccession + '/'

outputFolder = segdir + 'segOutput/'
genomeDataFile = segdir + 'files.genomedata'

# unzip the .bed file
#command = 'gunzip %ssegway.bed.gz %ssegway.bed' %(outputFolder, outputFolder)
#os.system(command)

bedFile = outputFolder + 'segway.bed'

#----- segtools signal-distribution
#'segtools-signal-distribution {} {} --outdir={}'.format(segbed, gd, outdir))
#signalDistFolder = '%scall_signal_distribution/' %(outputFolder)
#os.mkdir(signalDistFolder)
#command = 'segtools-signal-distribution %s %s --outdir=%s' %(bedFile, genomeDataFile, signalDistFolder)
#os.system(command)

#print('sigdist done')

#----- segtools-aggregation command
# 'segtools-aggregation --normalize --mode=gene {} {} --outdir={}'.format(segbed, gtf, featureAggFolder)
featureAggFolder = '%scall_feature_aggregation/' %(outputFolder)
os.mkdir(featureAggFolder)
gtfFile = '/home/mfarahbo/scratch/gencode.v29.primary_assembly.annotation_UCSC_names.gtf'
command = 'segtools-aggregation --normalize --mode=gene --flank-bases=10000 %s %s --outdir=%s' %(bedFile, gtfFile, featureAggFolder)
os.system(command)

#signal distribution: bed file and genomedata
#aggregation: bed file and GTF file

# chdir to the folder containing the posterior folder
os.chdir(outputFolder)

#----- segtools-length-distribution command
#os.system('segtools-length-distribution segway.bed')

#----- segtools-gmtk-parameters command
os.system('segtools-gmtk-parameters params.params')

# Draft
########################################
# segtools-gmtk-parameters command
#commnad = 'segtools-gmtk-parameters %s/params.params' %(outputFolder)
#os.system(command)
