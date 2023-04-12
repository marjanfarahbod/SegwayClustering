# wrapper for Segway, that is use the command, one line at a time from a commnad line.
# make the genomedata command, run it with array. Each job for getting the genomedata command needs 20G. Once it is done, the log will print that it is done.

import math
import os
import shutil
import pickle
import glob
from pathlib import Path

parameterIndex = int(os.environ['SLURM_ARRAY_TASK_ID'])
#parameterIndex = 3
print(parameterIndex)

### cedar add TODO: change to project address - DONE, but nevermind, project has not much space
# dataFolder = '/home/mfarahbo/projects/def-maxwl/mfarahbo/segway112Samples/'
dataFolder = '/home/mfarahbo/scratch/segway112Samples/'

# load the address file
inputFile = dataFolder + 'accessionList.pkl'
with open(inputFile, 'rb') as f:
    accessionList = pickle.load(f)

# Get the Segway accession with the parameter index from the list
segwayAccession = accessionList[int(parameterIndex)]

# sampleFolder
segdir = dataFolder + segwayAccession + '/'
previousBed = Path(segdir + 'segOutput/segway.bed.gz')

### TODO here, check to see if the segway output folder contains the .bed file, we skip the sample, otherwise we go on - DONE
print(segdir)
if previousBed.is_file():
    print('bed file exists')
    exit()

# get genomedata address
genomeDataFile = segdir + 'files.genomedata'

# where the output of Segway training process goes
# TODO: here we check, if the folder exists, we delete it, we make it again
trainFolder = segdir + 'call_segway_train/'
if os.path.isdir(trainFolder):
    print('deleting previous train folders')
    shutil.rmtree(trainFolder)

# making new train folder
os.system('mkdir %s' %(trainFolder))

# where the output of Segway annotation process goes
# TODO: here we check, if the folder exists, we delete it, we make it again
annotateFolder = segdir + 'call_segway_annotate/'
if os.path.isdir(annotateFolder):
    print('deleting previous annot folder')
    shutil.rmtree(annotateFolder)
    
os.system('mkdir %s' %(annotateFolder))

segwayOutputFolder = segdir + 'segOutput/'
os.system('mkdir %s' %(segwayOutputFolder))

# get the label count
trackFile = segdir + 'trackname_assay.txt'
trackCount = 0
with open(trackFile, 'r') as f:
    for line in f:
        trackCount +=1

labelCount = 10 + 2* int(math.sqrt(trackCount))
trainCommand = 'segway train --max-train-rounds=25 --num-instances=10 --track-weight=0.01 --segtransition-weight-scale=1 --ruler-scale=100 --prior-strength=1 --resolution=100 --minibatch-fraction=0.01 --num-labels=%d %s %s' %(labelCount, genomeDataFile, trainFolder)

os.system(trainCommand)

# segway posterior command
postCommand = 'segway posterior %s %s %s' %(genomeDataFile, trainFolder, annotateFolder)

os.system(postCommand)

# collect results
os.system('cp %s/segway.bed.gz %s/' %(annotateFolder, segwayOutputFolder))
os.system('cp %s/log/segway.sh %s/' %(trainFolder, segwayOutputFolder))
os.system('cp %s/segway.str %s/' %(trainFolder, segwayOutputFolder))
os.system('cp %s/params/params.params %s' %(trainFolder, segwayOutputFolder))
#os.system('cp -r %s/triangulation %s') %(trainFolder, segwayOutputFolder)
#os.system('cp -r %s/auxiliary %s') %(trainFolder, segwayOutputFolder)
#os.system('cp -r %s/likelihood %s') %(trainFolder, segwayOutputFolder)

# clean up
shutil.rmtree(trainFolder)
shutil.rmtree(annotateFolder)

print(segwayAccession)

