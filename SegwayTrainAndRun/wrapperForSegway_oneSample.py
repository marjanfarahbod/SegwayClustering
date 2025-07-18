# wrapper for Segway, that is use the command, one line at a time from a commnad line.
# make the genomedata command, run it with array. Each job for getting the genomedata command needs 20G. Once it is done, the log will print that it is done.

import math
import os
import shutil
import pickle
import glob

#parameterIndex = int(os.environ['SLURM_ARRAY_TASK_ID'])
parameterIndex = 1
print(parameterIndex)

### cedar add
dataFolder = '/home/mfarahbo/scratch/segway112Samples/'

# load the address file
inputFile = dataFolder + 'hg_accessionList.pkl'
with open(inputFile, 'rb') as f:
    accessionList = pickle.load(f)

# Get the Segway accession with the parameter index from the list
segwayAccession = accessionList[int(parameterIndex)]

# sampleFolder
segdir = dataFolder + segwayAccession + '/'

# get genomedata address
genomeDataFile = segdir + 'files.genomedata'

trainFolder = segdir + 'call_segway_train/'
os.system('mkdir %s' %(trainFolder))

outputFolder = segdir + 'call_segway_annotate/'
os.system('mkdir %s' %(outputFolder))


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
postFolder = '%sposterior' %(trainFolder)
postCommand = 'segway posterior %s %s %s' %(genomeDataFile, trainFolder, postFolder)

os.system(postCommand)

# collect results
#os.system('cp %s/segway.bed.gz %s/' %(postFolder, outputFolder))
#os.system('cp %s/log/segway.sh %s/' %(trainFolderAll, outputFolder))
#os.system('cp %s/segway.str %s/' %(trainFolder, outputFolder))
#os.system('cp %s/params/params.params %s' %(trainFolder, outputFolder))

# # # # doing the clean up
# shutil.rmtree(trainFolder)


