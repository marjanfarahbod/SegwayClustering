# wrapper for the genomedata command, that is use the command, one line at a time from a commnad line.
# make the genomedata command, run it with array. Each job for getting the genomedata command needs 20G. Once it is done, the log will print that it is done.

import math
import os
import shutil
import pickle
import glob

#parameterIndex = int(os.environ['SLURM_ARRAY_TASK_ID'])
parameterIndex = 0
print(parameterIndex)

### cedar add
dataFolder = '/home/mfarahbo/scratch/segway112Samples/'

# TODO load the address file
inputFile = dataFolder + 'accessionList.pkl'
with open(inputFile, 'rb') as f:
    accessionList = pickle.load(f)

# TODO get the Segway accession with the parameter index from the list
segwayAccession = accessionList[int(parameterIndex)]

# TODO: get the list of beggraph
segdir = dataFolder + segwayAccession + '/'

# get list of files ending with .BedGraph
bedGraphList = glob.glob('%s*.bedGraph' %(segdir))

# match the tarcks with the accession list
print(bedGraphList)

# make the command
command = 'genomedata-load -s h38.chrom.sizes.txt --sizes'
for bgfile in bedGraphList:
    fileID = bgfile.split('/')[-1].split('.')[0]
    command += ' -t %s=%s' %(fileID, bgfile)

command +=  ' --verbose ' + segdir + 'files.genomedata'

os.system(command)

print('printed genomedata')
