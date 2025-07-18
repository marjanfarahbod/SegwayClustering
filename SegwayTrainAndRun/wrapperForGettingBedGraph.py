# running from cedar scratch
import math
import os
import shutil
import glob
import pickle


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

# TODO: get the list of bigwig files in the directory

segdir = dataFolder + segwayAccession + '/'

# get list of files ending with .bigwig

bigwigList = glob.glob('%s*.bigWig' %(segdir))

print(bigwigList)

# for each file, run the bedgraph comment
os.chdir('/home/mfarahbo/scratch/')
for bigwigFile in bigwigList:
    fileID = bigwigFile.split('.')[0]
    command = './bigWigToBedGraph %s %s.bedGraph' %(bigwigFile, fileID)
    os.system(command)
    print(command)

print(segwayAccession)
    

