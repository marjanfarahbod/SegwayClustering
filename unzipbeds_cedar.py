# for a given list of samples, run the gzip command for the bed files.
import os
import pickle

dataFolder = '/home/mfarahbo/projects/def-maxwl/mfarahbo/segway112Samples/'

# load the address file
inputFile = dataFolder + 'hg_accessionList.pkl'
with open(inputFile, 'rb') as f:
    accessionList = pickle.load(f)

for accession in accessionList:
    print(accession)
    bedgzFile = dataFolder + accession + '/segOutput/segway.bed.gz'

    command = 'gunzip %s' %(bedgzFile)
    print(command)
    os.system(command)
