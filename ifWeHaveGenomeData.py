
import pickle
import os
import shutil

# load the address file
dataFolder = '/home/mfarahbo/scratch/segway112Samples/'
inputFile = dataFolder + 'accessionList.pkl'
with open(inputFile, 'rb') as f:
    accessionList = pickle.load(f)

# Get the Segway accession with the parameter index from the list

c = 0
for ann in accessionList:
    segdir = dataFolder + ann + '/'
    
    # if annfolder exists:
    folderList = os.listdir(segdir)
    if not(('files.genomedata') in folderList):
        print('file is not')
        print(ann)
        print(c)

    
