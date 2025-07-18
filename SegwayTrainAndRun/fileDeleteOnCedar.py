
import pickle
import os
import shutil

# load the address file
dataFolder = '/home/mfarahbo/scratch/segway112Samples/'
inputFile = dataFolder + 'accessionList.pkl'
with open(inputFile, 'rb') as f:
    accessionList = pickle.load(f)

# Get the Segway accession with the parameter index from the list

for ann in accessionList:
    print(ann)

    segdir = dataFolder + ann + '/'
    # where the output of Segway training process goes
    trainFolder = segdir + 'call_segway_train/'
    annotateFolder = segdir + 'call_segway_annotate/'
    

    # if annfolder exists:
    folderList = os.listdir(segdir)
    if ('call_segway_train') in folderList:
        print('deleting')
        shutil.rmtree(trainFolder)
        print('first folder')
        shutil.rmtree(annotateFolder)
        print('files deleted')



