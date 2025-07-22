

# get the sample info for the 105 runs
inputFile = dataFolder + dataSubFolder + runIDFile
with open(inputFile, 'rb') as f:
    runIDs = pickle.load(f)


# get the sample info for the 95 filesp

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
dataSubFolder = 'testBatch_May112022/'

folder = dataFolder + dataSubFolder

import os

book = os.listdir(folder)

for i in range(len(book)):
    print(str(i) + ':' + book[i])

May92RunIDs = book
outputFile = dataFolder + dataSubFolder + 'runID_list.pkl'
with open(outputFile, 'wb') as f:
    pickle.dump(May92RunIDs, f)


