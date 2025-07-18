#

import pickle

inputFile = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/the112batch/accessionList.pkl'
with open(inputFile, 'rb') as f:
    accessionList = pickle.load(f)

textFile = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/the112batch/accessionList.txt'
with open(textFile, 'w') as f:
    for item in accessionList:
        f.write(item+'\n')

file = dataFolder + dataSubFolder + 'hg_accessionList.pkl'
with open(file, 'rb') as f:
    hg_accessionList = pickle.load(f)

textFile = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/the112batch/hg_accessionList.txt'
with open(textFile, 'w') as f:
    for item in hg_accessionList:
        f.write(item+'\n')

