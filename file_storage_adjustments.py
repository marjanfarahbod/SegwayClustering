# for each of the .bed files in the folder, zip them and then delete the unzipped one.

import os
import shutil

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
dataSubFolder = 'testBatch105/fromAPI/'


inputFile = dataFolder + dataSubFolder + 'metaInfo.pkl'
with open(inputFile, 'rb') as f:
    annMeta = pickle.load(f)


for annAccession in annMeta.keys():

    print(annAccession)
    annFolder = dataFolder + dataSubFolder + annAccession + '/'

    allFiles = os.listdir(annFolder)

    for file in allFiles:
        if file.endswith('.bed'):
            print(file)
            # zip it
            #outputName = annFolder + file + '.gz'
            inputName = annFolder + file
            os.system('gzip %s' %(inputName))


