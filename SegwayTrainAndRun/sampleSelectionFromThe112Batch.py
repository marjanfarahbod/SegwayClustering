# This is to get the accession for human samples from the last 112 batch.

import util # this is for features_from_segtools_dir
import gzip
import pickle
import pandas as pd

# the file

file = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/MetaSheet310Samples.tsv'

df = pd.read_table(file)
for i,header in enumerate(list(df)):
    print('%d %s' %(i, header))

# the 112 IDs
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
dataSubFolder = 'the112batch/'

# load the address file
inputFile = dataFolder + dataSubFolder + 'accessionList.pkl'
with open(inputFile, 'rb') as f:
    accessionList = pickle.load(f)

# getting the rows with the 112 IDs and without the mouse in their columns
headers = list(df)
allAccessions = df[headers[15]] # header with the accession

# getting list of mouse accessions
mouse_accessions = [] 
for i,target in enumerate(df['targets']):
    if 'mouse' in target:
        print(i)
        mouse_accessions.append(allAccessions[i].split('/')[2])

hg_accessionList = []
for accession in accessionList:
    if not(accession in mouse_accessions):
        hg_accessionList.append(accession)

file = dataFolder + dataSubFolder + 'hg_accessionList.pkl'
with open(file, 'wb') as f:
    pickle.dump(hg_accessionList,f)
# save it, get the runs for it!

# how many are in the first 15 that we have already done the processing on
for accession in hg_accessionList:
    if accession in accessionList[0:16]:
        print(accession)

# 5, we have already done the analysis on

        
