# Get some paper data and supplementary figures

import gzip
import linecache
import pickle
import os
import numpy as np
import pandas as pd
import subprocess as sp
import seaborn as sns
import matplotlib.pyplot as plt


# General data folder
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/'

segwayStates = ['Promoter', 'PromoterFlanking', 'Enhancer', 'EnhancerLow', 'Bivalent', 'CTCF', 'Transcribed', 'K9K36', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent']
print(len(segwayStates))
segwayStateCount = len(segwayStates)

inputFileName = 'all235Annot_meta_corrected.pkl'
#inputFileName = 'all235Annot_meta.pkl'
inputFile = dataFolder +  inputFileName
with open(inputFile, 'rb') as f:
    allMeta = pickle.load(f)

accessionList = list(allMeta.keys())
accessionList.remove(accessionList[205])

########################################
# 1. count of track data
# get the address of the GMTK files for each of the samples

totalCount = 0
countPerSample = np.zeros(len(accessionList))
extraTrackCount = np.zeros(5)  # CTCF, DNase-seq, ATAC-seq, POLR2A, EP300
extraTrackList = ['CTCF', 'DNase-seq', 'ATAC-seq', 'POLR2A', 'EP300']
for a, accession in enumerate(accessionList):

    annotation = allMeta[accession]
    annotationFolder = annotation['folder']

    mapping_file = annotationFolder + 'trackname_assay.txt'
    # read the mapping_file
    track_assay_map = {}
    sampleTrack_list = [] # list of the track type for the sample 
    with open(mapping_file) as inputFile:
        for line in inputFile:
            fields = line.strip().split()
            track_assay_map[fields[0]] = fields[1]
            sampleTrack_list.append(fields[1])
            if fields[1] in extraTrackList:
                extraTrackCount[extraTrackList.index(fields[1])] +=1
    countPerSample[a] = len(sampleTrack_list)
    totalCount += len(sampleTrack_list)


# plots
# bar plot for the count of tracks per sample
barData = np.zeros(6)
for c in countPerSample:
    barData[int(c-6)] +=1

plt.bar(range(6), barData, color=[0,0,0,1], width=.75)
plt.title('count of samples with count of tracks 6-11')
figureFolder = 'suppFigures/'
figFile = plotFolder + figureFolder + 'sample_trackCount.pdf'
plt.savefig(figFile)
plt.close('all')
#plt.show()

# bar plot for count of samples with specific tracks
plt.bar(extraTrackList, extraTrackCount, color=[0,0,0,1], width=.75)
plt.title('count of samples with extra tracks')
figureFolder = 'suppFigures/'
figFile = plotFolder + figureFolder + 'track_sampleCount.pdf'
plt.savefig(figFile)
plt.close('all')
# plt.show()




