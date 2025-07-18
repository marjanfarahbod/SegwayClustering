# This is the code to get the heatmap for the samples for coverage of the labels for a specific region of the genome. The last update is that I only get the heatmap with the right colors, for the given region. I can modify it to use bedtools to retrieve the code (right now I get a subsample bed file for a chromosome, and then go through that, bedtools should be much faster). A missing part from the code is ordering the sampels so that it looks nice, plus if the sampels are more, things like that. 


from urllib.request import urlopen
import urllib
import json
import http.cookiejar
import requests
import re
import subprocess
import os
import pickle
import numpy as np
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'

dataSubFolder = 'testBatch105/fromAPI/'

inputFile = dataFolder + dataSubFolder + 'metaInfo.pkl'
with open(inputFile, 'rb') as f:
    annMeta = pickle.load(f)

sample_count = len(annMeta)

annAccessionList = list(annMeta.keys())
annAccession = annAccessionList[104]
print(annAccession)

inputFileName = 'all_annInfo_list.pkl'
inputFile = dataFolder + dataSubFolder + inputFileName
with open(inputFile, 'rb') as f:
    ann_info_list = pickle.load(f)

# get the mapping from accession to index
accession_index_map = {}
for ann in ann_info_list:
   accession_index_map[ ann['accession'] ] = ann['index']


annAccession = 'ENCSR566HBT'
index = accession_index_map[annAccession]
print(index)


# 1. Get chr19 files
#########################################

#coverage = np.zeros((max_region_count, int(region_length/10)+40))
# create the en file
for annAccession in annAccessionList[7:]:
    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
    print(sampleFolderAdd)

    index = accession_index_map[annAccession]
    print(index)

    # get the .bed file chr19 extract. 
    originalBedFile = annMeta[annAccession]['bedFile']
    originalBed_accession = originalBedFile.split('.')[0]
    segwayFile = dataFolder + dataSubFolder + annAccession + '/' + originalBed_accession + '_filteredSorted.bed.gz'

    os.system('gunzip %s' %(segwayFile))
    command = "grep -E 'chr21' %s" %(segwayFile[0:-3])
    print(command)
    
    out = sampleFolderAdd + 'chr21_whole.bed'
    f = open(out, 'w')
    subprocess.run(command, shell=True, stdout=f)

    os.system('gzip %s' %(segwayFile[0:-3]))

    print('chr21 enhancer regions selected')


s = 38748000
e = 38828000
myMat = np.zeros((105, int((e-s)/100)))
annCount = 0
for annAccession in annAccessionList:

    # fill the matrix
    print(annCount)
    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
    bedFile = sampleFolderAdd + 'chr21_whole.bed'

    with open(bedFile, 'r') as bed:
        for line in bed:
            start = int(line.strip().split()[1])
            end = int(line.strip().split()[2])
            if start > s and end < e:
                label = int(line.strip().split()[3].split('_')[0])
                matIndexStart = int((start - s)/100)
                matIndexEnd = int((end - s)/100)
                myMat[annCount, matIndexStart:matIndexEnd] = label

    annCount +=1


segwayLabels = ['EnhancerLow', 'Enhancer', 'PromoterFlanking', 'Promoter', 'Transcribed', 'CTCF', 'K9K36', 'Bivalent', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent']

label_color_mapping = {}
label_color_mapping[3] = [255, 0, 0] # promoter
label_color_mapping[2] = [200, 10, 0] # promoter
label_color_mapping[1] =  [180, 160, 0] # enhancer
label_color_mapping[0] = [220, 210, 0] # enhancer
label_color_mapping[4] = [0, 170, 0]# transcribed
label_color_mapping[5] =  [78, 111, 249]# CTCF
label_color_mapping[6] =  [159, 17, 193]# K9K36
label_color_mapping[7] = [205, 92, 90] # bivalent
label_color_mapping[8] = [128, 128, 128] # FacHet
label_color_mapping[9] = [16, 196, 191] # consHet
label_color_mapping[10] = [0, 0, 0] # Quis

colormap = label_color_mapping[0]
colorlist = []
colorlist.append(np.array(label_color_mapping[0])/255)
for label in range(1,11):
    colormap = np.vstack([colormap, label_color_mapping[label]])
    colorlist.append(np.array(label_color_mapping[label])/255)

import matplotlib.colors as colors
cmap = colors.ListedColormap(colorlist)
boundaries = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
norm = colors.BoundaryNorm(boundaries, cmap.N, clip=True)

myColorMat = np.zeros((105, int((e-s)/100)))
annIndex = 0
for annAccession in annAccessionList:

    label_term_mapping = {}
    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
    mnemonics_file = sampleFolderAdd + 'mnemonics_v04.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            print(line)
            label = int(line.strip().split()[0])
            term = line.strip().split()[1]
            label_term_mapping[label] = term

    for i, label in enumerate(myMat[annIndex,:]):
        mylabel = int(myMat[annIndex, i])
        myterm = label_term_mapping[mylabel]
        color = segwayLabels.index(myterm)
        myColorMat[annIndex, i] = color
    
    annIndex +=1

import seaborn as sns
figFile = dataFolder + 'ISMBFigure.pdf'
sns.heatmap(myColorMat, cmap=cmap, norm=norm)
plt.savefig(figFile)
import matplotlib.pyplot as plt
plt.show()


