# Plots for paper Figure 1
# 0. Initials
# 1. The heatmap for assay presence DONE
# 2. The actual annotation for samples sorted - I donno what this is
# 3. The long heatmap with all things sorted for all the signal. For genes without transcriptomic data, just the genomic regions(supplement). (so I want each "interpretation term" to be clustered, then them presented in a long heatmap)
# 4. The prediction probabilities DONE
# *. the big fat genome regional plot
# *. The prediction some sort of filtering?

#
#########################################
# 0. Initials
#########################################
from numpy import random
import pickle
import re
import numpy as np
import pandas as pd
import subprocess as sp
import seaborn as sns
import matplotlib.pyplot as plt

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/'
figureFolder = 'figure01/'

# load the meta file 
inputFileName = 'all235Annot_meta_corrected.pkl'
#inputFileName = 'all235Annot_meta.pkl'
inputFile = dataFolder +  inputFileName
with open(inputFile, 'rb') as f:
    allMeta = pickle.load(f)

# Segway states:
segwayStates = ['Promoter', 'PromoterFlanking', 'Enhancer', 'EnhancerLow', 'Bivalent', 'CTCF', 'Transcribed', 'K9K36', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent']
print(len(segwayStates))
segwayStateCount = len(segwayStates)

trackList = ['H3K36me3', 'H3K27me3', 'H3K4me3', 'H3K9me3', 'H3K27ac', 'H3K4me1', 'CTCF', 'DNase-seq', 'ATAC-seq', 'POLR2A', 'EP300']

########################################
# 1. The heatmap for assay presence 
########################################

'''
open up the matrix of 11 tracks 
'''
trackCount = len(trackList)
sampleCount = len(allMeta)
trackMat = np.zeros((sampleCount, trackCount))

accessionList = list(allMeta.keys())

for i,accession in enumerate(accessionList):

    print(accession)
    print(i)
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']

    # get the track file
    assay_file = annotationFolder + 'trackname_assay.txt'
    with open(assay_file, 'r') as assays:
        for line in assays:
            assay = line.strip().split()[1]
            if assay in trackList:
                index = trackList.index(assay)
                trackMat[i, index] = 1

sns.heatmap(trackMat)
plt.show()

book = np.sum(trackMat, axis=0)
kado = np.sum(trackMat, axis=1)

sortedInds = np.argsort(kado)[::-1]

sortedTrackMat = trackMat[sortedInds,:]
cmap = sns.color_palette("rocket_r", as_cmap=True)
 
sns.heatmap(sortedTrackMat, xticklabels=trackList, cmap=cmap)

plotFolder_add = plotFolder + figureFolder
figFile = plotFolder_add + 'heatmap_tracksForSamples.pdf'
print(figFile)
plt.tight_layout()
plt.savefig(figFile)
plt.close('all')

########################################
# 4. The prediction probabilities, filter
########################################
# this was done in code prob_filter.py

# 4.0 the boxplot for distributions
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
accessionList = list(allMeta.keys())
activeLabels = segwayLabels[0:5]

medAll = np.zeros(len(accessionList))
medAllActive = np.zeros(len(accessionList))
medAllActiveSum = np.zeros(len(accessionList))
enhancerProb = np.zeros(len(accessionList))
for ai,accession in enumerate(accessionList):

    if accession == 'ENCSR424JDX':
        continue

    print(accession)
    print(ai)
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']

    # the mnemonics file
    label_term_mapping = {}
    mnemonics_file = annotationFolder + 'mnemonics_v04.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term

    labelCount = len(label_term_mapping)

    # get the track file
    prob_file = annotationFolder + 'probs_v04.csv'
    df = pd.read_csv(prob_file)

    # fetch all the probs from the df
    allProbs = np.zeros(labelCount)
    for i in range(len(label_term_mapping)): # for each label
        l = str(i)
        term = label_term_mapping[l] # to find the column of the df
        allProbs[i] = df._get_value(i, term)
        

    # fetch all the from the active labels (Enhancer, EnhancerLow, Promoter, PromoterFlank, Transcribed)
    activeProbs = []
    for i in range(labelCount):
        term = label_term_mapping[str(i)]
        if term in activeLabels:
            activeProbs.append(allProbs[i])
        

    # fetch the active probs from the df, but sum up the probs of same group
    activeProbsSum = []
    enhancerSum = []
    for i in range(labelCount):
        term = label_term_mapping[str(i)]
        if term in activeLabels:
            print(term)
            if 'Enhancer' in term:
                sumProb = df._get_value(i, 'Enhancer') + df._get_value(i, 'EnhancerLow')
                print(df._get_value(i, 'Enhancer'))
                print(df._get_value(i, 'EnhancerLow'))
                print(sumProb)
                activeProbsSum.append(sumProb)
                enhancerSum.append(sumProb)
            if 'Promoter' in term:
                sumProb = df._get_value(i, 'Promoter') + df._get_value(i, 'PromoterFlanking')
                print(sumProb)
                activeProbsSum.append(sumProb)
            if not('Promoter' in term) and not('Enhancer' in term):
                activeProbsSum.append(allProbs[i])


    medAll[ai] = np.median(allProbs)
    medAllActive[ai] = np.median(activeProbs)
    medAllActiveSum[ai] = np.median(activeProbsSum)
    if len(enhancerSum)>0:
        enhancerProb[ai] = np.max(enhancerSum)

# the boxplot of the probabilities

plt.boxplot([medAll, medAllActive, medAllActiveSum, enhancerProb])
plt.show()
plotFolder_add = plotFolder + figureFolder
figFile = plotFolder_add + 'boxplot_probsMed.pdf'
print(figFile)
plt.tight_layout()
plt.savefig(figFile)
plt.close('all')

# 4.1 The distribution plot
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# this plot is for each term, the distribution of the probabilities for all other terms
# make a dictionary of tupples where each tuple is the index of the matrix, for which we have a list of values


probKeyList = []
for i in range(len(segwayStates)):
    for j in range(len(segwayStates)):
        probKeyList.append((i,j))

probDistDict = {}
for key in probKeyList:
    probDistDict[key] = []

for ai,accession in enumerate(accessionList):
 
    if accession == 'ENCSR424JDX':
        continue

    print(accession)
    print(ai)
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']

    # the mnemonics file
    label_term_mapping = {}
    mnemonics_file = annotationFolder + 'mnemonics_v04.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term

    labelCount = len(label_term_mapping)

    # get the track file
    prob_file = annotationFolder + 'probs_v04.csv'
    df = pd.read_csv(prob_file)

    for i in range(labelCount):
        # row: is the highest probability
        rowInd = segwayStates.index(label_term_mapping[str(i)])
        #thisTerm = label_term_mapping[str(i)]
        print(i)
        print(rowInd)
        #print(thisTerm)
        #print(df._get_value(i, thisTerm))
        #print('-----')

        # colomn: to the segway label
        for term in segwayStates:
            myProb = df._get_value(i, term)
            #print(term)
            #print(myProb)
            colInd = segwayStates.index(term)
            probDistDict[(rowInd, colInd)].append(myProb)

from matplotlib import colors
fig, axs = plt.subplots(11, 11, figsize=(12,12))
 
#cmap = sns.diverging_palette(245, 15, s=80, l=40, sep=2, n = 15, as_cmap=True, center='dark')
#cmap = sns.diverging_palette(245, 15, s=80, l=40, sep=2, n = 15, as_cmap=True, center='dark')
cmap = plt.cm.get_cmap('Oranges')

bgcolors = []
for i in range(120):
    bgcolors.append(cmap(i*5))

binCount = 15
barx = np.zeros(binCount)
for k in range(binCount):
    barx[k]= np.mean([bins[k], bins[k+1]])

meanBinMat = np.zeros((11,11))
for i in range(11):
    for j in range(11):
        print(i, j)
        #n, bins, patches = plt.hist(probDistDict[(i,j)], bins=binCount, density=True)
        bins = np.linspace(0, 1, 16)
        hist = np.histogram(probDistDict[(i,j)], bins, density=True)

        #meanBin = np.mean(hist[0])
        sib = [hist[0][i]*hist[1][i] for i in range(len(hist[0]))]
        meanBin = np.mean(sib)
        meanBinMat[i,j] = meanBin

        color = bgcolors[int(meanBin*100)]

        axs[i,j].set_facecolor(color)
        axs[i,j].bar(barx, hist[0], width=.06, color='black')
        #axs[i,j].bar(barx, n, width=.06, color='black')
        #axs[i,j].hist(probDistDict[(i,j)], bins=15, density=True, color=['blue', 'green', 'black'])
        axs[i,j].set_xlim(-.1,1.1)
        axs[i,j].set_ylim(0,8)
        #axs[i,j].set_xticks([0, 1])
        #axs[i,j].set_yticks([0, 4, 8])
        axs[i,j].set_xticks([])
        axs[i,j].set_yticks([])
        axs[i,j].set_yticklabels([])
        axs[i,j].set_xticklabels([])

#plt.show()
plotFolder_add = plotFolder + figureFolder
figFile = plotFolder_add + 'histograms_probsSqure.pdf'
print(figFile)
plt.tight_layout()
plt.savefig(figFile)
plt.close('all')

sns.heatmap(meanBinMat, cmap=cmap)
plotFolder_add = plotFolder + figureFolder
figFile = plotFolder_add + 'histograms_probsSqure_heatmapForColorBar.pdf'
print(figFile)
plt.tight_layout()
plt.savefig(figFile)
plt.close('all')

plt.show()
            
# do the histogram plot
book = probDistDict[(0,0)]

plt.hist(book, bins=15, density=True)
plt.bar(barx, n, width=.06, color=colors)
plt.show()

########################################
# *. the big fat genome regional plot
########################################
# Look around active genes - pull up some papers. Start with the same region.

# define the matrix for collecting
annotMat = np.zeros((234, 20000))
# for each of the accessions
for accession in accessionList:


    annotation = allMeta[accession]
    annFile = annotation['bedFile']

    if annFile.endswith('gz'):
        os.system('gunzip %s' %(annFile))
        annFile = annFile[0:-3]
    else:
        os.system('gunzip %s.gz' %(annFile))

    annLineCount = int(sp.getoutput('wc -l %s' %(annFile)).split()[0])

    line = linecache.getline(annFile, )
    

# pull up the .bed file.

# in the .bed file, reach out to the region

chr15, 73065000 + 20k



