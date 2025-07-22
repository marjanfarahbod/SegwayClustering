# 0. Initials
# 1. Three panel transcriptomic average
# 2. The ROC curves for the AUC values


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
import linecache

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/'
figureFolder = 'figure02/'

# load the meta file 
inputFileName = 'all235Annot_meta_corrected.pkl'
#inputFileName = 'all235Annot_meta.pkl'
inputFile = dataFolder +  inputFileName
with open(inputFile, 'rb') as f:
    allMeta = pickle.load(f)

accessionList = list(allMeta.keys())
accessionList.remove(accessionList[205])
# Segway states:
segwayStates = ['Promoter', 'PromoterFlanking', 'Enhancer', 'EnhancerLow', 'Bivalent', 'CTCF', 'Transcribed', 'K9K36', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent']
print(len(segwayStates))
segwayStateCount = len(segwayStates)

trackList = ['H3K36me3', 'H3K27me3', 'H3K4me3', 'H3K9me3', 'H3K27ac', 'H3K4me1', 'CTCF', 'DNase-seq', 'ATAC-seq', 'POLR2A', 'EP300']

# get the address of the GMTK files for each of the samples
runID_file = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/testBatch105/' + 'runID_accession_map_105run.pkl'
with open(runID_file, 'rb') as pickledFile:
    runID_accession105 = pickle.load(pickledFile)

accession_runID = {}
for runID in list(runID_accession105.keys()):
    ac = runID_accession105[runID]
    accession_runID[ac] = runID


length_files = {}
for accession in accessionList:
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']
    if accession == 'ENCSR121RCG':
        continue

    if '38batch' in annotationFolder:
        length_files[accession] = annotationFolder + 'segOutput/length_distribution/segment_sizes.tab'


    if 'May11' in annotationFolder:
        length_files[accession] = annotationFolder + 'call-segtools/segment_sizes.tab'

    if 'Batch105' in annotationFolder:
        runID = accession_runID[accession]
        length_files[accession] = dataFolder + 'testBatch105/all_segtools/' + runID +  '/length_distribution/segment_sizes.tab'


########################################
# 1. Three panel transcriptomic average
########################################

# first: sum the column to the labels for each sample and close the file

count = 0
for accession in accessionList:

    annotation = allMeta[accession]
    #if ((('38batch' in annotation['folder']) or ('May11' in annotation['folder'])) and not(annotation['RNAseqFile'] == 'none')):
    if((annotation['RNAseqFile'] != 'none')):
        count +=1
        annotationFolder = annotation['folder']
        expFile = annotationFolder + 'defaultExp_5kg_expSummary_newSort_Q30_regen.pkl'
        print(count, expFile)
        with open(expFile, 'rb') as file:
            expMats = pickle.load(file)

        # load the mnemonics 

        # 1.1 get the mnemonics
        # >>>>>>>>>>>>>>>>>>>>>>>>>>
        label_term_mapping = {}
        mnemonics_file = annotationFolder + 'mnemonics_v04.txt'
        with open(mnemonics_file, 'r') as mnemonics:
            for line in mnemonics:
                print(line)
                label = line.strip().split()[0]
                term = line.strip().split()[1]
                label_term_mapping[label] = term
        labelCount = len(label_term_mapping)

        newMatList = []
        for mat in expMats['clusterMats']:
            newMat = np.zeros((segwayStateCount, mat.shape[1]))
            for label in range(labelCount):
                termInd = segwayStates.index(label_term_mapping[str(label)])
                newMat[termInd, :] += mat[label,:]

            newMatList.append(newMat)

        expTermFile = annotationFolder + 'defaultExp_5kg_expSummary_newSort_Q30_toTerms.pkl'
        print(expTermFile)
        with open(expTermFile, 'wb') as f:
            pickle.dump(newMatList, f)

########################################
# get the intermediate data for the matrix
########################################
            
rnaAccessionList = []
for accession in accessionList:

    annotation = allMeta[accession]
    #if ((('38batch' in annotation['folder']) or ('May11' in annotation['folder'])) and not(annotation['RNAseqFile'] == 'none')):
    if((annotation['RNAseqFile'] != 'none')):
        rnaAccessionList.append(accession)

allBarValues = np.zeros((94, 3, 11, 160)) # a, accession; s, segwayStates; m, expMat; b, bars: log(obs/expected) [note:fractions]
allObsMats = np.zeros((94, 3, 11, 160)) # a, accession; s, segwayStates; m, expMat; b, bars: obs [note:fractions]
geneCountBars = np.zeros((94, 3, 11, 160)) # a, accession; s, segwayStates; m, expMat; b, bars: obs [just the count]
barPresence =  np.zeros((94, 11)) # some samples don't have the some labels
for a, accession in enumerate(rnaAccessionList):
    print(a)

    annotation = allMeta[accession]

    annotationFolder = annotation['folder']
    expFile = annotationFolder + 'defaultExp_5kg_expSummary_newSort_Q30_toTerms.pkl'
    print(expFile)
    with open(expFile, 'rb') as file:
        expMats = pickle.load(file)

    # 1.1 get the mnemonics - we need this fraction coverage
    # >>>>>>>>>>>>>>>>>>>>>>>>>>
    label_term_mapping = {}
    mnemonics_file = annotationFolder + 'mnemonics_v04.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term
            
    labelCount = len(label_term_mapping)

    # 1.2 get the fraction of coverage for each basepair
    # >>>>>>>>>>>>>>>>>>>>>>>>>>
    lfile = length_files[accession]
    fractions = np.zeros(segwayStateCount)
    with open(lfile, 'r') as inputFile:
        header = inputFile.readline()
        header = inputFile.readline()
        for line in inputFile:
            fields = line.strip().split()
            label = fields[0]
            stateInd = segwayStates.index(label_term_mapping[str(label)])
            coverage = fields[-1]
            #print(label, coverage)
            fractions[stateInd] += float(coverage)

    barPresence[a, fractions>0] = 1  # filling the presence
    fractions[fractions == 0] = .0001 # for the calculations only

    # fig, axs = plt.subplots(segwayStateCount, 3, figsize=(12,8))
                
    # xticks = [30, 130]
    # xticksLabels = ['TSS', 'TTS']

    indexList = np.array(list(range(160)))
    for m in range(3):

        thisMat = expMats[m]
        
        geneCountBars[a, m, :, :] = np.copy(thisMat)
        
        # make it the fraction
        thisMat = thisMat / np.sum(thisMat, axis=0)[np.newaxis,:]
        # just saving the observed here (thisMat so far is the observed)
        observedMat = np.copy(thisMat)
        allObsMats[a, m, :, :] = observedMat
        
        # versus expected
        thisMat = thisMat / fractions[:, None]
        logMat = np.log10(thisMat, out=np.zeros_like(thisMat), where=(thisMat!=0))
        allBarValues[a, m, :, :] = logMat

freqRatios = np.zeros((94, 11, 160))
# for the allObsMat, get the high to low ratio
for a in range(len(rnaAccessionList)):
    for s in range(11):
        if(barPresence[a,s] == 1):
            freqRatios[a , s, :] = np.divide(allObsMats[a, 2, s, : ], allObsMats[a, 0, s, :])

meanFRQ = np.zeros((11, 160)) # m, expMat; s, segwayState; b, bars, this is the mean for the fractions

for s in range(11): # for each label
    thisBarPresence = barPresence[:, s]
    for b in range(160): # for each bar
        barValues = freqRatios[:, s, b]
        values = barValues[thisBarPresence ==1]
        meanFRQ[s,b] = np.mean(values)
        
# now for each bar, we get the values that are present, and then we get the median, 75 and 25 percentile.
# allBarValues[a, m, s, b]
# barPresence[a, s]

medianBars = np.zeros((3, 11, 160)) # m, expMat; s, segwayState; b, bars
firstQBars = np.zeros((3, 11, 160)) # m, expMat; s, segwayState; b, bars
thirdQBars = np.zeros((3, 11, 160)) # m, expMat; s, segwayState; b, bars
meanBars = np.zeros((3, 11, 160)) # m, expMat; s, segwayState; b, bars
meanFRQ = np.zeros((1, 11, 160)) # m, expMat; s, segwayState; b, bars, this is the mean for the fractions

for m in range(3): # for each matrix
    for s in range(11): # for each label
        thisBarPresence = barPresence[:, s]
        for b in range(160): # for each bar
            barValues = allBarValues[:, m, s, b]
            values = barValues[thisBarPresence ==1]
            medianBars[m, s, b] = np.quantile(values, .5)
            firstQBars[m, s, b] = np.quantile(values, .25)
            thirdQBars[m, s, b] = np.quantile(values, .75)
            meanBars[m,s,b] = np.mean(values)
        
# the first figure is the transcriptomic figure
# now the plot, first the median bars

### for median and qartiles
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
fig, axs = plt.subplots(segwayStateCount, 3, figsize=(12,8))

xticks = [30, 130]
xticksLabels = ['TSS', 'TTS']

indexList = np.array(list(range(160)))

#thisMat = expMats['clusterMats'][0]

# make it the ratio
#thisMat = thisMat / np.sum(thisMat, axis=0)[np.newaxis,:]
# versus expected
#thisMat = thisMat / fractions[:, None]
#logMat = np.log10(thisMat, out=np.zeros_like(thisMat), where=(thisMat!=0))
for m in range(3):
    thisMedBars = medianBars[m, :, :]
    thisQ1Bars = firstQBars[m, :, :]
    thisQ3Bars = thirdQBars[m, :, :]
    for s in range(segwayStateCount):
        positiveInds = indexList[thisMedBars[s,:] >= 0]
        negativeInds = indexList[thisMedBars[s,:] < 0]
        posBar = np.copy(thisMedBars[s, :])
        posBar[negativeInds]=0
        axs[s,m].bar(range(160), posBar, color=[249/256, 80/256, 97/256,1], width=1)
        negBar = np.copy(thisMedBars[s, :])
        negBar[positiveInds]=0
        axs[s,m].bar(range(160), negBar, color=[50/256,164/256,249/256,1], width=1)
        axs[s,m].set_ylim((-2,2))
        ylabel = segwayStates[s]
        #axs[i,0].text(60, .55, ylabel, fontsize=8)
        axs[s,m].set_xticks(xticks)
        axs[s,m].set_yticks([-1, 1])
        
        #axs[s,m].plot(indexList, thisQ1Bars[s, :], 'k')
        #axs[s,m].plot(indexList, thisQ3Bars[s, :], 'k')

        if m == 0:
            axs[s,m].set_yticklabels([-1, 1], fontsize=8)
            axs[s,m].set_ylabel(ylabel, rotation=0, fontsize=10, labelpad=3, ha='right', va='center')

    axs[segwayStateCount-1,m].set_xticklabels(xticksLabels)
    axs[segwayStateCount-1,m].set_xlabel('Position relative to gene')


axs[0,0].set_title('Enrichment of labels at genes with \nzero expression (log10(observed/expected))', fontsize = 9)
axs[0,1].set_title('Enrichment of labels at genes with \nbottom 30% expression (log10(observed/expected))', fontsize = 9)
axs[0,2].set_title('Enrichment of labels at genes with \ntop 70% expression (log10(observed/expected))', fontsize = 9)

#figFile = plotFolder + figureFolder + 'transcriptomicPlot_median_QartileLines.pdf'
figFile = plotFolder + figureFolder + 'transcriptomicPlot_median.pdf'
plt.savefig(figFile)
plt.close('all')

### for means
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
fig, axs = plt.subplots(segwayStateCount, 3, figsize=(12,8))

xticks = [30, 130]
xticksLabels = ['TSS', 'TTS']

indexList = np.array(list(range(160)))

#thisMat = expMats['clusterMats'][0]

# make it the ratio
#thisMat = thisMat / np.sum(thisMat, axis=0)[np.newaxis,:]
# versus expected
#thisMat = thisMat / fractions[:, None]
#logMat = np.log10(thisMat, out=np.zeros_like(thisMat), where=(thisMat!=0))
for m in range(3):
    thisMeanBars = meanBars[m, :, :]
    for s in range(segwayStateCount):
        positiveInds = indexList[thisMeanBars[s,:] >= 0]
        negativeInds = indexList[thisMeanBars[s,:] < 0]
        posBar = np.copy(thisMeanBars[s, :])
        posBar[negativeInds]=0
        axs[s,m].bar(range(160), posBar, color=[249/256, 80/256, 97/256,1], width=1)
        negBar = np.copy(thisMeanBars[s, :])
        negBar[positiveInds]=0
        axs[s,m].bar(range(160), negBar, color=[50/256,164/256,249/256,1], width=1)
        axs[s,m].set_ylim((-2,2))
        ylabel = segwayStates[s]
        #axs[i,0].text(60, .55, ylabel, fontsize=8)
        axs[s,m].set_xticks(xticks)
        axs[s,m].set_yticks([-1, 1])
        
        if m == 0:
            axs[s,m].set_yticklabels([-1, 1], fontsize=8)
            axs[s,m].set_ylabel(ylabel, rotation=0, fontsize=10, labelpad=3, ha='right', va='center')

    axs[segwayStateCount-1,m].set_xticklabels(xticksLabels)
    axs[segwayStateCount-1,m].set_xlabel('Position relative to gene')


axs[0,0].set_title('Enrichment of labels at genes with \nzero expression (log10(observed/expected))', fontsize = 9)
axs[0,1].set_title('Enrichment of labels at genes with \nbottom 30% expression (log10(observed/expected))', fontsize = 9)
axs[0,2].set_title('Enrichment of labels at genes with \ntop 70% expression (log10(observed/expected))', fontsize = 9)

#figFile = plotFolder + figureFolder + 'transcriptomicPlot_median_QartileLines.pdf'
figFile = plotFolder + figureFolder + 'transcriptomicPlot_mean.pdf'
plt.savefig(figFile)
plt.close('all')


# the plot for the fraction fraction: meanFRQ
# >>>>>>>>>>>>>>>>>>>>
fig, axs = plt.subplots(segwayStateCount, 1, figsize=(4,8))

xticks = [30, 130]
xticksLabels = ['TSS', 'TTS']

indexList = np.array(list(range(160)))

thisMeanBars = np.log(meanFRQ)
for s in range(segwayStateCount):
    positiveInds = indexList[thisMeanBars[s,:] >= 0]
    negativeInds = indexList[thisMeanBars[s,:] < 0]
    posBar = np.copy(thisMeanBars[s, :])
    posBar[negativeInds]=0
    axs[s].bar(range(160), posBar, color=[249/256, 80/256, 97/256,1], width=1)
    negBar = np.copy(thisMeanBars[s, :])
    negBar[positiveInds]=0
    axs[s].bar(range(160), negBar, color=[50/256,164/256,249/256,1], width=1)
    axs[s].set_ylim((-3,5))
    ylabel = segwayStates[s]
    #axs[i,0].text(60, .55, ylabel, fontsize=8)
    axs[s].set_xticks(xticks)
    axs[s].set_yticks([-2, 4])
    axs[s].set_ylabel(ylabel, rotation=0, fontsize=8, labelpad=3, ha='right', va='center')
    if s < 10:
        axs[s].tick_params(axis='x', labelbottom=False)

#plt.show()
axs[segwayStateCount-1].set_xticklabels(xticksLabels)
axs[segwayStateCount-1].set_xlabel('Position relative to gene')
#plt.figure(constrained_layout=True)
plt.tight_layout()
plt.show()


axs[0,0].set_title('Enrichment of labels at genes with \nzero expression (log10(observed/expected))', fontsize = 9)
axs[0,1].set_title('Enrichment of labels at genes with \nbottom 30% expression (log10(observed/expected))', fontsize = 9)
axs[0,2].set_title('Enrichment of labels at genes with \ntop 70% expression (log10(observed/expected))', fontsize = 9)

#figFile = plotFolder + figureFolder + 'transcriptomicPlot_median_QartileLines.pdf'
figFile = plotFolder + figureFolder + 'transcriptomicPlot_mean_freqThing.pdf'
plt.savefig(figFile)
plt.close('all')

########################################
# rerunning that thing because I messed up
########################################


from transcription_overlap import SegwayTranscriptionEnrichment

# GTF data structure
fileName = dataFolder + '/geneLists.pkl'
with open(fileName, 'rb') as pickleFile:
    geneListsAndIDs = pickle.load(pickleFile)
    
geneList = geneListsAndIDs[0]
geneIDList = geneListsAndIDs[1]
del geneListsAndIDs


count = 0
for accession in accessionList:
    annotation = allMeta[accession]
    count+=1

    #if ((('38batch' in annotation['folder']) or ('May11' in annotation['folder'])) and not(annotation['RNAseqFile'] == 'none')):
    if (not(annotation['RNAseqFile'] == 'none')):
        #RNAFile = annotation['RNAseqFile'][0]
        
        if len(annotation['RNAseqFile'][0]) > 1:
            RNAFile = annotation['RNAseqFile'][0]
        else:
            RNAFile = annotation['RNAseqFile']
        print(RNAFile)

        print(count)
        print(accession)
        print(RNAFile)

        annotationFolder = annotation['folder']
        print(annotationFolder)

        expAccession = RNAFile[-15:-4]
        expFile = annotationFolder +  'geneExp_dict_' + expAccession + '.pkl'
        annFile = annotation['bedFile']
        mnemFile = annotationFolder + 'mnemonics_v04.txt'
        #geneList
        extension = 3000
        SegwayTranscriptionEnrichment(annotationFolder, annFile, expFile, extension, geneList, geneIDList, mnemFile)

########################################
# 2. The ROC curves for the AUC values
########################################

