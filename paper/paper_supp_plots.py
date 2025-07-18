# Get some paper data and supplementary figures
# 1. count of track data
# 2. label-length distribution plot
# 3. count of samples with different activites
# 4. coverage of the labels between the two annotations, 4.1 ChromHMM 4.2 Segway
# 5. ROC curve for the samples

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

book =  'EnhA1 EnhA2 EnhBiv EnhG1 EnhG2 EnhWk Het Quies ReprPC ReprPCWk TssA TssBiv TssFlnk TssFlnkD TssFlnkU Tx TxWk ZNF/Rpts'
chromLabels = book.split()
print(len(chromLabels))
chromLabels_reordered = ['TssA', 'TssFlnk', 'TssFlnkD', 'TssFlnkU', 'EnhA1', 'EnhA2', 'EnhG1', 'EnhG2',
                         'EnhWk', 'TssBiv', 'EnhBiv', 'Tx', 'TxWk', 'ZNF/Rpts', 'ReprPC', 'ReprPCWk' , 'Het', 'Quies']

inputFileName = 'all235Annot_meta_corrected.pkl'
#inputFileName = 'all235Annot_meta.pkl'
inputFile = dataFolder +  inputFileName
with open(inputFile, 'rb') as f:
    allMeta = pickle.load(f)

accessionList = list(allMeta.keys())
accessionList.remove(accessionList[205])

########################################
# 1. count of track data
########################################
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

########################################
# 2. label-length distribution plot
########################################

# get the mean and median of the length distribution for each of our labels.
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

    if '38batch' in annotationFolder:
        length_files[accession] = annotationFolder + 'segOutput/length_distribution/segment_sizes.tab'


    if 'May11' in annotationFolder:
        length_files[accession] = annotationFolder + 'call-segtools/segment_sizes.tab'

    if 'Batch105' in annotationFolder:
        runID = accession_runID[accession]
        length_files[accession] = dataFolder + 'testBatch105/all_segtools/' + runID +  '/length_distribution/segment_sizes.tab'


meanLengthLists = np.zeros((1000, 11)) # count of labels for each sample, count of segway states
medLengthLists = np.zeros((1000, 11)) # count of labels for each sample, count of segway states
counts = np.zeros(11) # count of the labels in each group
for accession in accessionList:

    annotation = allMeta[accession]
    annotationFolder = annotation['folder']

    # get the length_dist file
    lengthFile = length_files[accession]

    # the mnemonics file
    label_term_mapping = {}
    mnemonics_file = annotationFolder + 'mnemonics_v04.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term

    labelCount = len(label_term_mapping)

    with open(lengthFile, 'r') as f:
        line = f.readline()
        line = f.readline()
        for line in f:
            fields = line.strip().split()
            label = fields[0]
            meanL = np.float(fields[2])
            medL = np.float(fields[3])
            termInd = segwayStates.index(label_term_mapping[fields[0]])
            medLengthLists[int(counts[termInd]), termInd] = medL
            meanLengthLists[int(counts[termInd]), termInd] = meanL
            counts[termInd] += 1
        

stateMeans = {}
stateMeds = {}
for i, state in enumerate(segwayStates):
    stateMeans[state] = meanLengthLists[0:int(counts[i]), i]
    stateMeds[state] = medLengthLists[0:int(counts[i]), i]

for i, state in enumerate(segwayStates):
    plt.boxplot(stateMeans.values())

# note: ended up not using this plot

plt.show()

########################################
# 3. count of samples with different activities
########################################

groupStates = ['Promoter', 'Enhancer', 'Bivalent', 'CTCF', 'Transcribed', 'K9K36', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent']
labelCounts = np.zeros(9)
for accession in accessionList:

    annotation = allMeta[accession]
    annotationFolder = annotation['folder']

    # the mnemonics file
    label_term_mapping = {}
    mnemonics_file = annotationFolder + 'mnemonics_v04.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term
            
    tempCounts = np.zeros(9)
    for label in label_term_mapping:
        term = label_term_mapping[label]

        for i, state in enumerate(groupStates):
            if term in state:
                tempCounts[i] = 1

    if tempCounts[0] == 0 or tempCounts[1] ==0 or tempCounts[4] == 0:
        print(accession)
        print(annotationFolder)
        if accession in pSamples:
            print('---------')
        
    labelCounts += tempCounts

plt.grid(axis = 'y')
plt.bar(range(9), (labelCounts/234)*100, color=[0,0,0,1], width=.75)
plt.xticks(range(9))
plt.ylim([0, 100])
plt.show()

########################################
# 4. coverage of the labels between the two annotations
########################################

# 4.1 ChromHMM
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

accessionList = list(allMeta.keys())
clusterCount = len(chromLabels)
sampleCount = len(accessionList)

labelCoverage = np.zeros((sampleCount, 18))
for a, accession in enumerate(accessionList):
    annotation = allMeta[accession]
    print(a)

    sampleGenomeCoverage = 0

    #if ((('38batch' in annotation['folder']) or ('May11' in annotation['folder'])) and not(annotation['RNAseqFile'] == 'none')):
    if (annotation['chromFile'] == 'none'):
        print('no chrom file', accession)
        continue


    annotationFolder = annotation['folder']
    print(annotationFolder)

    chmmFileName = annotation['chromFile'].split('/')[-1]

    if (('38batch' in annotation['folder']) or ('May11' in annotation['folder'])):
        annFile = annotationFolder + 'sorted_' + chmmFileName
    else:
        annFile = annotationFolder + chmmFileName

    # prepare ann
    if annFile.endswith('.gz'):
        os.system('gunzip %s' %(annFile))
        annFile = annFile[0:-3]
    else:
        os.system('gunzip %s.gz' %(annFile))

    annLineCount = int(sp.getoutput('wc -l %s' %(annFile)).split()[0])

    annLineInd = 1 # this keeps the annotation line index, the first line is empty, thus the index starts at 1
    line = linecache.getline(annFile, annLineInd)
    with open(annFile, 'r') as f:
        f.readline()

        for line in f:
            fields = line.strip().split()

            ann_start = int(fields[1])
            ann_end = int(fields[2])
            annLabel = fields[3]
            sampleLabelCoverage = ann_end - ann_start
    
            sampleGenomeCoverage += sampleLabelCoverage
            
            labelInd = chromLabels.index(annLabel)
    
            labelCoverage[a, labelInd] += sampleLabelCoverage


    os.system('gzip %s' %(annFile))


# coverage of each label for ChromHMM + genome length covered

output = {}
output['labelCoverage'] = labelCoverage
output['accessionList'] = accessionList
output['chromLabels'] = chromLabels
outputFile = dataFolder + 'chromLabelCoverage.pkl'
with open(outputFile, 'wb') as f:
    pickle.dump(output, f)

inputFile = dataFolder + 'chromLabelCoverage.pkl'
with open(inputFile, 'rb') as f:
    chromLabelCoverage = pickle.load(f)

# coverage of each label for chrom
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# for each chrom label
labelCoverage = chromLabelCoverage['labelCoverage']
nonZero = justSums > 0

labelCoverage_per = labelCoverage[nonZero,:] / np.sum(labelCoverage[nonZero,:], axis=1)[:, np.newaxis] # chrom

book = pd.DataFrame(labelCoverage_per, columns = chromLabels)
book = book[chromLabels_reordered]
#book.mask(book==0)
boxplot = book.boxplot(rot = 90)
book.mask(book == 0).boxplot(rot=90)
#plt.grid(False)
#plt.boxplot(labelCoverageDel)
plt.ylim([-0.04,1])
plt.yticks([0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1])
#plt.grid(False)
plt.grid(axis='x')
plt.title('coverage of labels - chrom')
plt.title('coverage of labels - zeros removed - chrom')
plt.tight_layout()
#plt.show()

figFile = plotFolder + 'chromLabelCoverage_zerosRemoved.pdf'
figFile = plotFolder + 'chromLabelCoverage_allValues.pdf'
print(figFile)
plt.savefig(figFile)
plt.close('all')

# coverage of Segway groupings for Segway
# we don't have this
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
labelCoverageDel = labelCoverage[nonZero, :]
df = pd.DataFrame(labelCoverageDel, columns=chromLabels)
df = df[chromLabels_reordered]
labelCoverage_reordered = df.to_numpy()
activeCover = np.sum(labelCoverage_reordered[:, 0:14], axis=1)
activeCover_per = activeCover / np.sum(labelCoverage_reordered, axis=1)
silencedCover = np.sum(labelCoverage_reordered[:, 14:], axis=1)
silencedCover_per = silencedCover / np.sum(labelCoverage_reordered, axis=1)

columns = {'active':activeCover_per, 'silent':silencedCover_per}
df = pd.DataFrame(data=columns)

df.boxplot()
plt.ylim([-0.04,1])
plt.yticks([0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1])

# plt.grid(False)
plt.grid(axis = 'x')
plt.title('Chrom silent and active region coverage')
plt.tight_layout()
figFile = plotFolder + 'Chrom_active_silent_coverage.pdf'
print(figFile)
plt.savefig(figFile)
plt.close('all')
plt.show()

# coverage of Segway groupings for ChromHMM
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
df = pd.DataFrame(labelCoverage, columns=chromLabels)
df = df[chromLabels_reordered]
labelCoverage_reordered = df.to_numpy()

# promoter 0
pf = np.sum(labelCoverage_reordered[:, 1:4], axis=1) # promoterFlank 1:4
en = np.sum(labelCoverage_reordered[:, 4:8], axis=1) # enhancer 4:9
# enhancerWeak 8
bv = np.sum(labelCoverage_reordered[:, 9:11], axis=1) # bivalent 9:11
tr = np.sum(labelCoverage_reordered[:, 11:13], axis=1) # transcribed
# k9k36 13
fh = np.sum(labelCoverage_reordered[:, 14:16], axis=1) # facultativeHet
# constitutiveHet 16
# quiescent 17

chromDataSegwayLabels = np.matrix([labelCoverage_reordered[:, 0],
                                   pf,
                                   en,
                                   labelCoverage_reordered[:, 8],
                                   bv,
                                   tr,
                                   labelCoverage_reordered[:, 13],
                                   fh,
                                   labelCoverage_reordered[:, 16],
                                   labelCoverage_reordered[:,17]])

chromDataSegwayLabels = chromDataSegwayLabels.transpose()
segwayStatesNoCT = segwayStates.copy()
segwayStatesNoCT.remove('CTCF')
justSums = np.sum(labelCoverage, axis=1)
nonZero = justSums > 0

labelCoverage_per =  chromDataSegwayLabels[nonZero,:]/ np.sum(labelCoverage[nonZero,:], axis=1)[:, np.newaxis]
chromDataSegwayLabels_df = pd.DataFrame(labelCoverage_per, columns=segwayStatesNoCT )

chromDataSegwayLabels_df.boxplot(rot=90)
#plt.grid(False)
plt.grid(axis='x')
plt.ylim([-0.04,1])
plt.yticks([0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1])
plt.title('Chrom, segway label coverage')
plt.tight_layout()
figFile = plotFolder + 'Chrom_segwayLabelCoverage.pdf'
print(figFile)
plt.savefig(figFile)
plt.close('all')
plt.show()

# coverage of active and silent regions for ChromHMM

# 4.2 Segway
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

accessionList = list(allMeta.keys())
clusterCount = len(segwayStates)
sampleCount = len(accessionList)

genomeCoverage = np.zeros(sampleCount)
labelCoverage = np.zeros((sampleCount, 11))
for a, accession in enumerate(accessionList[206:]):
    a = a + 206
    annotation = allMeta[accession]
    print(a)

    sampleGenomeCoverage = 0

    #if ((('38batch' in annotation['folder']) or ('May11' in annotation['folder'])) and not(annotation['RNAseqFile'] == 'none')):
    #if (annotation['chromFile'] == 'none'):
    #    print('no chrom file', accession)
    #    continue

    annotationFolder = annotation['folder']
    print(annotationFolder)

    
    # get the mnemonics
    mnemFile = annotationFolder + 'mnemonics_v04.txt'
    label_term_mapping = {}
    with open(mnemFile, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term


    annFile = annotation['bedFile']
    
    # prepare ann
    if annFile.endswith('.gz'):
        os.system('gunzip %s' %(annFile))
        annFile = annFile[0:-3]
    else:
        os.system('gunzip %s.gz' %(annFile))

    annLineCount = int(sp.getoutput('wc -l %s' %(annFile)).split()[0])

    annLineInd = 1 # this keeps the annotation line index, the first line is empty, thus the index starts at 1
    line = linecache.getline(annFile, annLineInd)

    while(annLineInd < annLineCount -1):

        line = linecache.getline(annFile, annLineInd)
        annLineInd +=1
        
        fields = line.strip().split()

        ann_start = int(fields[1])
        ann_end = int(fields[2])
        annLabel = fields[3].split('_')[0]
        sampleLabelCoverage = ann_end - ann_start
    
        sampleGenomeCoverage += sampleLabelCoverage
            
        labelInd = segwayStates.index(label_term_mapping[annLabel])
    
        labelCoverage[a, labelInd] += sampleLabelCoverage

    linecache.clearcache()

    os.system('gzip %s' %(annFile))

labelCoverageDel = np.delete(labelCoverage, 205, 0)
accessionList.remove(accessionList[205])

output = {}
output['labelCoverage'] = labelCoverageDel
output['accessionList'] = accessionList
output['chromLabels'] = segwayStates
outputFile = dataFolder + 'segwayLabelCoverage.pkl'
with open(outputFile, 'wb') as f:
    pickle.dump(output, f)

inputFile = dataFolder + 'segwayLabelCoverage.pkl'
with open(inputFile, 'rb') as f:
    segwayLabelCoverage = pickle.load(f)

labelCoverageDel = segwayLabelCoverage['labelCoverage']

labelCoverageDel_per = labelCoverageDel / np.sum(labelCoverageDel, axis=1)[:, np.newaxis] # segway
    
# coverage of each label for Segway
book = pd.DataFrame(labelCoverageDel_per, columns = segwayStates)
book.mask(book==0)
book.replace(0.0,np.nan)
boxplot = book.boxplot(rot=90)
book.mask(book == 0).boxplot(rot=90)

#plt.boxplot(labelCoverageDel)
# plt.grid(False)
plt.grid(axis='x')
plt.ylim([-0.04,1])
plt.yticks([0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1])
plt.title('coverage of labels - zeros removed - segway')
plt.title('coverage of labels a- segway')
plt.tight_layout()
#plt.show()

figFile = plotFolder + 'segwayLabelCoverage_zerosRemoved.pdf'
figFile = plotFolder + 'segwayLabelCoverage_allValues.pdf'
print(figFile)
plt.savefig(figFile)
plt.close('all')

# coverage of Segway groupings for Segway
# we don't have this

activeCover = np.sum(labelCoverageDel[:, 0:8], axis=1)
activeCover_per = activeCover / np.sum(labelCoverageDel, axis=1)
silencedCover = np.sum(labelCoverageDel[:, 8:], axis=1)
silencedCover_per = silencedCover / np.sum(labelCoverageDel, axis=1)

columns = {'active':activeCover_per, 'silent':silencedCover_per}
df = pd.DataFrame(data=columns)

df.boxplot()
#plt.grid(False)
plt.grid(axis='x')
plt.ylim([-0.04,1])
plt.yticks([0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1])
plt.title('Segway silent and active region coverage')
figFile = plotFolder + 'Segway_active_silent_coverage.pdf'
print(figFile)
plt.savefig(figFile)
plt.close('all')
plt.show()


# 5. ROC curve for the samples
########################################

# get it from 5.1 


