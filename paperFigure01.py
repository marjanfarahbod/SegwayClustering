# Plots for paper Figure 1
# 0. Initials
# 1. The heatmap for assay presence DONE
# 2. The actual annotation for samples sorted - I donno what this is
# 3. The long heatmap with all things sorted for all the signal. For genes without transcriptomic data, just the genomic regions(supplement). (so I want each "interpretation term" to be clustered, then them presented in a long heatmap)
# 4. The prediction probabilities DONE
# *. the big fat genome regional plot
# *. The prediction some sort of filtering?
# *. The other tracks bar thing plot

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

accessionList = list(allMeta.keys())
accessionList.remove(accessionList[205])
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

        color = bgcolors[int(meanBin*80)]

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
import linecache
# define the matrix for collecting
chrList = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
           'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']

W = 300
annotMat = np.zeros((234, W))
scor = 44900796 # chr19 44905796..44909393 APOE
scor = 234040702 # SPP2, Chr2, 30000
scor = 24280190 # chr7
scor = 72772659 # chr15
regionChr = 'chr19' # 'chr15'
regionChrInd = chrList.index(regionChr)
# for each of the accessions
for accession in accessionList:

    if accession == 'ENCSR424JDX':
        print(accessionList.index(accession))
        continue
    
    a = accessionList.index(accession)

    print(a, accession)
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']
    annFile = annotation['bedFile']

    if annFile.endswith('gz'):
        os.system('gunzip %s' %(annFile))
        annFile = annFile[0:-3]
    else:
        os.system('gunzip %s.gz' %(annFile))

    annLineCount = int(sp.getoutput('wc -l %s' %(annFile)).split()[0])
    
    mnemFile = annotationFolder + 'mnemonics_v04.txt'
    
    label_term_mapping = {}
    with open(mnemFile, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term

    slineInd = 1
    elineInd = annLineCount
    lineInd = int(elineInd/2)
    line = linecache.getline(annFile, lineInd)
    chr = line.split()[0]
    while chrList.index(chr) != regionChrInd:
    
        if chrList.index(chr) < regionChrInd:
            slineInd = lineInd
            lineInd = int(slineInd + ((elineInd-slineInd)/2))
            line = linecache.getline(annFile, lineInd)
            chr = line.split()[0]
            print(chr)
        else:
            if chrList.index(chr) > regionChrInd:
                elineInd = lineInd
                lineInd = int(slineInd + ((elineInd-slineInd)/2))
                line = linecache.getline(annFile, lineInd)
                chr = line.split()[0]

    sind = int(line.split()[1])
    if sind > scor: # staring coordinate
        eind = int(line.split()[2])
        while (sind > scor) or (eind < scor):
            lineInd = lineInd - 1
            line = linecache.getline(annFile, lineInd)
            sind = int(line.split()[1])
            eind = int(line.split()[2])
    else:
        if sind < scor:
            eind = int(line.split()[2])
            while (sind > scor) or (eind < scor):
                lineInd = lineInd + 1
                line = linecache.getline(annFile, lineInd)
                sind = int(line.split()[1])
                eind = int(line.split()[2])
            

    walker = 0
    steps = int((eind - scor)/100)
    color = label_term_mapping[line.split()[3].split('_')[0]]
    colorInd = segwayStates.index(color)
    annotMat[a, walker:walker+steps-1] = colorInd
    walker += steps
    print(steps)
    while walker < W:
        print(line)
        print(steps)
        lineInd +=1
        print(walker)
        print(color)
        line = linecache.getline(annFile, lineInd)
        sind = int(line.split()[1])
        eind = int(line.split()[2])
        steps = min(int((eind-sind)/100), W-walker)
        color = label_term_mapping[line.split()[3].split('_')[0]]
        colorInd = segwayStates.index(color)
        print(colorInd)
        annotMat[a, walker:walker+steps] = colorInd
        walker += steps
        
    linecache.clearcache()
    os.system('gzip %s' %(annFile))

    
# book = np.delete(annotMat, 205, 0) # not needed if filtered at accessionList
annotMat = annotMat[0:234, :]
print(annotMat.shape)
outputFileName = 'annotationMap234Cells_chr19_APOE.pkl'
outputFile = dataFolder + outputFileName
with open(outputFile, 'wb') as f:
    pickle.dump(annotMat, f)

#inputFileName = 'annotationMap234Cells_chr15.pkl'
inputFileName = 'annotationMap234Cells_chr7.pkl'
inputFile = dataFolder + inputFileName
with open(inputFile, 'rb') as f:
    annotMat = pickle.load(f)


sns.heatmap(annotMat)
plt.show()
# in the .bed file, reach out to the region
#scor = 73065000
scor = 72772659
chr15, 73065000 + 20k
 n
# colors
label_color_mapping = {'Promoter': [255,0,0],
 'PromoterFlanking': [255,68,0],
 'Enhancer': [255,195,77],
 'EnhancerLow': [255,255,0],
 'Bivalent': [189,183,107],
 'CTCF': [196,225,5],
 'Transcribed': [0,128,0],
 'K9K36': [102,205,170],
 'FacultativeHet': [128,0,128],
 'ConstitutiveHet': [138,145,208],
 'Quiescent': [255,255,255]}

mycolors = [[255,0,0],
          [255,68,0],
          [255,195,77],
          [255,255,0],
          [189,183,107],
          [196,225,5],
          [0,128,0],
          [102,205,170],
          [128,0,128],
          [138,145,208],
          [255,255,255]]

colorList = []
for color in mycolors:
    colorList.append(np.array(color)/256)


import matplotlib.colors as colors
cmap = colors.ListedColormap(colorList)
boundaries = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
norm = colors.BoundaryNorm(boundaries, cmap.N, clip=True)

# we have to do a change in sib.
sibplus = np.copy(sib)

for i,label in reversed(list(enumerate(label_order))):
    book = np.where(sib == label)
    sibplus[book] = i


df = pd.DataFrame(annotMat)

sib = sns.clustermap(df, method = 'ward', col_cluster=False)
rowInds = sib.dendrogram_row.reordered_ind

sns.heatmap(annotMat[rowInds,:], cmap=cmap, norm=norm)
plt.show()
print(figFile)
figFile = plotFolder + figureFolder + 'annotMat_chr7.pdf'
plt.savefig(figFile)
plt.close('all')
plt.show()

# removing the 205 from rowInds
#accessionList.pop(205)
for k, i in enumerate(rowInds):
    accession = accessionList[i]
    annotation = allMeta[accession]
    print(k,annotation['tissueInfo'][0])


########################################
# *. the sample feature plot - average 
########################################
feature_names_plotOrder_renamed = ['initial exon',  "5' flanking (1-1000 bp)",'initial intron', 'internal exons','internal introns', 'terminal exon',  'terminal intron', "3' flanking (1-1000 bp)" , "5' flanking (1000-10000 bp)","3' flanking (1000-10000 bp)",'H3K4me3','H3K27ac','H3K4me1', 'H3K36me3','H3K27me3','H3K9me3']

inputFile = dataFolder + 'allsamples_featureMat.pkl'
with open(inputFile, 'rb') as f:
    allFeatureMetaData = pickle.load(f)

meanMat = np.zeros((11, 16))
wholeMat = allFeatureMetaData['mat']
from scipy.stats import zscore
plot_data_z = zscore(wholeMat, axis = 0)
# plot_data_z_thr = np.where(plot_data_z > 1, 1.1, plot_data_z)

stateIndex = np.array(allFeatureMatData['termIndex'])
for i in range(11):
    book = np.where(stateIndex == i)[0]
    subMat = plot_data_z_thr[book, :]
    meanMat[i,:] = np.mean(subMat, axis=0)

meanMat_thr = np.where(meanMat >2, 2, meanMat)
cmap = plt.cm.coolwarm
norm = colors.BoundaryNorm(np.arange(-1, 2, .5), cmap.N)
plotDf = pd.DataFrame(meanMat_thr, segwayStates, feature_names_plotOrder_renamed)
sns.heatmap(plotDf, center = 0,cmap=cmap, norm=norm, linewidths=.1, linecolor='white')
#plt.show()
#sns.heatmap(meanMat, center = 0,cmap=cmap, norm=norm)
plt.title('average z-score of the feature values among the samples')
plt.tight_layout()

# plt.show()

figFile = plotFolder + figureFolder + 'averageSignal_heatmap_6sec.pdf'
print(figFile)
plt.savefig(figFile)
plt.close('all')

#################################################
# *. the sample feature plot - whole - supplement
################################################

# sort the segway state line and get the index of sort

feature_names_plotOrder_renamed = ['initial exon',  "5' flanking (1-1000 bp)",'initial intron', 'internal exons','internal introns', 'terminal exon',  'terminal intron', "3' flanking (1-1000 bp)" , "5' flanking (1000-10000 bp)","3' flanking (1000-10000 bp)",'H3K4me3','H3K27ac','H3K4me1', 'H3K36me3','H3K27me3','H3K9me3']

inputFile = dataFolder + 'allsamples_featureMat.pkl'
with open(inputFile, 'rb') as f:
    allFeatureMetaData = pickle.load(f)

termIndex = np.asarray(allFeatureMetaData['termIndex'])
termInd_sortInd = np.argsort(termIndex)

# get the colors
mycolors = [[255,0,0],
          [255,68,0],
          [255,195,77],
          [255,255,0],
          [189,183,107],
          [196,225,5],
          [0,128,0],
          [102,205,170],
          [128,0,128],
          [138,145,208],
          [255,255,255]]

termInd_sort = np.sort(termIndex)
plot_colors = []
for ind in termInd_sort:
    color = np.array(mycolors[ind])/255
    plot_colors.append(color)


meanMat = np.zeros((11, 16))
wholeMat = allFeatureMetaData['mat']
from scipy.stats import zscore
plot_data_z = zscore(wholeMat, axis = 0)

wholeMat_sorted = plot_data_z[termInd_sortInd, :]
wholeMat_zsthr = np.where(wholeMat_sorted > 3, 3, wholeMat_sorted)

wholeMat_2PosThr = np.where(wholeMat_sorted > 1.75, 1.75, wholeMat_sorted)
wholeMat_2PosNegThr = np.where(wholeMat_2PosThr < -1.75, -1.75, wholeMat_2PosThr)
plotDf = pd.DataFrame(wholeMat_2PosNegThr,  columns=feature_names_plotOrder_renamed)

plotDf = pd.DataFrame(wholeMat_zsthr,  columns=feature_names_plotOrder_renamed)

#sib = sns.clustermap(wholeMat_zsthr, figsize=(6, 12), cmap=cmap,row_colors= plot_colors, row_cluster=False, col_cluster=False)

sib = sns.clustermap(plotDf, figsize=(6, 8), cmap=cmap,row_colors= plot_colors, row_cluster=False, col_cluster=False)
plt.show()

plt.title('classifier input data sorted by predicted label - zscore for plot, threshold < 3')
plt.tight_layout()

# plt.show()
figFile = plotFolder + figureFolder + 'classifier_signal_supplement_175Filter.pdf'
print(figFile)
plt.savefig(figFile)
plt.close('all')


#################################################
# *. feature plot for one and multiple samples - for the diagram
################################################

inputFile = dataFolder + 'plot15samples_termSortedFeatureMat.pkl'
with open(inputFile, 'rb') as f:
    allFeatureMatData = pickle.load(f)


wholeMat = allFeatureMatData['mat']
from scipy.stats import zscore
plot_data_z = zscore(wholeMat, axis = 0)
plot_data_z_thr = np.where(plot_data_z > 3, 3, plot_data_z)

cmap = plt.cm.coolwarm
sns.heatmap(plot_data_z_thr[0:108,], center = 0,cmap=cmap, linewidths=.1, linecolor='white')
plt.show()
#sns.heatmap(meanMat, center = 0,cmap=cmap, norm=norm)
plt.title('feature for a few samples - zscore')
plt.tight_layout()

# plt.show()

figFile = plotFolder + figureFolder + 'featureHeatmap_15samplesSection.pdf'
print(figFile)
plt.savefig(figFile)
plt.close('all')

########################################
# * get the one line color from the heatmap
########################################

inputFileName = 'annotationMap234Cells_chr15.pkl'
inputFile = dataFolder + inputFileName
with open(inputFile, 'rb') as f:
    annotMat = pickle.load(f)


for accession in accessionList:
    annotation = allMeta[accession]
    if annotation['tissueInfo'][0] == 'K562':
        print(accession)
    
ind = accessionList.index('ENCSR019DPG')
book = annotMat[ind, :].reshape((200,1))

sns.heatmap(book, cmap=cmap, norm=norm)

figFile = plotFolder + figureFolder + 'oneSample_annotColor_K562.pdf'
print(figFile)
plt.savefig(figFile)
plt.close('all')

plt.show()

########################################
# * length distribution plots
########################################


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


labelCoverageLists = {} # put the labels with the same state individually
stateCoverageLists = {} # sum the labels with the same state
for s in segwayStates:
    print(s)
    labelCoverageLists[s] = []
    stateCoverageLists[s] = []

for accession in accessionList:
    annotation = allMeta[accession]
    annFolder = annotation['folder']

    mnemFile = annFolder + 'mnemonics_v04.txt'
    
    label_term_mapping = {}
    with open(mnemFile, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term

    sampleStateCoverage = np.zeros(11)
    lfile = length_files[accession]     
    with open(lfile, 'r') as f:
        for i, line in enumerate(f):
            if i < 2:
                continue
            #print(line)
            label = line.split()[0]
            coverage = float(line.split()[6])
            term = label_term_mapping[label]
            #print(term)
            labelCoverageLists[term].append(coverage)
            termInd = segwayStates.index(term)
            #print(termInd)
            #print(coverage)
            sampleStateCoverage[termInd]+= coverage

    for i,state in enumerate(segwayStates):
        stateCoverageLists[state].append(sampleStateCoverage[i])


plotData = []
for i in range(11):
    print(i)
    plotData.append(stateCoverageLists[segwayStates[i]])

plt.boxplot(plotData)
plt.grid()
figFile = plotFolder + figureFolder + 'genomeCoverage.pdf'
print(figFile)
plt.savefig(figFile)
plt.close('all')

plt.show()

########################################
# *. The other tracks bar thing plot
########################################
from scipy.stats import zscore

runID_file = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/testBatch105/' + 'runID_accession_map_105run.pkl'
with open(runID_file, 'rb') as pickledFile:
    runID_accession105 = pickle.load(pickledFile)

accession_runID = {}
for runID in list(runID_accession105.keys()):
    ac = runID_accession105[runID]
    accession_runID[ac] = runID

signal_files = {}
for accession in accessionList:
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']
    #if accession == 'ENCSR121RCG':
    #    continue

    if '38batch' in annotationFolder:
        signal_files[accession] = annotationFolder + 'segOutput/call_signal_distribution/signal_distribution.tab'

    if 'May11' in annotationFolder:
        signal_files[accession] = annotationFolder + 'call-segtools/signal_distribution.tab'

    if 'Batch105' in annotationFolder:
        runID = accession_runID[accession]
        signal_files[accession] = dataFolder + 'testBatch105/all_segtools/' + runID +  '/signal_distribution/signal_distribution.tab'

#smallTrackList =  ['CTCF', 'DNase-seq', 'ATAC-seq']
smallTrackList = trackList[6:]
# we do the zscore along the column FIRST, and then do the bar for each of the labels

trackLabelKeyList = []
for i in range(len(segwayStates)):
    for j in range(len(smallTrackList)):
        trackLabelKeyList.append((i,j))

values_dict = {}
for key in trackLabelKeyList:
    values_dict[key] = []

# for each of the tracks, we will have a list of values and segway labels.
stlValue_dict = {}
for track in smallTrackList:
    info = {}
    info['values'] = []
    info['labels'] = []
    stlValue_dict[track] = info

for accession in accessionList:

    print(accessionList.index(accession))
    print(accession)
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

    labelCount = len(label_term_mapping)

    #  get the mapping
    # >>>>>>>>>>>>>>>>>>>>>>>>>>
    mapping_file = annotationFolder + 'trackname_assay.txt'
    # read the mapping_file
    track_assay_map = {}
    sampleTrack_list = [] # list of the track type for the sample 
    with open(mapping_file) as inputFile:
        for line in inputFile:
            fields = line.strip().split()
            track_assay_map[fields[0]] = fields[1]
            sampleTrack_list.append(fields[1])
    sampleTrack_count = len(sampleTrack_list)

    sampleTrack_list_sorted = []
    for track in trackList:
        for sampleTrack in sampleTrack_list:
            if sampleTrack == track:
                sampleTrack_list_sorted.append(sampleTrack)

    # ge the signal dist file
     # >>>>>>>>>>>>>>>>>>>>>>>>>>
    signal_file = signal_files[accession]

    signal_dist_mat = np.zeros((labelCount, sampleTrack_count))
    with open(signal_file, 'r') as inputFile:
        header = inputFile.readline()
        for line in inputFile:
            fields = line.strip().split()
            #print(line)
            track = track_assay_map[fields[1]]
            track_ind = sampleTrack_list_sorted.index(track) # track list order
            cluster_ind = int(fields[0]) # cluster list order
            signal_dist_mat[cluster_ind][track_ind] = round(float(fields[2]), 4)
            
    z_signal_dist_mat = zscore(signal_dist_mat, axis=0)

    for t, track in enumerate(smallTrackList):
        if track in sampleTrack_list:
            trackInd = sampleTrack_list_sorted.index(track)
            for label in range(labelCount):
                value = z_signal_dist_mat[label, trackInd]
                segwayStateInd = segwayStates.index(label_term_mapping[str(label)])
                stlValue_dict[track]['values'].append(value)
                stlValue_dict[track]['labels'].append(segwayStateInd)
                values_dict[(segwayStateInd, t)].append(value)

# find the max and the mean
for track in smallTrackList:
    myVals = np.array(stlValue_dict[track]['values'])
    print(np.min(myVals))
    print(np.max(myVals))

# the range for the x axis is -4 to 4

# figure 11 * 3

# for each plot, fetch the track, and then labels. 
from matplotlib import colors
fig, axs = plt.subplots( len(segwayStates), len(smallTrackList), figsize=(4,8))
 
#cmap = sns.diverging_palette(245, 15, s=80, l=40, sep=2, n = 15, as_cmap=True, center='dark')
#cmap = sns.diverging_palette(245, 15, s=80, l=40, sep=2, n = 15, as_cmap=True, center='dark')
cmap = plt.cm.get_cmap('Oranges')
cmap = plt.cm.get_cmap('coolwarm')

bgcolors = []
for i in range(100):
    bgcolors.append(cmap(i*5))
    
bins = np.linspace(-3, 3, 21)
binCount = 20
barx = np.zeros(binCount)
for k in range(binCount):
    barx[k]= np.mean([bins[k], bins[k+1]])

meanBinMat = np.zeros((segwayStateCount, len(smallTrackList)))
for i in range(segwayStateCount):
    for j in range(len(smallTrackList)):
        print(i, j)
        #n, bins, patches = plt.hist(probDistDict[(i,j)], bins=binCount, density=True)
        #bins = np.linspace(0, 1, 16)
        hist = np.histogram(values_dict[(i,j)], bins, density=True)

        #meanBin = np.mean(hist[0])
        sib = [hist[0][i]*hist[1][i] for i in range(len(hist[0]))]
        meanBin = np.mean(sib)
        meanBinMat[i,j] = meanBin
        print(meanBin)
        #if meanBin > .18:
        #    meanBin = .18

        #print('....', int(meanBin + np.abs(meanBin*300)))
        #color = bgcolors[int((meanBin + .179)*140)]
        
        axs[i,j].set_facecolor(color)
        axs[i,j].bar(barx, hist[0], width=.5, color='black')
        #axs[i,j].bar(barx, n, width=.06, color='black')
        #axs[i,j].hist(probDistDict[(i,j)], bins=15, density=True, color=['blue', 'green', 'black'])
        axs[i,j].set_xlim(-3,3)
        axs[i,j].set_ylim(0,2)
        #axs[i,j].set_xticks([0, 1])
        #axs[i,j].set_yticks([0, 4, 8])
        axs[i,j].set_xticks([])
        axs[i,j].set_yticks([])
        axs[i,j].set_yticklabels([])
        axs[i,j].set_xticklabels([])
        myMean = np.mean(values_dict[(i,j)])
        axs[i,j].plot([myMean, myMean], [0, 2], 'w', lineWidth ='1.1')

#plt.show()
plotFolder_add = plotFolder + figureFolder
figFile = plotFolder_add + 'sixTracks_dist_heatmap_meanLine.pdf'
print(figFile)
plt.tight_layout()
plt.savefig(figFile)
plt.close('all')

# just the redo of the above plot
########################################
# for each plot, fetch the track, and then labels. 
from matplotlib import colors
fig, axs = plt.subplots( len(segwayStates), len(smallTrackList), figsize=(4,8))
 
#cmap = sns.diverging_palette(245, 15, s=80, l=40, sep=2, n = 15, as_cmap=True, center='dark')
#cmap = sns.diverging_palette(245, 15, s=80, l=40, sep=2, n = 15, as_cmap=True, center='dark')
cmap = plt.cm.get_cmap('Oranges')
cmap = plt.cm.get_cmap('coolwarm')

bgcolors = []
for i in range(170):
    bgcolors.append(cmap(i*5))
    
bins = np.linspace(-3.5, 3.5, 21)
binCount = 20
barx = np.zeros(binCount)
for k in range(binCount):
    barx[k]= np.mean([bins[k], bins[k+1]])

myMeanMat = np.zeros((segwayStateCount, len(smallTrackList)))
for i in range(segwayStateCount):
    for j in range(len(smallTrackList)):
        print(i, j)
        #n, bins, patches = plt.hist(probDistDict[(i,j)], bins=binCount, density=True)
        #bins = np.linspace(0, 1, 16)
        hist = np.histogram(values_dict[(i,j)], bins, density=True)
        myMeanMat[i,j] = np.mean(values_dict[i,j])

minMyMean = np.min(myMeanMat)
maxMyMean = np.max(myMeanMat)
dis = maxMyMean - minMyMean
fact = 1/dis

adjustedMean = np.zeros((segwayStateCount, len(smallTrackList)))
for i in range(segwayStateCount):
    for j in range(len(smallTrackList)):
        adjustedMean[i,j] = (myMeanMat[i,j] - minMyMean) * fact

meanBinMat = np.zeros((segwayStateCount, len(smallTrackList)))
for i in range(segwayStateCount):
    for j in range(len(smallTrackList)):
        print(i, j)
        #n, bins, patches = plt.hist(probDistDict[(i,j)], bins=binCount, density=True)
        #bins = np.linspace(0, 1, 16)
        hist = np.histogram(values_dict[(i,j)], bins, density=True)
        #myMean = np.mean(values_dict[i,j])

        #meanBinMat[i,j] = myMean
        #print(meanBin)
        #if meanBin > .18:
        #    meanBin = .18

        #print('....', int(meanBin + np.abs(meanBin*300)))
        color = bgcolors[int(adjustedMean[i,j]*70)]

        axs[i,j].set_facecolor(color)
        axs[i,j].bar(barx, hist[0], width=.5, color='black')
        #axs[i,j].bar(barx, n, width=.06, color='black')
        #axs[i,j].hist(probDistDict[(i,j)], bins=15, density=True, color=['blue', 'green', 'black'])
        axs[i,j].set_xlim(-3.5,3.5)
        axs[i,j].set_ylim(0,1.5)
        #axs[i,j].set_xticks([0, 1])
        #axs[i,j].set_yticks([0, 4, 8])
        axs[i,j].set_xticks([])
        axs[i,j].set_yticks([])
        axs[i,j].set_yticklabels([])
        axs[i,j].set_xticklabels([])
        myMean = np.mean(values_dict[(i,j)])
        axs[i,j].plot([0, 0], [0, 2], color = [11/256, 249/256, 21/256, 1], lineWidth ='1.3')
        axs[i,j].plot([0, 0], [0, 2], 'w', lineWidth ='1.3', linestyle='dashed')
        #axs[i,j].plot([0, 0], [0, 2], 'w', lineWidth ='1.1')
        
#plt.show()
plotFolder_add = plotFolder + figureFolder
figFile = plotFolder_add + 'sixTracks_dist_heatmap_meanLineV06.pdf'
print(figFile)
plt.tight_layout()
plt.savefig(figFile)
plt.close('all')

                
