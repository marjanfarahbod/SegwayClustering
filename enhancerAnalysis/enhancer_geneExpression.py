# so here we are.
# for one file, the coverage of enhancer label in 2k, 5k, 10k around genes 3 matrices
#
# 0. Initials
# 1. For a given sample get the AUCs based on the coverage of Enhancer and Enhancer low and together, for all samples together 
# 1.1 for each sample with transcriptomic data, run the function for both chrom and Segway
# 1.2 for each sample with transcriptomic data, open the file and get the AUC
# 2. For the set of example samples, get the bigwig file
# 2.1 get the accession name for the files
# 2.2 run the function for each of the files
# 2.3 get the output
# 3. The bar plots
# 3.1.0 Segway: get the chr20 section of the .bed file
# 3.1.1 Chrom: get the chr20 section of the .bed file
# 3.2.0 Segway colors
# 3.2.1 chrom colors
# 3.3.0 Segway, get the plot
# 3.3.1 chrom, get the plot
# 4. The three section plot
# 5. intensity of H3K4me1 within a distance from the genes
# 6. Between the two samples, get the enhancers between the genes that are expressed vs genes that are not

########################################
# 0. Initials
########################################
import linecache
import pickle
import re
import numpy as np
import pandas as pd
import subprocess as sp
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from QC_transcriptionComparison_util import Gene, Exon, Annotation, AnnotationClass
import pyBigWig

import glob
 
# General data folder
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/'
figureFolder = 'figure02/'


inputFileName = 'all235Annot_meta_corrected.pkl'
#inputFileName = 'all235Annot_meta.pkl'
inputFile = dataFolder +  inputFileName
with open(inputFile, 'rb') as f:
    allMeta = pickle.load(f)

# GTF data structure
fileName = dataFolder + '/geneLists.pkl'
with open(fileName, 'rb') as pickleFile:
    geneListsAndIDs = pickle.load(pickleFile)

book =  'EnhA1 EnhA2 EnhBiv EnhG1 EnhG2 EnhWk Het Quies ReprPC ReprPCWk TssA TssBiv TssFlnk TssFlnkD TssFlnkU Tx TxWk ZNF/Rpts'
chromLabels = book.split()

segwayStates = ['Promoter', 'PromoterFlanking', 'Enhancer', 'EnhancerLow', 'Bivalent', 'CTCF', 'Transcribed', 'K9K36', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent']
print(len(segwayStates))
segwayStateCount = len(segwayStates)

geneList = geneListsAndIDs[0]
geneIDList = geneListsAndIDs[1]
del geneListsAndIDs

accessionList = list(allMeta.keys())
accessionList.remove(accessionList[205])

########################################
# 1. For a given sample get the AUCs based on the coverage of Enhancer and Enhancer low and together, for all samples together 
########################################

# 1.1 for each sample with transcriptomic data, run the function for both chrom and Segway
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

for accession in accessionList:
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']

    if not(annotation['RNAseqFile'] == 'none') and not(annotation['chromFile'] == 'none'):

        # >>>>>>>>>> Segway, ev = 2000
        annFile = annotation['bedFile']
        mnemFile = annotationFolder + 'mnemonics_v04.txt'
        ev = 2000
        side = 'right'
        EPLabelCover_Segway_extensionValue(annotationFolder,
                                           annFile, mnemFile, geneList, geneIDList, ev, side)

        annFile = annotation['bedFile']
        mnemFile = annotationFolder + 'mnemonics_v04.txt'
        ev = 2000
        side = 'left'
        EPLabelCover_Segway_extensionValue(annotationFolder,
                                           annFile, mnemFile, geneList, geneIDList, ev, side)

        print(ev, accession)
        # >>>>>>>>>> Segway, ev = 5000
        annFile = annotation['bedFile']
        mnemFile = annotationFolder + 'mnemonics_v04.txt'
        ev = 5000
        side = 'right'
        EPLabelCover_Segway_extensionValue(annotationFolder,
                                           annFile, mnemFile, geneList, geneIDList, ev, side)

        annFile = annotation['bedFile']
        mnemFile = annotationFolder + 'mnemonics_v04.txt'
        ev = 5000
        side = 'left'
        EPLabelCover_Segway_extensionValue(annotationFolder,
                                           annFile, mnemFile, geneList, geneIDList, ev, side)

        # >>>>>>>>>> Segway, ev = 10000
        annFile = annotation['bedFile']
        mnemFile = annotationFolder + 'mnemonics_v04.txt'
        ev = 10000
        side = 'right'
        EPLabelCover_Segway_extensionValue(annotationFolder,
                                           annFile, mnemFile, geneList, geneIDList, ev, side)

        annFile = annotation['bedFile']
        mnemFile = annotationFolder + 'mnemonics_v04.txt'
        ev = 10000
        side = 'left'
        EPLabelCover_Segway_extensionValue(annotationFolder,
                                           annFile, mnemFile, geneList, geneIDList, ev, side)

        
        # >>>>>>>>>> Chrom, ev = 2000
        chmmFileName = annotation['chromFile'].split('/')[-1]

        if (('38batch' in annotation['folder']) or ('May11' in annotation['folder'])):
            annFile = annotationFolder + 'sorted_' + chmmFileName
        else:
            annFile = annotationFolder + chmmFileName
            
        ev = 2000
        side = 'right'
        EPLabelCover_chromHMM_extensionValue(annotationFolder, annFile, chromLabels, geneList, geneIDList, ev, side)
        # chromhmm has missing regions in the genome. Its coverage is ...
                    
        ev = 2000
        side = 'left'
        EPLabelCover_chromHMM_extensionValue(annotationFolder, annFile, chromLabels, geneList, geneIDList, ev, side)

        # >>>>>>>>>> Chrom, ev = 5000
                    
        ev = 5000
        side = 'right'
        EPLabelCover_chromHMM_extensionValue(annotationFolder, annFile, chromLabels, geneList, geneIDList, ev, side)
        # chromhmm has missing regions in the genome. Its coverage is ...
                    
        ev = 5000
        side = 'left'
        EPLabelCover_chromHMM_extensionValue(annotationFolder, annFile, chromLabels, geneList, geneIDList, ev, side)

        # >>>>>>>>>> Chrom, ev = 10000
                    
        ev = 10000
        side = 'right'
        EPLabelCover_chromHMM_extensionValue(annotationFolder, annFile, chromLabels, geneList, geneIDList, ev, side)
        # chromhmm has missing regions in the genome. Its coverage is ...
                    
        ev = 10000
        side = 'left'
        EPLabelCover_chromHMM_extensionValue(annotationFolder, annFile, chromLabels, geneList, geneIDList, ev, side)

# 1.2 for each sample with transcriptomic data, open the file and get the AUC
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# regionMode, expMode, regMode
switchLists = [['segway','chrom'], [2000, 5000, 10000]]

modeKeyList = []
for i in range(2):
    model = switchLists[0][i]
    for j in range(3):
        ev = switchLists[1][j]
        myKey = (model, ev)
        modeKeyList.append(myKey)

import random

rnaList = []
for a, accession in enumerate(accessionList):
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']
    if not(annotation['RNAseqFile'] == 'none') and not(annotation['chromFile'] == 'none'):
        rnaList.append(a)

self = True # true or false it can be
segwayAUCs = np.zeros((88, 3))
chromAUCs = np.zeros((88, 3))
count = 0
evList = [2000, 5000, 10000]
aucCurves = {}
rnaAccessionList = []
for a, accession in enumerate(accessionList):
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']
    aucOutput = {}
    if not(annotation['RNAseqFile'] == 'none') and not(annotation['chromFile'] == 'none'):
        print(a)
        rnaAccessionList.append(accession)
        # Segway
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        label_term_mapping = {}
        mnemonics_file = annotationFolder + 'mnemonics_v04.txt'
        with open(mnemonics_file, 'r') as mnemonics:
            for line in mnemonics:
                #print(line)
                label = line.strip().split()[0]
                term = line.strip().split()[1]
                label_term_mapping[label] = term

        # >>>>>>>>>>>>>>>>>>>>>>>>>>
        if self:
            rnaAccession = accession
            rnaAnnotation = annotation
        else:
            temp = random.randint(0,87)
            rnaAInd = rnaList[temp]
            if rnaAInd == a:
                rnaAInd = rnaList[temp+1]
            rnaAccession = accessionList[rnaAInd]
            rnaAnnotation = allMeta[rnaAccession]
            
        rnaFolder = rnaAnnotation['folder']
                
        if len(rnaAnnotation['RNAseqFile'][0]) > 1:
            RNAFile = rnaAnnotation['RNAseqFile'][0]
        else:
            RNAFile = rnaAnnotation['RNAseqFile']
            print(RNAFile)
            
        expAccession = RNAFile[-15:-4]
        expFile = rnaFolder +  'geneExp_dict_' + expAccession + '.pkl'
        print(expFile)
        print(annotationFolder)

        with open(expFile, 'rb') as pickledFile:
            expression = pickle.load(pickledFile) # load the expression file

        expArray = np.asarray([expression[x] for x in geneIDList])[0:24700] # genes in the transcriptomic data that are in the gene map file
        
        # >>>>>>>>>> Segway evList
        model = 'segway'
        for e, ev in enumerate(evList):
            side = 'right'
            inputFile = annotationFolder  + 'extendedGene_labelCover_%d_%s.pkl' %(ev,side)
            with open(inputFile, 'rb') as f:
                myMat_right = pickle.load(f)

            side = 'left'
            inputFile = annotationFolder  + 'extendedGene_labelCover_%d_%s.pkl' %(ev,side)
            with open(inputFile, 'rb') as f:
                myMat_left = pickle.load(f)

            myMat = myMat_left + myMat_right# coverage of labels

            enhancerCover = np.zeros(myMat.shape[0]) # total coverage of enhancer regions
            enhancerLowCover = np.zeros(myMat.shape[0]) # total coverage of enhancer regions
            for label in range(len(label_term_mapping)):
                if label_term_mapping[str(label)] == 'Enhancer':
                    enhancerCover += myMat[:, label]
                if label_term_mapping[str(label)] == 'EnhancerLow':
                    enhancerLowCover += myMat[:,label]


            #enhSortedArg = np.argsort(enhancerCover)[::-1]
            enhSortedArg = np.argsort(enhancerLowCover+enhancerCover)[::-1] # sum of both coverage
        
            xwalk = 0
            ywalk = 0
            area = 0
            aucc = np.zeros((len(expArray), 2))
            for j in range(len(expArray)):
                if expArray[enhSortedArg[j]] == 0:
                    xwalk +=1
                    area+= ywalk
                else:
                    ywalk +=1

                aucc[j, 0] = xwalk
                aucc[j, 1] = ywalk

            segwayAUCs[count, e] = area / (xwalk*ywalk)
            aucOutput[('segway', ev)] = aucc

        # >>>>>>>>>> Chrom, ev = 2000
        model = 'chrom'
        chmmFileName = annotation['chromFile'].split('/')[-1]

        if (('38batch' in annotation['folder']) or ('May11' in annotation['folder'])):
            annFile = annotationFolder + 'sorted_' + chmmFileName
        else:
            annFile = annotationFolder + chmmFileName
            
        for e, ev in enumerate(evList):
            side = 'right'

            inputFile = annotationFolder  + 'chrom_extendedGene_labelCover_%d_%s.pkl' %(ev,side)
            with open(inputFile, 'rb') as f:
                myMat_right = pickle.load(f)

            side = 'left'
            inputFile = annotationFolder  + 'chrom_extendedGene_labelCover_%d_%s.pkl' %(ev,side)
            with open(inputFile, 'rb') as f:
                myMat_left = pickle.load(f)

            myMat = myMat_left + myMat_right# coverage of labels
            book = np.sum(myMat, axis=1)

            enhancerCover = np.sum(myMat[:, [0,1,3,4,5]], axis=1)
            #enhancerLowCover = myMat[:, 5]

            enhSortedArg = np.argsort(enhancerCover[0:24700])[::-1]
            #enhSortedArg = np.argsort(enhancerLowCover[0:24700])[::-1]
        
            xwalk = 0
            ywalk = 0
            area = 0
            aucc = np.zeros((len(expArray), 2))
            for j in range(len(expArray)):
                if expArray[enhSortedArg[j]] == 0:
                    xwalk +=1
                    area+= ywalk
                else:
                    ywalk +=1

                aucc[j, 0] = xwalk
                aucc[j, 1] = ywalk

            chromAUCs[count, e] = area / (xwalk*ywalk)
            aucOutput[('chrom', ev)] = aucc

        count+=1

    aucCurves[accession] = aucOutput

plt.boxplot(segwayAUCs)
plt.ylim([.45, .85])
plt.title('segway enhancerLow expression prediction - 2k, 5k, 10k')
figFile = plotFolder + 'segwayEnhancerLow_exp_prediction.pdf'
print(figFile)
plt.tight_layout()
plt.savefig(figFile, bbox_inches='tight')
plt.close('all')

plt.boxplot(chromAUCs)
plt.ylim([.45, .85])
plt.title('chrom enhancerLow expression prediction - 2k, 5k, 10k')
figFile = plotFolder + 'chromEnhancerLow_exp_prediction.pdf'
print(figFile)
plt.tight_layout()
plt.savefig(figFile, bbox_inches='tight')
plt.close('all')

plt.show()

result = {'chromAUCs': chromAUCs, 'segwayAUCs': segwayAUCs, 'aucCurves': aucCurves}
file = dataFolder + 'enhancerLow_geneExp_output.pkl'
with open(file, 'wb') as f:
    pickle.dump(result, f)

file = dataFolder + 'enhancerLowAndEnhancer_geneExp_output.pkl'
with open(file, 'wb') as f:
    pickle.dump(result, f)

annotationResults = result

myCurve = aucOutput[('segway', 2000)]
plt.plot()

# loading the enhancer results
file = dataFolder + 'enhancer_geneExp_output.pkl'
with open(file, 'rb') as f:
    annotationResults = pickle.load(f)

segAUCs = annotationResults['segwayAUCs']
segAUCs = result['segwayAUCs']
print(np.median(segAUCs, axis=0))
print(np.max(segAUCs, axis=0))

chromAUCs = annotationResults['chromAUCs']
chromAUCs = result['chromAUCs']
print(np.median(chromAUCs, axis=0))
print(np.max(chromAUCs, axis=0))


# loading the enhancer and enhancer low results
file = dataFolder + 'enhancerLowAndEnhancer_geneExp_output.pkl'
with open(file, 'rb') as f:
    annotationResults_low = pickle.load(f)

segAUCs_low = annotationResults_low['segwayAUCs']
print(np.median(segAUCs, axis=0))
print(np.max(segAUCs, axis=0))

chromAUCs_low = annotationResults_low['chromAUCs']
print(np.median(chromAUCs, axis=0))
print(np.max(chromAUCs, axis=0))

plt.boxplot(segAUCs_low)
plt.ylim([.45, .85])
plt.title('segway enhancerLow expression prediction - 2k, 5k, 10k')
plt.show()
plt.close('all')

plt.boxplot(chromAUCs_low)
plt.ylim([.45, .85])
plt.title('chrom enhancerLow expression prediction - 2k, 5k, 10k')
plt.show()
plt.close('all')

# getting all the pairwise pvalues:

from scipy.stats import ttest_ind
from itertools import combinations

data = np.hstack((segAUCs, segAUCs_low, chromAUCs, chromAUCs_low))

column_indices = range(data.shape[1])
pairs = list(combinations(column_indices, 2))

pvalues = {}
for i, j in pairs:
    stat, p = ttest_ind(data[:, i], data[:, j], equal_var=False)
    pvalues[(i,j)] = p

n_cols = data.shape[1]

# Compute pairwise Welch’s t-tests
p_matrix = np.full((n_cols, n_cols), np.nan)
for i, j in combinations(range(n_cols), 2):
    _, p = ttest_ind(data[:, i], data[:, j], equal_var=False)
    p_matrix[i, j] = p
    p_matrix[j, i] = p  # symmetry

groups = ['seg', 'seg_low', 'chr', 'chr_low']
subgroups = [2, 5, 10]
# Create column labels 
labels = [f'{grp}{i}' for grp in groups for i in subgroups]

# Plot heatmap
plt.figure(figsize=(12, 10))
ax = sns.heatmap(p_matrix, 
                 annot=True, fmt=".2e", cmap="viridis", 
                 xticklabels=labels, yticklabels=labels,
                 cbar_kws={'label': 'p-value'},
                 linewidths=0.5, linecolor='gray')

# Set titles and layout
ax.set_title("Pairwise Welch's t-test p-values", fontsize=14)
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0)
plt.tight_layout()
plt.show()

figFile = plotFolder + 'seg_chromEnhExp_pvalues.pdf'
print(figFile)
plt.tight_layout()
plt.savefig(figFile, bbox_inches='tight')
plt.close('all')

########################################
# 2. For the set of example samples, get the bigwig file
########################################

# 2.1 get the accession name for the files
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
for i in indexList:
    print(accessionList[i])


#indexList = [1, 13, 52, 86, 94, 119, 123, 129, 132, 139, 233, 112, 67, 56 ,11, 9, 5]
#indexList = [2, 3, 4, 8, 11, 15, 16]
indexList = [1, 13, 52, 86, 94, 119, 123, 129, 132, 139, 233, 112, 67, 56 ,11, 9, 5, 2, 3, 4, 8, 11, 15, 16, 32, 64, 70, 119, 148, 153]
indexList = [1, 13, 52, 86, 94, 119, 123, 129, 132, 139, 233, 112, 67, 56 ,11, 9, 5, 2, 3, 4, 8, 11, 15, 16,  32, 64, 70, 119, 148, 153, 56, 94, 12, 127, 135, 144, 158, 163, 172, 176, 181, 182, 194, 196, 201, 210]
for accessionIndex in indexList:
    accession = accessionList[accessionIndex]
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']
    
    # >>>>>>>>>>>>>>>>>>>>>>>>>>
    label_term_mapping = {}
    mnemonics_file = annotationFolder + 'mnemonics_v04.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term

    # get the track file - H3K4me1
    assay_file = annotationFolder + 'trackname_assay.txt'
    with open(assay_file, 'r') as assays:
        print('------')
        for line in assays:
            if line.split()[1] == 'H3K4me1':
                print(line)
                h3k4 = line.split()[0]
                print(accession)
            if line.split()[1] == 'H3K27ac':
                print(line)
                h3k27 = line.split()[0]
                print(accession)


# 2.2 run the function for each of the files
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

for accessionIndex in indexList[6:]:
    accession = accessionList[accessionIndex]
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']

    if not(annotation['RNAseqFile'] == 'none') and not(annotation['chromFile'] == 'none'):

        assay_file = annotationFolder + 'trackname_assay.txt'
        with open(assay_file, 'r') as assays:
            for line in assays:
                if line.split()[1] == 'H3K4me1':
                    print(line)
                    h3k4 = line.split()[0]
                    print(accession)

        histoneFile = '/Users/marjanfarahbod/Downloads/%s.bigWig' %(h3k4)
        print(histoneFile)

        # >>>>>>>>>> Segway, ev = 2000
        annFile = annotation['bedFile']
        mnemFile = annotationFolder + 'mnemonics_v04.txt'
        ev = 2000
        side = 'right'
        EPLabelCover_Segway_extensionValue_histone(annotationFolder,
                                           annFile, mnemFile, geneList, geneIDList, ev, side, histoneFile)

        annFile = annotation['bedFile']
        mnemFile = annotationFolder + 'mnemonics_v04.txt'
        ev = 2000
        side = 'left'
        EPLabelCover_Segway_extensionValue_histone(annotationFolder,
                                                   annFile, mnemFile, geneList, geneIDList, ev, side, histoneFile)

        print(ev, accession)
        # >>>>>>>>>> Segway, ev = 5000
        annFile = annotation['bedFile']
        mnemFile = annotationFolder + 'mnemonics_v04.txt'
        ev = 5000
        side = 'right'
        EPLabelCover_Segway_extensionValue_histone(annotationFolder,
                                                   annFile, mnemFile, geneList, geneIDList, ev, side, histoneFile)

        annFile = annotation['bedFile']
        mnemFile = annotationFolder + 'mnemonics_v04.txt'
        ev = 5000
        side = 'left'
        EPLabelCover_Segway_extensionValue_histone(annotationFolder,
                                           annFile, mnemFile, geneList, geneIDList, ev, side, histoneFile)

        # >>>>>>>>>> Segway, ev = 10000
        annFile = annotation['bedFile']
        mnemFile = annotationFolder + 'mnemonics_v04.txt'
        ev = 10000
        side = 'right'
        EPLabelCover_Segway_extensionValue_histone(annotationFolder,
                                                   annFile, mnemFile, geneList, geneIDList, ev, side, histoneFile)

        annFile = annotation['bedFile']
        mnemFile = annotationFolder + 'mnemonics_v04.txt'
        ev = 10000
        side = 'left'
        EPLabelCover_Segway_extensionValue_histone(annotationFolder,
                                                   annFile, mnemFile, geneList, geneIDList, ev, side, histoneFile)

        
        # >>>>>>>>>> Chrom, ev = 2000
        chmmFileName = annotation['chromFile'].split('/')[-1]

        if (('38batch' in annotation['folder']) or ('May11' in annotation['folder'])):
            annFile = annotationFolder + 'sorted_' + chmmFileName
        else:
            annFile = annotationFolder + chmmFileName
            
        ev = 2000
        side = 'right'
        EPLabelCover_chromHMM_extensionValue_histone(annotationFolder, annFile,
                                             chromLabels, geneList, geneIDList, ev, side, histoneFile)
        # chromhmm has missing regions in the genome. Its coverage is ...
                    
        ev = 2000
        side = 'left'
        EPLabelCover_chromHMM_extensionValue_histone(annotationFolder, annFile,
                                                     chromLabels, geneList, geneIDList, ev, side, histoneFile)

        # >>>>>>>>>> Chrom, ev = 5000
                    
        ev = 5000
        side = 'right'
        EPLabelCover_chromHMM_extensionValue_histone(annotationFolder, annFile,
                                                     chromLabels, geneList, geneIDList, ev, side, histoneFile)
        # chromhmm has missing regions in the genome. Its coverage is ...
                    
        ev = 5000
        side = 'left'
        EPLabelCover_chromHMM_extensionValue_histone(annotationFolder, annFile,
                                                     chromLabels, geneList, geneIDList, ev, side, histoneFile)

        # >>>>>>>>>> Chrom, ev = 10000
                    
        ev = 10000
        side = 'right'
        EPLabelCover_chromHMM_extensionValue_histone(annotationFolder, annFile,
                                                     chromLabels, geneList, geneIDList, ev, side, histoneFile)
        # chromhmm has missing regions in the genome. Its coverage is ...
                    
        ev = 10000
        side = 'left'
        EPLabelCover_chromHMM_extensionValue_histone(annotationFolder, annFile,
                                                     chromLabels, geneList, geneIDList, ev, side, histoneFile)


# 2.3 get the output
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

switchLists = [['segway','chrom'], [2000, 5000, 10000]]

modeKeyList = []
for i in range(2):
    model = switchLists[0][i]
    for j in range(3):
        ev = switchLists[1][j]
        myKey = (model, ev)
        modeKeyList.append(myKey)

count = 0
self = True
segwayAUCs = np.zeros((len(indexList), 3))
chromAUCs = np.zeros((len(indexList), 3))
aucCurves = {}
for accessionIndex in indexList:
    accession = accessionList[accessionIndex]
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']

    aucOutput = {}
    # >>>>>>>>>>>>>>>>>>>>>>>>>>
    if self:
        rnaAccession = accession
        rnaAnnotation = annotation
    else:
        temp = random.randint(0,87)
        rnaAInd = rnaList[temp]
        if rnaAInd == accessionIndex:
            rnaAInd = rnaList[temp+1]

        rnaAccession = accessionList[rnaAInd]
        rnaAnnotation = allMeta[rnaAccession]
        
    rnaFolder = rnaAnnotation['folder']
                
    if len(rnaAnnotation['RNAseqFile'][0]) > 1:
        RNAFile = rnaAnnotation['RNAseqFile'][0]
    else:
        RNAFile = rnaAnnotation['RNAseqFile']
        print(RNAFile)
            
    expAccession = RNAFile[-15:-4]
    expFile = rnaFolder +  'geneExp_dict_' + expAccession + '.pkl'
    print(expFile)
    print(annotationFolder)

    with open(expFile, 'rb') as pickledFile:
        expression = pickle.load(pickledFile) # load the expression file

    expArray = np.asarray([expression[x] for x in geneIDList])[0:24700] # genes in the transcriptomic data that a

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # segway
    label_term_mapping = {}
    mnemonics_file = annotationFolder + 'mnemonics_v04.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term
                
    for e, ev in enumerate(evList):
        side = 'right'
        inputFile = annotationFolder  + 'extendedGene_labelCover_%d_%s_hvals.pkl' %(ev,side)
        with open(inputFile, 'rb') as f:
            myMat_right = pickle.load(f)

        side = 'left'
        inputFile = annotationFolder  + 'extendedGene_labelCover_%d_%s_hvals.pkl' %(ev,side)
        with open(inputFile, 'rb') as f:
            myMat_left = pickle.load(f)

        myMat = myMat_left['histoneVals'] + myMat_right['histoneVals']# coverage of labels
        normMat = myMat_left['expMat'] + myMat_right['expMat']

        enhancerValue = np.zeros(myMat.shape[0]) # total coverage of enhancer regions
        enhancerCover = np.zeros(myMat.shape[0]) # total coverage of enhancer regions
        enhancerLowValue = np.zeros(myMat.shape[0]) # total coverage of enhancer regions
        enhancerLowCover = np.zeros(myMat.shape[0]) # total coverage of enhancer regions
        for label in range(len(label_term_mapping)):
            if label_term_mapping[str(label)] == 'Enhancer':
                enhancerValue += myMat[:, label]
                enhancerCover += normMat[:, label]
            if label_term_mapping[str(label)] == 'EnhancerLow':
                enhancerLowValue += myMat[:,label]
                enhancerLowCover += normMat[:,label]

        enhCoverNorm = np.zeros(len(enhancerCover))
        for k in range(len(enhCoverNorm)):
            enhCoverNorm[k] = enhancerValue[k] / (enhancerCover[k]+1)
            
        enhSortedArg = np.argsort(enhancerCover)#[::-1]
        enhNormSortedArg = np.argsort(enhCoverNorm)[::-1]
        expSortedArg = np.argsort(expArray)

        xwalk = 0
        ywalk = 0
        area = 0
        aucc = np.zeros((len(expArray), 2))
        for j in range(len(expArray)):
            #if expArray[enhSortedArg[j]] == 0:
            if expArray[enhNormSortedArg[j]] == 0:
                xwalk +=1
                area+= ywalk
            else:
                ywalk +=1

            aucc[j, 0] = xwalk
            aucc[j, 1] = ywalk

        segwayAUCs[count, e] = area / (xwalk*ywalk)
        aucOutput[('segway', ev)] = aucc

        # >>>>>>>>>> Chrom, ev = 2000
        chmmFileName = annotation['chromFile'].split('/')[-1]

        if (('38batch' in annotation['folder']) or ('May11' in annotation['folder'])):
            annFile = annotationFolder + 'sorted_' + chmmFileName
        else:
            annFile = annotationFolder + chmmFileName
            
        for e, ev in enumerate(evList):
            side = 'right'

            inputFile = annotationFolder  + 'chrom_extendedGene_labelCover_%d_%s_hvals.pkl' %(ev,side)
            with open(inputFile, 'rb') as f:
                myMat_right = pickle.load(f)

            side = 'left'
            inputFile = annotationFolder  + 'chrom_extendedGene_labelCover_%d_%s_hvals.pkl' %(ev,side)
            with open(inputFile, 'rb') as f:
                myMat_left = pickle.load(f)

            myMat = myMat_left['histoneVals'] + myMat_right['histoneVals']# coverage of labels
            normMat = myMat_left['expMat'] + myMat_right['expMat']

            enhancerValue = np.sum(myMat[:, [0,1,3,4]], axis=1)
            enhancerLowValue = myMat[:, 5]
            enhancerCover = np.sum(normMat[:, [0,1,3,4]], axis=1)
            enhancerLowCover = normMat[:, 5]
            
            enhCoverNorm = np.zeros(len(enhancerCover))
            for k in range(len(enhCoverNorm)):
                enhCoverNorm[k] = enhancerValue[k] / (enhancerCover[k]+1)
            
            enhNormSortedArg = np.argsort(enhCoverNorm[0:24700])[::-1]
            enhSortedArg = np.argsort(enhancerCover[0:24700])[::-1]
        
            xwalk = 0
            ywalk = 0
            area = 0
            aucc = np.zeros((len(expArray), 2))
            for j in range(len(expArray)):
                if expArray[enhNormSortedArg[j]] == 0:
                    xwalk +=1
                    area+= ywalk
                else:
                    ywalk +=1

                aucc[j, 0] = xwalk
                aucc[j, 1] = ywalk

            chromAUCs[count, e] = area / (xwalk*ywalk)
            aucOutput[('chrom', ev)] = aucc
    count+=1

    aucCurves[accession] = aucOutput

result = {'chromAUCs': chromAUCs, 'segwayAUCs': segwayAUCs, 'aucCurves': aucCurves}
file = dataFolder + 'enhancer_geneExp_output_histoneVals.pkl'
with open(file, 'wb') as f:
    pickle.dump(result, f)

annotationResults_hvals = result
plt.plot(annotationResults_hvals[('segway', 2000)])

# We have the plots for this
print(annotationResults_hvals['chromAUCs'])
print(annotationResults_hvals['segwayAUCs'])

print(annotationResults['chromAUCs'])
print(annotationResults['segwayAUCs'])

#
accessionIndex = 1
accession = accessionList[accessionIndex]

plt.plot(enhancerCover[expSortedArg])
plt.scatter(range(len(enhSortedArg)), enhancerCover[expSortedArg])
plt.scatter(range(len(enhNormSortedArg)), expArray[enhNormSortedArg])
plt.show()

aucCurves = aucCurves[accession][('segway', 2000)]
plt.plot(aucCurves[:, 0], aucCurves[:,1])
plt.show()

plt.boxplot(annotationResults_hvals['chromAUCs'])
plt.ylim([.45, .85])
plt.title('Chrom enhancer expression prediction - histone values - 2k, 5k, 10k')
figFile = plotFolder + 'ChromEnhancer_exp_prediction_histone.pdf'
print(figFile)
plt.tight_layout()
plt.savefig(figFile, bbox_inches='tight')
plt.close('all')

plt.boxplot(annotationResults_hvals['segwayAUCs'])
plt.ylim([.45, .85])
plt.title('segway enhancer expression prediction - histone values - 2k, 5k, 10k')
figFile = plotFolder + 'SegwayEnhancer_exp_prediction_histone.pdf'
print(figFile)
plt.tight_layout()
plt.savefig(figFile, bbox_inches='tight')
plt.close('all')


########################################
# 3. The bar plots
########################################

# 3.1.0 Segway: get the chr20 section of the .bed file
# did this for both chr19 and chr20
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# for the list of samples, get the chr19.bed file
for accessionIndex in indexList[24:]:
    accession = accessionList[accessionIndex]
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']

    annFile = annotation['bedFile']

    if annFile.endswith('.gz'):
        os.system('gunzip %s' %(annFile))
    else:
        os.system('gunzip %s.gz' %(annFile))

    command = "grep -E 'chr19.*' %s" %(annFile)
    print(command)
    
    out = annotation['folder'] + 'chr19.bed'
    f = open(out, 'w')
    import subprocess
    subprocess.run(command, shell=True, stdout=f)

    os.system('gzip %s' %(annFile))
 
# for the list of samples, get the chr20.bed and 17 file
for accessionIndex in indexList:
    accession = accessionList[accessionIndex]
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']

    annFile = annotation['bedFile']

    if annFile.endswith('.gz'):
        os.system('gunzip %s' %(annFile))
    else:
        os.system('gunzip %s.gz' %(annFile))

    command = "grep -E 'chr18.*' %s" %(annFile)
    print(command)
    
    out = annotation['folder'] + 'chr18.bed'
    f = open(out, 'w')
    import subprocess
    subprocess.run(command, shell=True, stdout=f)

    os.system('gzip %s' %(annFile))


# 3.1.1 Chrom: get the chr17 section of the .bed file
# I did it for chr19 too
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  
import subprocess
for accessionIndex in indexList:
    accession = accessionList[accessionIndex]
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']

    chmmFileName = annotation['chromFile'].split('/')[-1]
    print(chmmFileName)
    
    if (('38batch' in annotation['folder']) or ('May11' in annotation['folder'])):
        annFile = annotationFolder + 'sorted_' + chmmFileName
    else:
        annFile = annotationFolder + chmmFileName

    if annFile.endswith('.gz'):
        os.system('gunzip %s' %(annFile))

        command = "grep -E 'chr18.*' %s" %(annFile[0:-3])
        print(command)
    else:
        os.system('gunzip %s.gz' %(annFile))

        command = "grep -E 'chr18.*' %s" %(annFile)
        print(command)

    out = annotation['folder'] + 'chr18_chrom.bed'
    f = open(out, 'w')

    subprocess.run(command, shell=True, stdout=f)

    print(annFile)
    if annFile.endswith('.gz'):
        os.system('gzip %s' %(annFile[0:-3]))
    else:
        os.system('gzip %s' %(annFile))

# 3.2.0 Segway colors
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
colorStateMap = {'Bivalent' : np.asarray([189, 183, 107])/256,
                 'Transcribed' : np.asarray([0, 128, 0])/256,
                 'CTCF' : np.asarray([196, 225, 5])/256,
                 'K9K36' : np.asarray([102, 205, 170])/256,
                 'Enhancer' : np.asarray([255, 195, 77])/256,
                 'PromoterFlanking' : np.asarray([255, 68, 0])/256,
                 'ConstitutiveHet' : np.asarray([138, 145, 208])/256,
                 'FacultativeHet' : np.asarray([128, 0, 128])/256,
                 'Promoter' : np.asarray([255, 0, 0])/256,
                 'Quiescent' : np.asarray([190, 190, 190])/256,
                 'EnhancerLow' : np.asarray([255, 255, 0])/256}

# 3.2.1 chrom colors
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
colorStateMapChrom = {'11' : np.asarray([189, 183, 107])/256,
                 '15' : np.asarray([0, 128, 0])/256,
                      '16' : np.asarray([0, 200, 0])/256,
                 '17' : np.asarray([102, 205, 170])/256,
                 '0' : np.asarray([255, 195, 77])/256,
                '1' : np.asarray([255, 195, 77])/256,
                      '2' :np.asarray([189, 183, 107])/256,
                      '3' : np.asarray([255, 195, 77])/256,
                      '4' : np.asarray([255, 195, 77])/256,
                 '12' : np.asarray([255, 68, 0])/256,
                      '13' : np.asarray([255, 68, 0])/256,
                      '14' : np.asarray([255, 68, 0])/256,
                 '6' : np.asarray([138, 145, 208])/256,
                 '8' : np.asarray([128, 0, 128])/256,
                 '9' : np.asarray([128, 0, 128])/256,                      
                 '10' : np.asarray([255, 0, 0])/256,
                 '7' : np.asarray([180, 180, 180])/256,
                      '5' : np.asarray([255, 255, 0])/256,
                      '18': np.asarray([0,0,0])/256}


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> just getting list of histone files
hfileList= []
for accessionIndex in indexList:
    accession = accessionList[accessionIndex]
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']

    assay_file = annotationFolder + 'trackname_assay.txt'
    with open(assay_file, 'r') as assays:
        for line in assays:
            #if line.split()[1] == 'H3K4me1':
            if line.split()[1] == 'H3K27ac':
                print(line)
                #h3k4 = line.split()[0]
                h3k27 = line.split()[0]
                print(accession)

    #hfileList.append('/Users/marjanfarahbod/Downloads/%s.bigWig' %(h3k4))
    #hfileList.append('%sbigwigs/%s.bigWig' %(dataFolder, h3k4))
    hfileList.append('%sbigwigs/%s.bigWig' %(dataFolder, h3k27))
    #print(histoneFile)

# 3.3.0 Segway, get the plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# note : lines with #OOO is where the chrom number should change in the code

import os


endInd = 64000000 # chr20 # 20 does not have many genes 
endInd = 58000000 # chr19
endInd = 83000000 # chr17
endInd = 80000000 # chr18 

enhCoverage = np.zeros((len(indexList), 5))
enhLowCoverage = np.zeros((len(indexList), 5))
genomeHistoneCove = np.zeros((len(indexList), 5))
allActiveCoverage = np.zeros((len(indexList), 5))

for a, accessionIndex in enumerate(indexList):
    accession = accessionList[accessionIndex]
    print(accession)
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']

    assay_file = annotationFolder + 'trackname_assay.txt'
    with open(assay_file, 'r') as assays:
        for line in assays:
            #if line.split()[1] == 'H3K4me1':
            if line.split()[1] == 'H3K27ac':
                print(line)
                h3k27 = line.split()[0]
                #h3k4 = line.split()[0]
                #print(accession)

    #histoneFile = '/Users/marjanfarahbod/Downloads/%s.bigWig' %(h3k4)
    histoneFile = hfileList[a]
    print(histoneFile)
    #randInd = random.randint(0, len(indexList))
    #histoneFile = hfileList[randInd]
    #print(histoneFile)

    if not(os.path.exists(histoneFile)):
        print(histoneFile)
        print('file does not exisst')
        print(a)
        print(accessionIndex)
    
    bw = pyBigWig.open(histoneFile)
    vals = bw.values('chr19', 0, endInd) #OOO

    vals=np.asarray(vals)
    print('vals In')
    
    meanSig = np.zeros((int(endInd/100)))
    for i in range(len(meanSig-1)):
        #print(i*100, (i+1)*100)
        meanSig[i] = np.mean(vals[i*100:(i+1)*100])
    
    mnemFile = annotationFolder + 'mnemonics_v04.txt'
    label_term_mapping = {}
    with open(mnemFile, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term
            
    chr19File = annotation['folder'] + 'chr19.bed' #OOO
    segStates = np.zeros(len(meanSig))
    i = 0
    annLineInd = 1
    sumStates = np.zeros(len(label_term_mapping))
    while i < len(segStates):
        line = linecache.getline(chr19File, annLineInd)  #OOO
        annStart = int(line.split('\t')[1])
        annEnd = int(line.split('\t')[2])
        annState = int(line.split('\t')[3].split('_')[0])
        step = int((annEnd - annStart) / 100)
        segStates[i:i+step] = annState
        i = i + step
        annLineInd +=1
        sumStates[annState] += annEnd - annStart

    segStatesToTerms = np.zeros(len(meanSig))
    for i in range(len(segStates)):
        segStatesToTerms[i] = segwayStates.index(label_term_mapping[str(int(segStates[i]))])

    # color the scatter hahahaha
    colorString = []
    for state in segStatesToTerms:
        colorString.append(colorStateMap[label_term_mapping[str(int(state))]])

    # meanSig varies from 0 to 4 and greater. It is fold change, so we care about 0, 1, 2, 3, 4
    stateCoverage = np.zeros((5, 11))
    meanSigCopy = meanSig.copy()
    meanSigCopy[meanSigCopy >= 5] = 4.99
    for i in range(5):
        #inds = np.where(meanSigCopy > i & meanSigCopy<=(i+1)[0])
        inds = np.where(np.logical_and(meanSigCopy >= i, meanSigCopy<(i+1)))[0]
        myStates = segStatesToTerms[inds]
        for j in range(11):
            stateCoverage[i, j] = np.sum(myStates==j)

    stateCoverage = stateCoverage/np.sum(stateCoverage, axis=1)[:, np.newaxis]
    enhCoverage[a, :] = stateCoverage[:, 2]
    enhLowCoverage[a, :] = stateCoverage[:, 3]
    allActiveCoverage[a, :] = np.sum(stateCoverage[:, 0:7], axis=1)

    # plots 
    for i in range(11):
        print(i)
        if i >0:
            plt.bar(range(5), stateCoverage[:,i], bottom=np.sum(stateCoverage[:, 0:i], axis=1), color=colorStateMap[segwayStates[i]])
        else:
            plt.bar(range(5), stateCoverage[:,i], color=colorStateMap[segwayStates[i]])
        print(stateCoverage[:,i])

    titleString = ''
    for i in range(5):
        val = np.sum((np.logical_and(meanSigCopy >= i, meanSigCopy<(i+1)))/len(meanSigCopy))
        print(i, val)
        #ratioList.append(np.sum((np.logical_and(meanSigCopy >= i, meanSigCopy<(i+1)))/len(meanSigCopy)))
        titleString = titleString + '   %.3f' %(val)
        genomeHistoneCove[a, i] = val
    #plt.close('all')

    #figFile = plotFolder + 'segwayLabelCoverage_H3K4me1_%s_chr19_v02.pdf' %(accession) #OOO
    figFile = plotFolder + 'segwayLabelCoverage_H3K27ac_%s_chr19_v02.pdf' %(accession) #OOO
    print(figFile)
    plt.title(titleString)
    plt.tight_layout()
    plt.savefig(figFile, bbox_inches='tight')
    plt.close('all')

segRes = {}
segRes['enhCoverage'] = enhCoverage
segRes['enhLowCoverage'] = enhLowCoverage
segRes['allActiveCoverage'] = allActiveCoverage
#file = dataFolder + 'H3K4me1_enhCoverage_segway_extras_chr18.pkl'  #OOO
file = dataFolder + 'H3K27ac_enhCoverage_segway_extras_chr19.pkl'  #OOO
with open(file, 'wb') as f:
    pickle.dump(segRes, f)

file = dataFolder + 'H3K4me1_enhCoverage_segway.pkl' 
with open(file, 'rb') as f:
    segRes=pickle.load(f)

file = dataFolder + 'H3K4me1_enhCoverage_chrom.pkl'
with open(file, 'rb') as f:
    chromRes=pickle.load(f)

file = dataFolder + 'H3K4me1_enhCoverage_segway_extras.pkl'
with open(file, 'rb') as f:
    segResX=pickle.load(f)

file = dataFolder + 'H3K4me1_enhCoverage_chrom_extras.pkl'
with open(file, 'rb') as f:
    chromResX=pickle.load(f)


# do the t-test and print the pvalues
from scipy import stats
#stats.ttest_ind(segwayGeneAUC, chromGeneAUC, equal_var=False)

coverage = ['enhCoverage', 'enhLowCoverage', 'allActiveCoverage']
ps = np.zeros((3, 5))
for i,c in enumerate(coverage):
    for j in range(5):
        #ps[i, j] = stats.ttest_ind(chromRes[c][:,j], segRes[c][:,j], equal_var=False)[1]
        chromArr = np.concatenate((chromRes[c][:,j], chromResX[c][:,j]))
        segArr = np.concatenate((segRes[c][:,j], segResX[c][:,j]))

        #chromArr = chromRes[c][:,j]
        #segArr = segRes[c][:,j]
        ps[i, j] = stats.ttest_ind(chromArr, segArr)[1]

# save this to a file
myFile = dataFolder + 'histoneFigure3C_Epvalues.tsv'
np.savetxt(myFile, ps, delimiter='\t')


plt.boxplot(genomeHistoneCove)
#figFile = plotFolder + 'segwayLabelCoverage_H3K4me1Foldchange_genomeCoverage_chr19.pdf'
figFile = plotFolder + 'segwayLabelCoverage_H3K27acFoldchange_genomeCoverage_chr19.pdf' #OOO
print(figFile)
plt.tight_layout()
plt.savefig(figFile, bbox_inches='tight')
plt.close('all')

enhCoverage = segRes['enhCoverage']
plt.boxplot(enhCoverage)
plt.ylim((-.1, 1.1))
#figFile = plotFolder + 'segwayLabelCoverage_H3K4me1_chr18.pdf' #OOO
figFile = plotFolder + 'segwayLabelCoverage_H3K27ac_chr19.pdf' #OOO
print(figFile)
plt.tight_layout()
plt.savefig(figFile, bbox_inches='tight')
plt.close('all')

enhLowCoverage = segRes['enhLowCoverage']
plt.boxplot(enhLowCoverage)
plt.ylim((-.1, 1.1))
#figFile = plotFolder + 'segwayLabelCoverage_H3K4me1_chr18_enhLow.pdf'
figFile = plotFolder + 'segwayLabelCoverage_H3K27ac_chr19_enhLow.pdf' #OOO
print(figFile)
plt.tight_layout()
plt.savefig(figFile, bbox_inches='tight')
plt.close('all')

allActiveCoverage = segRes['allActiveCoverage']
plt.boxplot(allActiveCoverage)
plt.ylim((-.1, 1.1))
#figFile = plotFolder + 'segwayLabelCoverage_H3K4me1_chr18allActive.pdf'
figFile = plotFolder + 'segwayLabelCoverage_H3K27ac_chr19allActive.pdf' 
print(figFile)
plt.tight_layout()
plt.savefig(figFile, bbox_inches='tight')
plt.close('all')


# 3.3.1 Chrom, get the plot 
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# endInd = 64000000 # chr20
endInd = 58000000 # chr19
endInd = 80000000 # chr18
enhCoverage = np.zeros((len(indexList), 5))
enhLowCoverage = np.zeros((len(indexList), 5))
genomeHistoneCove = np.zeros((len(indexList), 5))
allActiveCoverage = np.zeros((len(indexList), 5))

for a, accessionIndex in enumerate(indexList):
    accession = accessionList[accessionIndex]
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']

    # get histone file
    # >>>>>>>>>>
    #histoneFile = '/Users/marjanfarahbod/Downloads/%s.bigWig' %(h3k4)
    histoneFile = hfileList[a]
    #randInd = random.randint(0, len(indexList))
    #histoneFile = hfileList[randInd]
    print(histoneFile)
    
    bw = pyBigWig.open(histoneFile)
    vals = bw.values('chr19', 0, endInd) #OOO

    vals = np.asarray(vals)
    print('vals In')
    
    meanSig = np.zeros((int(endInd/100)))
    for i in range(len(meanSig-1)):
        meanSig[i] = np.mean(vals[i*100:(i+1)*100])

    # get annotation file
    # >>>>>>>>>>
    chr19FileChrom = annotation['folder'] + 'chr19_chrom.bed' #OOO

    i = 0
    annLineInd = 1
    chromStates = np.zeros(int(endInd/100)) + 18
    line = linecache.getline(chr19FileChrom, annLineInd) #OOO
    annStart = int(line.split('\t')[1])
    annEnd = int(line.split('\t')[2])
    annState = chromLabels.index((line.split('\t')[3]))

    while int(annEnd/100) < len(chromStates):
        si = int(annStart/100)
        ei = int(annEnd/100)+1
        chromStates[si:ei] = annState
        annLineInd +=1
        line = linecache.getline(chr19FileChrom, annLineInd) #OOO
        annStart = int(line.split('\t')[1])
        annEnd = int(line.split('\t')[2])
        annState = chromLabels.index((line.split('\t')[3]))

    '''
    pend = 0
    while i < len(chromStates):
        #print(line)
        line = linecache.getline(chr19FileChrom, annLineInd)
        annStart = int(line.split('\t')[1])
        annEnd = int(line.split('\t')[2])
        annState = chromLabels.index((line.split('\t')[3]))
        #>>>>

        #<<<<
        
        #print(annState)
        if not(annStart == pend):
            if pend ==  7350695:
                break
            step = int((annStart-pend) / 100)
            chromStates[i:i+step] = 18
            i = i+step
            print(annStart, pend)
            print(annLineInd)
            print(step)
            print('----------------')
        
        step = int((annEnd - annStart) / 100)
        chromStates[i:i+step] = annState
        i = i + step
        annLineInd +=1
        pend = annEnd
    '''

    linecache.clearcache()

    # color the scatter hahahaha
    colorStringChrom = []
    for state in chromStates:
        colorStringChrom.append(colorStateMapChrom[str(int(state))])

    orderList = [10, 12, 13, 14, 0, 1, 3, 4, 5, 2, 11, 15, 16, 17, 8, 9, 6, 7, 18] # MAKE THE TODO LIST AND MAKE IT PRETTY

    stateCoverageChrom = np.zeros((5, 19))
    meanSigCopy = meanSig.copy()
    meanSigCopy[meanSigCopy >= 5] = 4.99
    for i in range(5):
        #inds = np.where(meanSigCopy > i & meanSigCopy<=(i+1)[0])
        inds = np.where(np.logical_and(meanSigCopy >= i, meanSigCopy<(i+1)))[0]
        myStates = chromStates[inds]
        for j in range(19):
            stateCoverageChrom[i, j] = np.sum(myStates==j)

    stateCoverageChrom = stateCoverageChrom/np.sum(stateCoverageChrom, axis=1)[:, np.newaxis]
    enhCoverage[a, :] = np.sum(stateCoverageChrom[:, [0,1,3,4]], axis=1)
    enhLowCoverage[a, :] = stateCoverageChrom[:, 5]
    allActiveCoverage[a, :] = np.sum(stateCoverageChrom[:, [0,1,2,3,4,5,10,11,12,13,14,15,16,17]], axis=1)

    for i in range(19):
        print(i)
        if i >0:
            bottom = np.sum(stateCoverageChrom[:, orderList[0:i]], axis=1)
            plt.bar(range(5), stateCoverageChrom[:,orderList[i]], bottom=bottom, color=colorStateMapChrom[str(int(orderList[i]))])

        else:
            plt.bar(range(5), stateCoverageChrom[:, orderList[i]], color=colorStateMapChrom[str(int(orderList[i]))])
            print(stateCoverageChrom[:,orderList[i]])

    titleString = ''
    for i in range(5):
        val = np.sum((np.logical_and(meanSigCopy >= i, meanSigCopy<(i+1)))/len(meanSigCopy))
        print(i, val)
        #ratioList.append(np.sum((np.logical_and(meanSigCopy >= i, meanSigCopy<(i+1)))/len(meanSigCopy)))
        titleString = titleString + '   %.3f' %(val)
        genomeHistoneCove[a, i] = val

    #figFile = plotFolder + 'ChromLabelCoverage_H3K4me1_%s_chr18_corrected.pdf' %(accession) #OOO
    figFile = plotFolder + 'ChromLabelCoverage_H3K27ac_%s_chr19_corrected.pdf' %(accession) #OOO
    print(figFile)
    plt.tight_layout()
    plt.title(titleString)
    plt.savefig(figFile, bbox_inches='tight')
    plt.close('all')
    #aplt.show()

chromRes = {}
chromRes['enhCoverage'] = enhCoverage
chromRes['enhLowCoverage'] = enhLowCoverage
chromRes['allActiveCoverage'] = allActiveCoverage
#file = dataFolder + 'H3K4me1_enhCoverage_chrom_extras_chr18.pkl' #OOO
file = dataFolder + 'H3K27ac_enhCoverage_chrom_extras_chr19.pkl' #OOO
with open(file, 'wb') as f:
    pickle.dump(chromRes, f)

    
plt.boxplot(genomeHistoneCove)
#figFile = plotFolder + 'chromLabelCoverage_H3K4me1Foldchange_genomeCoverage_chr18.pdf' #OOO
figFile = plotFolder + 'chromLabelCoverage_H3K27acFoldchange_genomeCoverage_chr19.pdf' #OOO
print(figFile)
plt.tight_layout()
plt.savefig(figFile, bbox_inches='tight')
plt.close('all')

plt.boxplot(enhCoverage)
plt.ylim((-.1, 1.1))
#figFile = plotFolder + 'chromLabelCoverage_H3K4me1_chr18.pdf' #OOO
figFile = plotFolder + 'chromLabelCoverage_H3K27ac_chr19.pdf' #OOO
print(figFile)
plt.tight_layout()
plt.savefig(figFile, bbox_inches='tight')
plt.close('all')

plt.boxplot(enhLowCoverage)
plt.ylim((-.1, 1.1))
#figFile = plotFolder + 'chromLabelCoverage_H3K4me1_chr18enhLow.pdf'  #OOO
figFile = plotFolder + 'chromLabelCoverage_H3K27ac_chr19enhLow.pdf'  #OOO
print(figFile)
plt.tight_layout()
plt.savefig(figFile, bbox_inches='tight')
plt.close('all')


plt.boxplot(allActiveCoverage)
plt.ylim((-.1, 1.1))
#figFile = plotFolder + 'chromLabelCoverage_H3K4me1_chr18_allActive.pdf' #OOO
figFile = plotFolder + 'chromLabelCoverage_H3K27ac_chr19_allActive.pdf' #OOO
print(figFile)
plt.tight_layout()
plt.savefig(figFile, bbox_inches='tight')
plt.close('all')


########################################
# 4. The three section plot
########################################


matFilters = {}
# creat the filter
for a, accessionIndex in enumerate(indexList):
    accession = accessionList[accessionIndex]
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']
    for e, ev in enumerate(evList):
        myKey = (accession, ev)
        side = 'right'
        inputFile = annotationFolder  + 'extendedGene_labelCover_%d_%s_hvals.pkl' %(ev,side)
        with open(inputFile, 'rb') as f:
            myMat_right = pickle.load(f)

        side = 'left'
        inputFile = annotationFolder  + 'extendedGene_labelCover_%d_%s_hvals.pkl' %(ev,side)
        with open(inputFile, 'rb') as f:
            myMat_left = pickle.load(f)

        #myMat = myMat_left['histoneVals'] + myMat_right['histoneVals']# histone values, sum over bins of 100
        segMat = myMat_left['expMat'] + myMat_right['expMat']# coverage of labels
        segFilter = (np.sum(segMat,axis=1) > 0)

        side = 'right'
        inputFile = annotationFolder  + 'chrom_extendedGene_labelCover_%d_%s_hvals.pkl' %(ev,side)
        with open(inputFile, 'rb') as f:
            myMat_right = pickle.load(f)

        side = 'left'
        inputFile = annotationFolder  + 'chrom_extendedGene_labelCover_%d_%s_hvals.pkl' %(ev,side)
        with open(inputFile, 'rb') as f:
            myMat_left = pickle.load(f)

        #myMat = myMat_left['histoneVals'] + myMat_right['histoneVals']# coverage of labels
        chromMat = myMat_left['expMat'] + myMat_right['expMat']
        chromFilter = np.sum(chromMat,axis=1) > 0

        myFilter = np.logical_and(chromFilter[0:24700] , segFilter[0:24700])
        matFilters[myKey] = myFilter

self = True # true or false it can be
count = 0
evList = [2000, 5000, 10000]
for a, accessionIndex in enumerate(indexList):
    accession = accessionList[accessionIndex]
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']

    if not(annotation['RNAseqFile'] == 'none') and not(annotation['chromFile'] == 'none'):
        print(a)
        # Segway
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        label_term_mapping = {}
        mnemonics_file = annotationFolder + 'mnemonics_v04.txt'
        with open(mnemonics_file, 'r') as mnemonics:
            for line in mnemonics:
                #print(line)
                label = line.strip().split()[0]
                term = line.strip().split()[1]
                label_term_mapping[label] = term

        # >>>>>>>>>>>>>>>>>>>>>>>>>>
        if self:
            rnaAccession = accession
            rnaAnnotation = annotation
        else:
            temp = random.randint(0,87)
            rnaAInd = rnaList[temp]
            if rnaAInd == a:
                rnaAInd = rnaList[temp+1]
            rnaAccession = accessionList[rnaAInd]
            rnaAnnotation = allMeta[rnaAccession]
            
        rnaFolder = rnaAnnotation['folder']
                
        if len(rnaAnnotation['RNAseqFile'][0]) > 1:
            RNAFile = rnaAnnotation['RNAseqFile'][0]
        else:
            RNAFile = rnaAnnotation['RNAseqFile']
            print(RNAFile)
            
        expAccession = RNAFile[-15:-4]
        expFile = rnaFolder +  'geneExp_dict_' + expAccession + '.pkl'
        print(expFile)
        print(annotationFolder)

        with open(expFile, 'rb') as pickledFile:
            expression = pickle.load(pickledFile) # load the expression file

        expArray = np.asarray([expression[x] for x in geneIDList])[0:24700] # genes in the transcriptomic data that are in the gene map file
        
        # >>>>>>>>>> Segway evList
        fig, axs = plt.subplots(4,3, figsize=(8,8))
        for e, ev in enumerate(evList):
            myFilter = matFilters[(accession, ev)]
            expArray = np.asarray([expression[x] for x in geneIDList])[0:24700] # genes in the transcriptomic data that are in the gene map file
            
            side = 'right'
            inputFile = annotationFolder  + 'extendedGene_labelCover_%d_%s_hvals.pkl' %(ev,side)
            with open(inputFile, 'rb') as f:
                myMat_right = pickle.load(f)

            side = 'left'
            inputFile = annotationFolder  + 'extendedGene_labelCover_%d_%s_hvals.pkl' %(ev,side)
            with open(inputFile, 'rb') as f:
                myMat_left = pickle.load(f)

            myMat = myMat_left['histoneVals'] + myMat_right['histoneVals']# histone values, sum over bins of 100
            normMat = myMat_left['expMat'] + myMat_right['expMat']# coverage of labels

            #enhancerValue = np.zeros(myMat.shape[0]) # total value of enhancer regions
            enhancerCover = np.zeros(myMat.shape[0]) # total coverage of enhancer regions
            enhancerLowValue = np.zeros(myMat.shape[0]) # total coverage of enhancer regions
            enhancerLowCover = np.zeros(myMat.shape[0]) # total coverage of enhancer regions
            for label in range(len(label_term_mapping)):
                if label_term_mapping[str(label)] == 'Enhancer':
                    #enhancerValue += myMat[:, label]
                    enhancerCover += normMat[:, label]
                if label_term_mapping[str(label)] == 'EnhancerLow':
                    #enhancerLowValue += myMat[:,label]
                    enhancerLowCover += normMat[:,label]

            #enhancerCover = enhancerCover + enhancerLowCover # only this one when we want both enhlow and enh
            #enhancerValue = enhancerValue + enhancerLowValue # only this one when we want both enhlow and enh
            enhancerValue = np.sum(myMat, axis=1)/(ev*2) # sum value of H3K4me1
            
            enhancerValue = enhancerValue[myFilter]
            expArray = expArray[myFilter]
            enhancerCover = enhancerCover[myFilter]
            enhancerLowCover = enhancerLowCover[myFilter]

            '''
            enhCoverNorm = np.zeros(len(enhancerCover))
            for k in range(len(enhCoverNorm)):
                enhCoverNorm[k] = enhancerValue[k] / (enhancerCover[k]+1)
            '''
            
            #enhNormSortedArg = np.argsort(enhCoverNorm[0:24700])
            #enhNormSorted = np.sort(enhCoverNorm[0:24700])
            #enhSortedArg = np.argsort(enhancerCover[0:24700])
            #enhSorted = np.sort(enhancerCover[0:24700])
            enhValueSortedArg = np.argsort(enhancerValue) # >>> this is what we use for sorting
        
            geneBinSize = 100
            geneBin = 0
            geneBinCount = int(len(expArray)/geneBinSize)
            expRatios = np.zeros(geneBinCount) # ratio of genes expressed >>>> Expression
            #expMeans =  np.zeros(geneBinCount) # mean expression value
            #expMeansZ =  np.zeros(geneBinCount) # mean expression for exp genes
            #expMediansZ = np.zeros(geneBinCount) # median expression for exp genes

            # we need signal, coverage, expression in the plot
            enhMeans = np.zeros(geneBinCount) # mean enhancer values >>>> Signal
            enhCovers = np.zeros(geneBinCount) # mean enhancer coverage >>>> Coverage
            enhLowCovers = np.zeros(geneBinCount) # mean enhancer coverage >>>> Coverage Low
            #enhNormMean = np.zeros(geneBinCount) # mean enhancer value for coverage
            while geneBin < geneBinCount:
                #genes_inds = enhNormSortedArg[geneBin*geneBinSize:(geneBin+1)*geneBinSize]
                genes_inds = enhValueSortedArg[geneBin*geneBinSize:(geneBin+1)*geneBinSize]
                myExps = expArray[genes_inds] # >>> Expression
                expRatios[geneBin] = np.sum(myExps>0) / geneBinSize # >>> Expression
                #expMeansZ[geneBin] = np.mean(myExps[myExps>0])
                #expMeans[geneBin] = np.mean(myExps) 
                myEnhVals = enhancerValue[genes_inds]
                #enhMean[geneBin] = np.mean(enhSorted[geneBin*geneBinSize:(geneBin+1)*geneBinSize]) /(ev*2)
                enhMeans[geneBin] = np.mean(myEnhVals) # >>> Signal

                myEnhCovers = enhancerCover[genes_inds] # >>> Coverage
                enhCovers[geneBin] = np.mean(myEnhCovers) # >>> Coverage

                myEnhLowCovers = enhancerLowCover[genes_inds] # >>> Coverage Low
                enhLowCovers[geneBin] = np.mean(myEnhLowCovers) # >>> Coverage Low

                #enhNormMean[geneBin] = np.mean(enhNormSorted[geneBin*geneBinSize:(geneBin+1)*geneBinSize])
                #expMediansZ[geneBin] = np.median(myExps[myExps>0])
                geneBin +=1


            axs[0,e].plot((enhMeans[0:246]/((ev*2))))
            #axs[0,e].set_title('ratio of genes expressed')
            #axs[0,e].set_ylim((-.1,10))            
            axs[1,e].plot(enhCovers/(ev*2))
            #axs[1,e].set_title('mean enh label coverage')
            axs[1,e].set_ylim((-.1,1.1))
            axs[2,e].plot(expRatios)
            #axs[2,e].set_title('mean h3k4me1 value for enh regions')
            axs[2,e].set_ylim((-.1,1))
            axs[3,e].plot(enhLowCovers/(ev*2))
            axs[3,e].set_ylim((-.1,1.1))

            '''
            if e ==1:
                axs[0,e].set_title('Ratio of genes expressed')
                axs[1,e].set_title('mean enh and enhLow label coverage')
                axs[2,e].set_title('mean h3k4me1 value for enh an enhLow regions')
            '''
            
            if e ==1:
                axs[0,e].set_title('mean H3K4me1 value')
                axs[1,e].set_title('coverage of the enhLabel')
                axs[2,e].set_title('portion of expressed genes')
                axs[3,e].set_title('coverage of enhlowlabel')

            #axs[3].plot(expMeansZ)
            #axs[3].set_title('mean exp value')
            #axs[4].plot(expMediansZ)
            #axs[4].set_title('median exp value')

        #plt.tight_layout()
        #plt.show()

        #figFile = plotFolder + 'SegwayH3K4me1_plot_oneSample_%s_enhLowAndEnh_v02.pdf' %(accession)
        figFile = plotFolder + 'SegwayH3K4me1_plot_oneSample_%s_enh_v02.pdf' %(accession)
        print(figFile)
        plt.tight_layout()
        plt.savefig(figFile, bbox_inches='tight')
        plt.close('all')


        # >>>>>>>>>> Chrom evList
        fig, axs = plt.subplots(4,3, figsize=(8,8))
        for e, ev in enumerate(evList):
            myFilter = matFilters[(accession, ev)]

            expArray = np.asarray([expression[x] for x in geneIDList])[0:24700] # genes in the transcriptomic data that are in the gene map file
            
            side = 'right'
            inputFile = annotationFolder  + 'chrom_extendedGene_labelCover_%d_%s_hvals.pkl' %(ev,side)
            with open(inputFile, 'rb') as f:
                myMat_right = pickle.load(f)

            side = 'left'
            inputFile = annotationFolder  + 'chrom_extendedGene_labelCover_%d_%s_hvals.pkl' %(ev,side)
            with open(inputFile, 'rb') as f:
                myMat_left = pickle.load(f)

            myMat = myMat_left['histoneVals'] + myMat_right['histoneVals']# coverage of labels
            normMat = myMat_left['expMat'] + myMat_right['expMat']

            #enhancerValue = np.sum(myMat[:, [0,1,3,4, 5]], axis=1)
            enhancerValue = np.sum(myMat, axis=1)[0:24700] # sum value of H3k4me1
            enhancerLowValue = myMat[:, 5]
            enhancerCover = np.sum(normMat[:, [0,1,3,4]], axis=1)[0:24700]
            enhancerLowCover = normMat[0:24700, 5]

            enhancerValue = enhancerValue[myFilter]
            expArray = expArray[myFilter]
            enhancerCover = enhancerCover[myFilter]
            enhancerLowCover = enhancerLowCover[myFilter]

            '''
            enhCoverNorm = np.zeros(len(enhancerCover))
            for k in range(len(enhCoverNorm)):
                enhCoverNorm[k] = enhancerValue[k] / (enhancerCover[k]+1)
            '''
            #enhNormSortedArg = np.argsort(enhCoverNorm[0:24700])
            #enhNormSorted = np.sort(enhCoverNorm[0:24700])
            #enhSortedArg = np.argsort(enhancerCover[0:24700])
            #enhSorted = np.sort(enhancerCover[0:24700])
            enhValueSortedArg = np.argsort(enhancerValue[0:24700]) # >>> this is what we use for sorting
            
            geneBinSize = 100
            geneBin = 0
            geneBinCount = int(len(expArray)/geneBinSize)
            expRatios = np.zeros(geneBinCount) # ratio of genes expressed
            #expMeans =  np.zeros(geneBinCount) # mean expression value
            #expMeansZ =  np.zeros(geneBinCount) # mean expression for exp genes
            #expMediansZ = np.zeros(geneBinCount) # median expression for exp genes
            
            # we need signal, coverage, expression in the plot
            enhMeans = np.zeros(geneBinCount) # mean enhancer coverage
            enhCovers = np.zeros(geneBinCount) # mean enhancer coverage >>>> Coverage
            enhLowCovers = np.zeros(geneBinCount) # mean enhancer coverage >>>> Coverage Low            
            #enhNormMean = np.zeros(geneBinCount) # mean enhancer value for coverage
            while geneBin < geneBinCount:
                # genes_inds = enhNormSortedArg[geneBin*geneBinSize:(geneBin+1)*geneBinSize]
                genes_inds = enhValueSortedArg[geneBin*geneBinSize:(geneBin+1)*geneBinSize]
                myExps = expArray[genes_inds] # >>> Expression
                expRatios[geneBin] = np.sum(myExps>0) / geneBinSize
                #expMeansZ[geneBin] = np.mean(myExps[myExps>0])
                #expMeans[geneBin] = np.mean(myExps)
                myEnhVals = enhancerValue[genes_inds]
                #enhMean[geneBin] = np.mean(enhSorted[geneBin*geneBinSize:(geneBin+1)*geneBinSize]) /(ev*2)
                #enhNormMean[geneBin] = np.mean(enhNormSorted[geneBin*geneBinSize:(geneBin+1)*geneBinSize])
                #expMediansZ[geneBin] = np.median(myExps[myExps>0])
                enhMeans[geneBin] = np.mean(myEnhVals) # >>> Signal

                myEnhCovers = enhancerCover[genes_inds] # >>> Coverage
                enhCovers[geneBin] = np.mean(myEnhCovers) # >>> Coverage

                myEnhLowCovers = enhancerLowCover[genes_inds] # >>> Coverage Low
                enhLowCovers[geneBin] = np.mean(myEnhLowCovers) # >>> Coverage Low
                
                geneBin +=1


            axs[0,e].plot((enhMeans[0:246]/((ev*2))))
            #axs[0,e].set_title('ratio of genes expressed')
            #axs[0,e].set_ylim((-.1,10))            
            axs[1,e].plot(enhCovers/(ev*2))
            #axs[1,e].set_title('mean enh label coverage')
            axs[1,e].set_ylim((-.1,1.1))
            axs[2,e].plot(expRatios)
            #axs[2,e].set_title('mean h3k4me1 value for enh regions')
            axs[2,e].set_ylim((-.1,1))
            axs[3,e].plot(enhLowCovers/(ev*2))
            axs[3,e].set_ylim((-.1,1.1))


            '''
            if e ==1:
                axs[0,e].set_title('Ratio of genes expressed')
                axs[1,e].set_title('mean enh and enhLow label coverage')
                axs[2,e].set_title('mean h3k4me1 value for enh an enhLow regions')
            '''
            
            if e ==1:
                axs[0,e].set_title('mean H3K4me1 value')
                axs[1,e].set_title('coverage of the enhLabel')
                axs[2,e].set_title('portion of expressed genes')
                axs[3,e].set_title('coverage of enhlowlabel')


        figFile = plotFolder + 'ChromH3K4me1_plot_oneSample_%s_enh_v02.pdf' %(accession)
        print(figFile)
        plt.tight_layout()
        plt.savefig(figFile, bbox_inches='tight')
        plt.close('all')


##############################
# find the fair  example
##############################

indexList # list of available accessions
sib = np.mean(chromAUCs, axis=1) # aucs 
inIndsSortIndex = np.argsort(sib) # sort index for inInds, based on aucs

inInds = []
for a, accession in enumerate(accessionList):
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']
    if not(annotation['RNAseqFile'] == 'none') and not(annotation['chromFile'] == 'none'):
        inInds.append(a)

# sorted accessions:
inIndsSorted = [] # sorted list of indices based on AUCs
for b in (inIndsSortIndex):
    inIndsSorted.append(inInds[b])

middleInds = inIndsSorted[38:55]
for i in indexList: # we want to see where they are in the inIndsSorted
    if i in middleInds:
        print(i)
    #print(i)
    print(accessionList[i])
    print(inIndsSorted.index(i)) # where is it in the sorted list
        
# for the three coverage system, sort by h3k4 values, plot the expression level of the genes based on that sort
# plot the coverage of segway enhancer label
# plot the coverage of chrom enhancer label

########################################
# # getting some items for plot 
########################################

 accession
'ENCSR388IAJ'


########################################
# 5. intensity of H3K4me1 within a distance from the genes
########################################
# get the meanSig from 3.3.0 for one sample

# get the list of genes index in chr19
cgi = 0
chr19List = []
while cgi<24700:
    myGene = geneIDList[cgi]
    if geneList[geneIDList[cgi]].chrom == 'chr19':
        chr19List.append(cgi)
    cgi+=1

chr19ExpList = []
for i in chr19List:
    if expArray[i] > 0:
        chr19ExpList.append(i)

chr19ExpListZero = []
for i in chr19List:
    if expArray[i] == 0:
        chr19ExpListZero.append(i)

# now for each gene, get the count of regions with any foldchange, right around the promoter, and at the rest of it.
signalCount = np.zeros((len(chr19ExpListZero), 15))
for c, cgi in enumerate(chr19ExpListZero):
    geneStart = geneList[geneIDList[cgi]].start
    geneEnd = geneList[geneIDList[cgi]].end
    strand = geneList[geneIDList[cgi]].strand
    if strand == '+':
        promRegionS = int((geneStart - 2000)/100)
        promRegionE = int((geneStart + 300)/100)

    if strand == '-':
        promRegionS = int((geneEnd - 300)/100)
        promRegionE = int((geneEnd + 2000)/100)

    # promoter region
    signalCount[c, 0] = np.sum(meanSig[promRegionS:promRegionE] < 1)
    signalCount[c, 1] = np.sum(meanSig[promRegionS:promRegionE] < 2) - signalCount[c,0]
    signalCount[c, 2] = np.sum(meanSig[promRegionS:promRegionE] < 3) - np.sum(signalCount[c,0:2])
    signalCount[c, 3] = np.sum(meanSig[promRegionS:promRegionE] < 4) - np.sum(signalCount[c,0:3])    
    signalCount[c, 4] = (promRegionE - promRegionS) - np.sum(signalCount[c,0:4])

    if strand == '+':
        rightS = int((geneStart + 300)/100)
        rightE = int((geneStart + 300 + 2000)/100)

    if strand == '-':
        rightS = int((geneEnd + 2000)/100)
        rightE = int((geneEnd + 4000)/100)

    # right region
    signalCount[c, 5] = np.sum(meanSig[rightS:rightE] < 1)
    signalCount[c, 6] = np.sum(meanSig[rightS:rightE] < 2) - signalCount[c,5]
    signalCount[c, 7] = np.sum(meanSig[rightS:rightE] < 3) - np.sum(signalCount[c,5:7])
    signalCount[c, 8] = np.sum(meanSig[rightS:rightE] < 4) - np.sum(signalCount[c,5:8])    
    signalCount[c, 9] = (rightE - rightS) - np.sum(signalCount[c,5:9])


    if strand == '+':
        leftS = int((geneStart - 4000)/100)
        leftE = int((geneStart - 2000)/100)

    if strand == '-':
        leftS = int((geneEnd - 2300)/100)
        leftE = int((geneEnd - 300)/100)

    # left region
    signalCount[c, 10] = np.sum(meanSig[leftS:leftE] < 1)
    signalCount[c, 11] = np.sum(meanSig[leftS:leftE] < 2) - signalCount[c,10]
    signalCount[c, 12] = np.sum(meanSig[leftS:leftE] < 3) - np.sum(signalCount[c,10:12])
    signalCount[c, 13] = np.sum(meanSig[leftS:leftE] < 4) - np.sum(signalCount[c,10:13])    
    signalCount[c, 14] = (leftE - leftS) - np.sum(signalCount[c,10:14])

book = signalCount[:, [1,2, 3, 4, 6,7, 8, 9, 11,12, 13,14]]
book = signalCount.copy()
bookSumCoverage = book[:, [5,6,7,8,9]] + book[:,[10,11,12,13,14]]
bookProm = book[:,[0,1,2,3,4]]
plt.boxplot(book)
plt.show()

halva = signalCount[:, [1,2, 3, 4, 6,7, 8, 9, 11,12, 13,14]]
halva = signalCount.copy()
halvaSumCoverage = halva[:, [5,6,7,8,9]] + halva[:,[10,11,12,13,14]]
halvaProm = halva[:, [0,1,2,3,4]]
plt.boxplot(halva)
plt.show()

fig, axs = plt.subplots(1,2, figsize=(8,8))

axs[0].boxplot(book)
axs[1].boxplot(halva)
plt.show()

fig, axs = plt.subplots(2,2, figsize=(8,8))

axs[0,0].boxplot(bookSumCoverage)
axs[0,1].boxplot(halvaSumCoverage)
axs[1,0].boxplot(bookProm)
axs[1,1].boxplot(halvaProm)

plt.show()

        
########################################
# 6. Between the two samples, get the enhancers between the genes that are expressed vs genes that are not
########################################

# tissue or sample specific gene expression in chromosome 19: I want a selection of genes that are only highly expressed in a few samples, and I want to plot the activity of the genome for various regions around them.

# for chromosome 20, how many genes are expressed in each sample, and pick a gene that is only expressed in one sample and the enhancer area around it with the scatter plot. 

# >>>>>>>>>>>>>>>>>>>>>>>>>>
label_term_mapping = {}
mnemonics_file = annotationFolder + 'mnemonics_v04.txt'
with open(mnemonics_file, 'r') as mnemonics:
    for line in mnemonics:
        print(line)
        label = line.strip().split()[0]
        term = line.strip().split()[1]
        label_term_mapping[label] = term


# get the track file
assay_file = annotationFolder + 'trackname_assay.txt'
with open(assay_file, 'r') as assays:
    for line in assays:
        if line.split()[1] == 'H3K4me1':
            print(line)
            h3k4 = line.split()[0]
            print(accession)


histoneFile = '/Users/marjanfarahbod/Downloads/%s.bigWig' %(h3k4)
print(file)

annotationFolder = annotation['folder']
print(annotationFolder)

annFile = annotation['bedFile']
mnemFile = annotationFolder + 'mnemonics_v04.txt'
#geneList
#geneIDList =
ev = 2000
side = 'right'
EPLabelCover_Segway_extensionValue(annotationFolder, annFile, mnemFile, geneList, geneIDList, ev, side, histoneFile)

ev = 2000
side = 'left'
inputFile = annotationFolder  + 'extendedGene_labelCover_%d_%s_hvals.pkl' %(ev,side)
with open(inputFile, 'rb') as f:
    myMat_left = pickle.load(f)

side = 'right'
inputFile = annotationFolder  + 'extendedGene_labelCover_%d_%s_hvals.pkl' %(ev,side)
with open(inputFile, 'rb') as f:
    myMat_right = pickle.load(f)

#myMat = myMat_right['expMat'] + myMat_left['expMat']
myMat = myMat_left['expMat'] # coverage of labels
hvall = myMat_left['histoneVals'] # sum value of regions with labels

enhArray_hisl = np.zeros(normExpMat.shape[0]) # sum of histone fold change values for enhancer regions
lowEnhArray_hisl = np.zeros(normExpMat.shape[0]) # sum of histone fold change values for enhancer low regions
enhancerCover = np.zeros(normExpMat.shape[0]) # total coverage of enhancer regions
enhancerLowCover = np.zeros(normExpMat.shape[0]) # total coverage of enhancer regions
for label in range(len(label_term_mapping)):
    if label_term_mapping[str(label)] == 'Enhancer':
        enhArray_hisl += hvall[:, label]
        enhancerCover += myMat[:, label]
    if label_term_mapping[str(label)] == 'EnhancerLow':
        lowEnhArray_hisl += hvall[:, label]
        enhancerLowCover += myMat[:,label]

meanHistoneLeft = np.zeros(hvall.shape[0])
for i,value in enumerate(enhArray_hisl):
    meanHistoneLeft[i] = value / (enhancerCover[i] + 1) # average histone value for regions with enhancer labels

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#myMat = myMat_right['expMat'] + myMat_left['expMat']
myMat = myMat_right['expMat'] # coverage of labels
hvalR = myMat_right['histoneVals'] # sum value of regions with labels

enhArray_hisR = np.zeros(normExpMat.shape[0]) # sum of histone fold change values for enhancer regions
lowEnhArray_hisR = np.zeros(normExpMat.shape[0]) # sum of histone fold change values for enhancer low regions
enhancerCover_R = np.zeros(normExpMat.shape[0]) # total coverage of enhancer regions
enhancerLowCover_R = np.zeros(normExpMat.shape[0]) # total coverage of enhancer regions
for label in range(len(label_term_mapping)):
    if label_term_mapping[str(label)] == 'Enhancer':
        enhArray_hisR += hvalR[:, label]
        enhancerCover_R += myMat[:, label]
    if label_term_mapping[str(label)] == 'EnhancerLow':
        lowEnhArray_hisR += hvalR[:, label]
        enhancerLowCover_R += myMat[:,label]

meanHistoneRight = np.zeros(hvall.shape[0])
for i,value in enumerate(enhArray_hisR):
    meanHistoneRight[i] = value / (enhancerCover_R[i] + 1) # average histone value for regions with enhancer labels

# 0.2 get the expression data
# >>>>>>>>>>>>>>>>>>>>>>>>>>
if len(annotation['RNAseqFile'][0]) > 1:
    RNAFile = annotation['RNAseqFile'][0]
else:
    RNAFile = annotation['RNAseqFile']
    print(RNAFile)

expAccession = RNAFile[-15:-4]
expFile = annotationFolder +  'geneExp_dict_' + expAccession + '.pkl'

with open(expFile, 'rb') as pickledFile:
    expression = pickle.load(pickledFile) # load the expression file


sumExp = myMat.sum(axis = 1)
normExpMat = myMat/sumExp[:,None] # normalizing the coverage to get ratio of the coverage by the 

enhArray = np.zeros(normExpMat.shape[0])
lowEnhArray = np.zeros(normExpMat.shape[0])
for label in range(len(label_term_mapping)):
    if label_term_mapping[str(label)] == 'Enhancer':
        enhArray += normExpMat[:, label]
    if label_term_mapping[str(label)] == 'EnhancerLow':
        lowEnhArray += normExpMat[:, label]
        
expArray = np.asarray([expression[x] for x in geneIDList])[0:24700] # genes in the transcriptomic data that are in the gene map file

book = np.argsort(expArray)

enhSortedArg = np.argsort(enhArray)[::-1]
enhSortedArg_histone = np.argsort(meanHistoneLeft)[::-1]
enhSortedArg_histoneR = np.argsort(meanHistoneRight)[::-1]

sumHis = meanHistoneLeft + meanHistoneRight
enhSortedArg_sumHis = np.argsort(sumHis)[::-1]

xwalk = 0
ywalk = 0
area = 0
aucc = np.zeros((len(expArray), 2))
for j in range(len(expArray)):
    if expArray[enhSortedArg_sumHis[j]] == 0:
        xwalk +=1
        area+= ywalk
    else:
        ywalk +=1

    aucc[j, 0] = xwalk
    aucc[j, 1] = ywalk

auc = area / (xwalk*ywalk)
print(auc)
plt.plot(aucc[:, 0], aucc[:,1])
pplot(enhSortedArray)
plt.show()

# get the 
########################################

for i in badInds[0]:
    print(rnaAccessionList[i])

