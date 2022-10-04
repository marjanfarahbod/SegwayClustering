# plots for the transcript comparison data

import linecache
import pickle
import re
import numpy as np
import pandas as pd
import subprocess as sp
import seaborn as sns
import matplotlib.pyplot as plt
from QC_transcriptionComparison_util import Gene, Exon, Annotation, AnnotationClass


# open the folder and look up the stuff
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/testBatch105/'

dataSubFolder = 'testBatch105/fromAPI/'

inputFile = dataFolder + dataSubFolder + 'metaInfo.pkl'
with open(inputFile, 'rb') as f:
    annMeta = pickle.load(f)

segwayLabels = ['Quiescent', 'ConstitutiveHet', 'FacultativeHet', 'Transcribed', 'Promoter', 'Enhancer', 'RegPermissive', 'Bivalent', 'LowConfidence']
    
# # with classes
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

inputFile = dataFolder + dataSubFolder + 'biosample_tissue_info.pkl'
with open(inputFile, 'rb') as f:
    tissue_info = pickle.load( f)

# GTF data structure
fileName = dataFolder + '/geneLists.pkl'
with open(fileName, 'rb') as pickleFile:
    geneListsAndIDs = pickle.load(pickleFile)
    
geneList = geneListsAndIDs[0]
geneIDList = geneListsAndIDs[1]
del geneListsAndIDs

annAccessionList = list(annMeta.keys())

for annAccession in annAccessionList:
    
    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
    print(sampleFolderAdd)
    
    expFile = 'defaultExp_5kg_expSummary.pkl'
    inputFile = sampleFolderAdd + expFile

    try:
        with open(inputFile, 'rb') as f:
            expMat = pickle.load(f)
    except FileNotFoundError:
        print('no transcript data for this file')
        continue

    # plot the expmat

    # get the labels for clusterMats

    # load segway cluster info
    segwayBedAccession = annMeta[annAccession]['bedFile'][0:-4]
    annSummaryFile = sampleFolderAdd + segwayBedAccession + '_annotationSummary.pkl'

    with open(annSummaryFile, 'rb') as f:
        summaryAnnotation = pickle.load(f)

    # sorting the cluster labels for segway (this is the order used for the cluster matrix, see QC transcript code file)
    clusters = summaryAnnotation['clusters']
    clusterList = list(clusters.keys())
    sortedClusterList = []
    for label in segwayLabels:
        #print(label)
        for item in clusterList:
            #print(item)
            if item.split('_')[1] == label:
                sortedClusterList.append(item)

    # we do three heatmap plots in one frame
    clusterCount = len(clusterList)
    fractions = np.zeros([clusterCount])
    total_bp = 0

    for i, cluster in enumerate(sortedClusterList):
        print(clusters[cluster].bp_count)
        fractions[i] = clusters[cluster].bp_count
        total_bp += clusters[cluster].bp_count

    fractions = fractions / total_bp

    # get the header: count of genes expressed zero, low and med/high
    expFileName = annMeta[annAccession]['expressionFile']
    expAccession = re.split('_|\.', expFileName)[2]
    inputFile = dataFolder + dataSubFolder + annAccession + '/geneExp_dict_' + expAccession + '.pkl'
    print(inputFile)
    with open(inputFile, 'rb') as pickledFile:
        expression = pickle.load(pickledFile)

    # get the list of gene IDs
    
    allExp = [expression[x] for x in geneIDList]
    zeroExpCount = sum(i == 0 for i in allExp)
    lowMedExpCount = sum(i > 0 and i < 1 for i in allExp)
    highExpCount = sum(i > 1 for i in allExp)


    fig, axs = plt.subplots(1, 3, figsize =(12, 6), gridspec_kw={'width_ratios':[1,1,1.2]})
    xticks = [30, 130]
    xticklabels = ['TSS', 'end']
    
    # heatmap1 : the not expressed
    thisMat = expMat['clusterMats'][0]

    # make it the ratio
    thisMat = thisMat / np.sum(thisMat, axis=0)[np.newaxis,:]
    # versus expected
    thisMat = thisMat / fractions[:, None]

    h1 = pd.DataFrame(np.log10(thisMat, out=np.zeros_like(thisMat), where=(thisMat!=0)), index = sortedClusterList)
    cmap = sns.diverging_palette(240, 10, s=80, l=30, as_cmap=True)
    # g1 = sns.heatmap(h1, center=0, cmap=cmap)
    g1 = sns.heatmap(h1, center=0, cmap=cmap, vmin=-2, vmax=2, ax=axs[0], cbar=False, xticklabels=xticklabels)
    g1.set_xticks(xticks)
    #sns.despine(fig=None, ax=None, top=False, right=False, left=False, bottom=False, offset=None, trim=False)
    g1.hlines(np.linspace(1,clusterCount, clusterCount), *g1.get_xlim(), color = 'grey')
    g1.set_title('zero expression\n %d genes' %(zeroExpCount), fontsize = 10)

    #plt.show()

    # heatmap 2: the low expressed
    thisMat = expMat['clusterMats'][1]

    # make it the ratio
    thisMat = thisMat / np.sum(thisMat, axis=0)[np.newaxis,:]
    # versus expected
    thisMat = thisMat / fractions[:, None]

    h2 = pd.DataFrame(np.log10(thisMat, out=np.zeros_like(thisMat), where=(thisMat!=0)), index = sortedClusterList)
    cmap = sns.diverging_palette(240, 10, s=80, l=30, as_cmap=True)
    # g1 = sns.heatmap(h1, center=0, cmap=cmap)
    g2 = sns.heatmap(h2, center=0, cmap=cmap, vmin=-2, vmax=2,  ax=axs[1], cbar=False, yticklabels=False, xticklabels=xticklabels)
    g2.set_xticks(xticks)
    g2.hlines(np.linspace(1,clusterCount, clusterCount), *g1.get_xlim(), color = 'grey')
    g2.set_title('low or medium expression\n %d genes' %(lowMedExpCount), fontsize = 10)

    # heatmap 3: the high expressed
    thisMat = expMat['clusterMats'][2]

    # make it the ratio
    thisMat = thisMat / np.sum(thisMat, axis=0)[np.newaxis,:]
    # versus expected
    thisMat = thisMat / fractions[:, None]

    h3 = pd.DataFrame(np.log10(thisMat, out=np.zeros_like(thisMat), where=(thisMat!=0)), index = sortedClusterList)
    cmap = sns.diverging_palette(240, 10, s=80, l=30, as_cmap=True)
    # g1 = sns.heatmap(h1, center=0, cmap=cmap)
    g3 = sns.heatmap(h3, center=0, cmap=cmap, vmin=-2, vmax=2,  ax=axs[2], yticklabels=False, xticklabels=xticklabels)
    g3.set_xticks(xticks)
    g3.hlines(np.linspace(1,clusterCount, clusterCount), *g1.get_xlim(), color = 'grey')
    g3.set_title('high expression \n %d genes' %(highExpCount), fontsize = 10)

    fig.suptitle(annAccession + ' - ' + tissue_info[annAccession][0] + ' - ' + tissue_info[annAccession][1], fontsize = 12)

    # plt.show()
    plotFolder_add = plotFolder + annAccession + '/'
    figFile = plotFolder_add + 'expression_transcript.pdf'
    print(figFile)
    plt.savefig(figFile)
    plt.close('all')

# with clusters
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# want to do the plot for exp summary files.
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/Mehdi_testRun/tomarjan/'
dataSubFolder = 'segway/'
sampleSubList = ['GM12878_rep1', 'K562_rep1_rs5', 'MCF-7_rep1', 'MCF-7_rep1_rs5', 'GM12878_rep1_rs7', 'K562_rep2', 'GM12878_rep2', 'MCF-7_rep2', 'MCF-7_rep1_rs7', 'K562_rep1_rs7', 'K562_rep1', 'GM12878_rep1_rs5']

# GTF data structure
fileName = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/geneLists.pkl'
with open(fileName, 'rb') as pickleFile:
    geneListsAndIDs = pickle.load(pickleFile)
    
geneList = geneListsAndIDs[0]
geneIDList = geneListsAndIDs[1]
del geneListsAndIDs

for item in sampleSubList:

    # get the list of expSummary files:
    expSum_folder = dataFolder + dataSubFolder + item + '/'
    expFileList = glob.glob(expSum_folder + '*_expSummary.pkl')

    # load the annotation summary file
    inputFile = expSum_folder + 'segway_annotationSummary.pkl'
    with open(inputFile, 'rb') as pickledFile:
        summaryAnnotation = pickle.load(pickledFile)
        
    clusters = summaryAnnotation # we only have clusters for these samples, no need to sort the clusters
    clusterList = list(clusters.keys())
    clusterCount = len(clusterList)

    clusterListIndMap = {}
    for i in range(len(clusterList)):
        clusterListIndMap[clusterList[i]] = i


    for expFile in expFileList:

        # get the accession of the expression file 
        expFileAccess = expFile.split('/')[-1].split('_')[0]
        # get the expression vector 
        RNAseqFolder = dataFolder + 'RNA_seq/' + item.split('_')[0] + '/'
        RNAseqFileList = glob.glob(expFolder + 'geneExp_dict_*.pkl')
        for file in RNAseqFileList:
            if expFileAccess in file:
                RNAseqFile = file

        with open(RNAseqFile, 'rb') as inputFile:
            expression = pickle.load(inputFile)

        allExp = [expression[x] for x in geneIDList]
        zeroExpCount = sum(i == 0 for i in allExp)
        lowMedExpCount = sum(i > 0 and i < 1 for i in allExp)
        highExpCount = sum(i > 1 for i in allExp)


        print(expFile)
        with open(expFile, 'rb') as pickledFile:
            expMat = pickle.load(pickledFile)
            
        # we do three heatmap plots in one frame
        clusterCount = len(clusterList)
        fractions = np.zeros([clusterCount])
        total_bp = 0

        for i, cluster in enumerate(clusterList):
            print(clusters[cluster].bp_count)
            fractions[i] = clusters[cluster].bp_count
            total_bp += clusters[cluster].bp_count

        fractions = fractions / total_bp
    
    
        fig, axs = plt.subplots(1, 3, figsize =(12, 6), gridspec_kw={'width_ratios':[1,1,1.2]})
        xticks = [30, 130]
        xticklabels = ['TSS', 'end']
    
        # heatmap1 : the not expressed
        thisMat = expMat[0]
        print(thisMat[12][25])

        # make it the ratio
        thisMat = thisMat / np.sum(thisMat, axis=0)[np.newaxis,:]
        # versus expected
        thisMat = thisMat / fractions[:, None]

        h1 = pd.DataFrame(np.log10(thisMat, out=np.zeros_like(thisMat), where=(thisMat!=0)), index = clusterList)
        cmap = sns.diverging_palette(240, 10, s=80, l=30, as_cmap=True)
        # g1 = sns.heatmap(h1, center=0, cmap=cmap)
        g1 = sns.heatmap(h1, center=0, cmap=cmap, vmin=-2, vmax=2, ax=axs[0], cbar=False, xticklabels=xticklabels)
        g1.set_xticks(xticks)
        #sns.despine(fig=None, ax=None, top=False, right=False, left=False, bottom=False, offset=None, trim=False)
        g1.hlines(np.linspace(1,clusterCount, clusterCount), *g1.get_xlim(), color = 'grey')
        g1.set_title('zero expression\n %d genes' %(zeroExpCount), fontsize = 10)

        #plt.show()

        # heatmap 2: the low expressed
        thisMat = expMat[1]

        # make it the ratio
        thisMat = thisMat / np.sum(thisMat, axis=0)[np.newaxis,:]
        # versus expected
        thisMat = thisMat / fractions[:, None]

        h2 = pd.DataFrame(np.log10(thisMat, out=np.zeros_like(thisMat), where=(thisMat!=0)), index = clusterList)
        cmap = sns.diverging_palette(240, 10, s=80, l=30, as_cmap=True)
        # g1 = sns.heatmap(h1, center=0, cmap=cmap)
        g2 = sns.heatmap(h2, center=0, cmap=cmap, vmin=-2, vmax=2,  ax=axs[1], cbar=False, yticklabels=False, xticklabels=xticklabels)
        g2.set_xticks(xticks)
        g2.hlines(np.linspace(1,clusterCount, clusterCount), *g1.get_xlim(), color = 'grey')
        g2.set_title('low or medium expression\n %d genes' %(lowMedExpCount), fontsize = 10)

        # heatmap 3: the high expressed
        thisMat = expMat[2]

        # make it the ratio
        thisMat = thisMat / np.sum(thisMat, axis=0)[np.newaxis,:]
        # versus expected
        thisMat = thisMat / fractions[:, None]

        h3 = pd.DataFrame(np.log10(thisMat, out=np.zeros_like(thisMat), where=(thisMat!=0)), index = sortedClusterList)
        cmap = sns.diverging_palette(240, 10, s=80, l=30, as_cmap=True)
        # g1 = sns.heatmap(h1, center=0, cmap=cmap)
        g3 = sns.heatmap(h3, center=0, cmap=cmap, vmin=-2, vmax=2,  ax=axs[2], yticklabels=False, xticklabels=xticklabels)
        g3.set_xticks(xticks)
        g3.hlines(np.linspace(1,clusterCount, clusterCount), *g1.get_xlim(), color = 'grey')
        g3.set_title('high expression \n %d genes' %(highExpCount), fontsize = 10)

        fig.suptitle(item + ' - expression accession: ' + expFileAccess, fontsize = 12)

        #plt.show()
        figFile = expSum_folder + expFileAccess + '_expression_transcript.pdf'
        print(figFile)
        plt.savefig(figFile)
        plt.close('all')




     


