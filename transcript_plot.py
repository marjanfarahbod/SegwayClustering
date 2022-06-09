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

dataSubFolder = 'testBatch105/fromAPI/'

inputFile = dataFolder + dataSubFolder + 'metaInfo.pkl'
with open(inputFile, 'rb') as f:
    annMeta = pickle.load(f)

segwayLabels = ['Quiescent', 'ConstitutiveHet', 'FacultativeHet', 'Transcribed', 'Promoter', 'Enhancer', 'RegPermissive', 'Bivalent', 'LowConfidence']
    
# 
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

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
    
    plt.show()


    # todo: get the count of genes as expressed, low expressed and high expressed

     


