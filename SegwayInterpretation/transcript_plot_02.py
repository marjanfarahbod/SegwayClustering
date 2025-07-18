# this is a modified version of transcript_plot_01.py. Instead of getting the labels from one label file, I get them from mnemonics file. 

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
runFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/runs/'

dataSubFolder = 'testBatch105/fromAPI/'

inputFile = dataFolder + dataSubFolder + 'metaInfo.pkl'
with open(inputFile, 'rb') as f:
    annMeta = pickle.load(f)

inputFileName = 'all_annInfo_list.pkl'
inputFile = dataFolder + dataSubFolder + inputFileName
with open(inputFile, 'rb') as f:
    ann_info_list = pickle.load(f)

# get the mapping from accession to index
accession_index_map = {}
for ann in ann_info_list:
   accession_index_map[ ann['accession'] ] = ann['index']

#classifier_input_file = runFolder + "model_296exp_reg0.067_auc0.77on.32test_classifier_labels.pickle"
#with open(classifier_input_file, "rb") as f:
    #classifier_labels = pickle.load(f)

#book = classifier_labels['cluster_list']
#interpretation_terms_set = set(classifier_labels['interpretation_term'])

#segwayLabels = ['Quiescent', 'ConstitutiveHet', 'FacultativeHet', 'Transcribed', 'Promoter', 'Enhancer', 'RegPermissive', 'Bivalent', 'LowConfidence']
segwayLabels = ['Enhancer_low', 'Enhancer', 'Promoter_flanking', 'Promoter', 'Transcribed', 'CTCF', 'K9K36', 'Bivalent', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent', 'Unclassified']

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

    index = accession_index_map[annAccession]
    
    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
    print(sampleFolderAdd)

    # get the label from label from mnemonics: a dictionary from labels to terms
    label_term_mapping = {}
    mnemonics_file = sampleFolderAdd + 'mnemonics_v02.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term
    
    expFile = 'defaultExp_5kg_expSummary.pkl'
    inputFile = sampleFolderAdd + expFile

    try:
        with open(inputFile, 'rb') as f:
            expMat = pickle.load(f)
    except FileNotFoundError:
        print('no transcript data for this file')
        continue


    # load segway cluster info
    segwayBedAccession = annMeta[annAccession]['bedFile'][0:-4]
    annSummaryFile = sampleFolderAdd + segwayBedAccession + '_annotationSummary.pkl'

    with open(annSummaryFile, 'rb') as f:
        summaryAnnotation = pickle.load(f)

    # sorting the cluster labels for segway (this is the order used for the cluster matrix, see QC transcript code file)
    clusters = summaryAnnotation['clusters']
    clusterList = list(clusters.keys())

    # loading the cluster list for the transcript data based on how it was saved in QC transcript code file
    order_file = sampleFolderAdd + 'defaultExp_5kg_expSummary_clusterOrder.pkl'
    with open(order_file, 'rb') as f:
         cluster_ordered_list = pickle.load(f)

    # get the terms for the clusters
    updatedClusterList = []
    for clusterID in cluster_ordered_list:
        term = label_term_mapping[clusterID]
        label_term = clusterID + '_' + term
        updatedClusterList.append(label_term)
    
    updatedClusterList_reordered = []
    #cluster_order = [] # for saving the cluster order, I did it once here and it is saved, so it is commented
    for label in segwayLabels:
        #print(label)
        for item in updatedClusterList:
            #print(item)
            if item.split('_')[1] == label:
                updatedClusterList_reordered.append(item)
                #cluster_order.append(item.split('_')[0]) # for saving the cluster order, I did it once here and it is saved, so it is commented

    print(updatedClusterList_reordered)

    # we do three heatmap plots in one frame
    clusterCount = len(clusterList)
    fractions = np.zeros([clusterCount])
    total_bp = 0

    clusterID_bp_count = {}
    for cluster in clusters:
        clusterID = cluster.split('_')[0]
        clusterID_bp_count[clusterID] = clusters[cluster].bp_count

    for i, cluster in enumerate(updatedClusterList):
        clusterID = cluster.split('_')[0]
        print(clusterID_bp_count[clusterID])
        fractions[i] = clusterID_bp_count[clusterID]
        total_bp += clusterID_bp_count[clusterID]

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

    #h1 = pd.DataFrame(np.log10(thisMat, out=np.zeros_like(thisMat), where=(thisMat!=0)), index = sortedClusterList)
    h1 = pd.DataFrame(np.log10(thisMat, out=np.zeros_like(thisMat), where=(thisMat!=0)), index = updatedClusterList)
    h1 = h1.reindex(updatedClusterList_reordered)
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

    h2 = pd.DataFrame(np.log10(thisMat, out=np.zeros_like(thisMat), where=(thisMat!=0)), index = updatedClusterList)
    h2 = h2.reindex(updatedClusterList_reordered)
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

    h3 = pd.DataFrame(np.log10(thisMat, out=np.zeros_like(thisMat), where=(thisMat!=0)), index = updatedClusterList)
    h3 = h3.reindex(updatedClusterList_reordered)
    cmap = sns.diverging_palette(240, 10, s=80, l=30, as_cmap=True)
    # g1 = sns.heatmap(h1, center=0, cmap=cmap)
    g3 = sns.heatmap(h3, center=0, cmap=cmap, vmin=-2, vmax=2,  ax=axs[2], yticklabels=False, xticklabels=xticklabels)
    g3.set_xticks(xticks)
    g3.hlines(np.linspace(1,clusterCount, clusterCount), *g1.get_xlim(), color = 'grey')
    g3.set_title('high expression \n %d genes' %(highExpCount), fontsize = 10)
    #plt.tight_layout()

    fig.suptitle(annAccession + ' - ' + tissue_info[annAccession][0] + ' - ' + tissue_info[annAccession][1], fontsize = 12)

    # plt.show()
    plotFolder_add = plotFolder + annAccession + '/'
    figFile = plotFolder_add + 'expression_transcript_03.pdf'
    print(figFile)
    plt.savefig(figFile, bbox_inches='tight')
    plt.close('all')
