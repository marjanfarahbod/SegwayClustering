# how much do they overlap (union intersect) for all the labels.
# how much do they overlap with the detected enhancers for the cell type - each.
# distribution of enhancer labels for the two samples. Samples investigated separately. (I can change the model if I like it)



########################################
# Phantom5 enhancer analyses is pending until I find stuff
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

import glob

# General data folder
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'

dataSubFolder = 'testBatch105/fromAPI/'

# GTF data structure
fileName = dataFolder + '/geneLists.pkl'
with open(fileName, 'rb') as pickleFile:
    geneListsAndIDs = pickle.load(pickleFile)
    
geneList = geneListsAndIDs[0]
geneIDList = geneListsAndIDs[1]
del geneListsAndIDs

inputFile = dataFolder + dataSubFolder + 'metaInfo.pkl'
with open(inputFile, 'rb') as f:
    annMeta = pickle.load(f)

annAccessionList = list(annMeta.keys())

inputFile = dataFolder + dataSubFolder + 'chmmFile_list_dict.pkl'
with open(inputFile, 'rb') as f:
    chmmFile_dict = pickle.load(f)

book =  'EnhA1 EnhA2 EnhBiv EnhG1 EnhG2 EnhWk Het Quies ReprPC ReprPCWk TssA TssBiv TssFlnk TssFlnkD TssFlnkU Tx TxWk ZNF/Rpts'
chromLabels = book.split()

# chromhmm segway comparison
########################################

segwayLabels = ['EnhancerLow', 'Enhancer', 'PromoterFlanking', 'Promoter', 'Transcribed', 'CTCF', 'K9K36', 'Bivalent', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent']

for annAccession in annAccessionList:

    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
    # get the name of the segway bed file from annMeta
    originalBedFile = annMeta[annAccession]['bedFile']
    bedAccession = originalBedFile.split('.')[0]


    chmmFile = chmmFile_dict[annAccession]
    chmmAccession = re.search('ENCFF.*\.', chmmFile)[0][:-1]
    
    inputFileName = 'overlap_segway_%s_chmm_%s.pkl' %(bedAccession, chmmAccession)
    inputFile = dataFolder + dataSubFolder + annAccession + '/' + inputFileName
    with open(inputFile, 'rb') as f:
        overlap = pickle.load(f)

    overlap_mat = overlap.to_numpy()

    # total bp count that is covered by both annotations
    total_bp = np.sum(overlap_mat)

    # fraction of base pairs in each label to the total bp count
    chmm_labelFraction = np.sum(overlap_mat, axis = 0) / total_bp
    segway_labelFraction = np.sum(overlap_mat, axis = 1) / total_bp

    # get the observed versus the expected fraction
    expectedFraction_mat = np.zeros(overlap_mat.shape)
    for i in range(len(segway_labelFraction)):
        for j in range(len(chmm_labelFraction)):
            expectedFraction_mat[i][j] = segway_labelFraction[i] * chmm_labelFraction[j]

    observedFraction_mat = np.zeros(overlap_mat.shape)
    for i in range(len(segway_labelFraction)):
        for j in range(len(chmm_labelFraction)):
            observedFraction_mat[i][j] = overlap_mat[i][j] / total_bp

    obs_exp = np.divide(observedFraction_mat, expectedFraction_mat)

    # get the normalized value
    overlap_mat_colNorm = overlap_mat / np.sum(overlap_mat, axis=0)[np.newaxis, :] # chmm
    overlap_mat_rowNorm = overlap_mat / np.sum(overlap_mat, axis=1)[:, np.newaxis] # segway

    
    # load the updated mnemonics
    label_term_mapping = {}
    mnemonics_file = sampleFolderAdd + 'mnemonics_v04.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term

    # change segway_cluster_list_updated based on the mnemonics
    segway_cluster_list = []
    for cluster in segway_cluster_list_old:
        label = cluster.split('_')[0]
        term = label_term_mapping[label]
        segway_cluster_list.append(label + '_' + term)

    # reorder the cluster list
    segway_cluster_list_reordered = []
    for label in segwayLabels:
        for cluster in segway_cluster_list:
            if cluster.split('_')[1] == label:
                segway_cluster_list_reordered.append(cluster)

                    # add the ratio to the axis labels
    segwayAxis_list = []
    for i, cluster in enumerate(segway_cluster_list):
        segwayAxis_list.append(cluster + '_' + str(round(segway_labelFraction[i], 4)))


    segwayAxis_list_reordered = []
    for label in segwayLabels:
        for axis in segwayAxis_list:
            if axis.split('_')[1] == label:
                segwayAxis_list_reordered.append(axis)

    chmm_class_list = list(overlap.columns.values)

    chmmAxis_list = []
    for i, chmmclass in enumerate(chmm_class_list):
        chmmAxis_list.append(chmmclass + '_' + str(round(chmm_labelFraction[i], 4)))

    fig, axs = plt.subplots(1, 3, figsize=(16,4.2))
    
    obs_exp_log = np.log10(obs_exp, out=np.zeros_like(obs_exp), where=(obs_exp!=0))
    obs_exp_log = np.where(obs_exp_log < 0, 0, obs_exp_log)

    h1 = pd.DataFrame(obs_exp_log, index = segwayAxis_list, columns = chmmAxis_list)
    h1 = h1.reindex(segwayAxis_list_reordered)
    #cmap = sns.diverging_palette(240, 10, s=100, l=30, as_cmap=True)
    cmap = sns.diverging_palette(240, 10, s=80, l=30, as_cmap=True)
    g1 = sns.heatmap(h1, center = 0,cmap=cmap, ax=axs[0])
    #g1.set_ylabel('this')
    #g1.set_xlabel('that')
    g1.set_title('ratio - log (observed to expected)')
    #g1.set_xticklabels(g1.get_xticklabels(), rotation=45)
    #plt.ylabel('segway')
    #plt.ylabel('chmm')
    #plt.subplots_adjust(left=.075, right=.95, top=.9, bottom=.25)
    sns.set(font_scale=.8)
    plt.tight_layout()
    #plt.title('ratio - observed vs expected')
    #plt.show()


    h1 = pd.DataFrame(overlap_mat_colNorm, index = segwayAxis_list, columns = chmmAxis_list)
    h1 = h1.reindex(segwayAxis_list_reordered)
    cmap = sns.diverging_palette(240, 10, s=80, l=30, as_cmap=True)
    g1 = sns.heatmap(h1, center = 0,cmap=cmap, ax=axs[1])
    #g1.set_ylabel('this')
    #g1.set_xlabel('that')
    g1.set_title('ratio of bp overlap - chmm')
    #plt.ylabel('segway')
    #plt.ylabel('chmm')
    #plt.subplots_adjust(left=.075, right=.95, top=.9, bottom=.25)
    plt.tight_layout()
    
    h2 = pd.DataFrame(overlap_mat_rowNorm, index = segwayAxis_list, columns = chmmAxis_list)
    h2 = h2.reindex(segwayAxis_list_reordered)
    cmap = sns.diverging_palette(240, 10, s=80, l=30, as_cmap=True)
    g2 = sns.heatmap(h2, center = 0,cmap=cmap, ax=axs[2])
    #g1.set_ylabel('this')
    #g1.set_xlabel('that')
    g2.set_title('ratio of bp overlap - Segway')
    #plt.ylabel('segway')
    #plt.ylabel('chmm')
    #plt.subplots_adjust(left=.075, right=.95, top=.9, bottom=.25)
    plt.tight_layout()

    plt.show()

# I am looking for a measurement between the samples that I can compare
# which percentage of the enhancer labels are covered with log10 = 1 overlap
# it is important to note that Segway identifies far more bs as enhancers as chromhmm
########################################

# what fraction of base-pairs from chmm enhancer regions were covered by Segway Enhacer/enhancerLow with more than log ratio .3

segwayLabels = ['EnhancerLow', 'Enhancer', 'PromoterFlanking', 'Promoter', 'Transcribed', 'CTCF', 'K9K36', 'Bivalent', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent']

# ITSELF
# bestmatch
bestMatch = {}
bestMatchValue = {}
for annAccession in annAccessionList[0:10]:

    print(annAccession)

    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
    # get the name of the segway bed file from annMeta
    originalBedFile = annMeta[annAccession]['bedFile']
    bedAccession = originalBedFile.split('.')[0]

    chmmFile = chmmFile_dict[annAccession]
    if chmmFile == 'none':
        continue
    chmmAccession = re.search('ENCFF.*\.', chmmFile)[0][:-1]
    inputFileName = 'overlap_segway_%s_chmm_%s.pkl' %(bedAccession, chmmAccession)
    inputFile = dataFolder + dataSubFolder + annAccession + '/' + inputFileName
    with open(inputFile, 'rb') as f:
        overlap = pickle.load(f)

    overlap_mat = overlap.to_numpy()

    # total bp count that is covered by both annotations
    total_bp = np.sum(overlap_mat)

    # fraction of base pairs in each label to the total bp count
    chmm_labelFraction = np.sum(overlap_mat, axis = 0) / total_bp
    segway_labelFraction = np.sum(overlap_mat, axis = 1) / total_bp

    # get the observed versus the expected fraction
    expectedFraction_mat = np.zeros(overlap_mat.shape)
    for i in range(len(segway_labelFraction)):
        for j in range(len(chmm_labelFraction)):
            expectedFraction_mat[i][j] = segway_labelFraction[i] * chmm_labelFraction[j]

    observedFraction_mat = np.zeros(overlap_mat.shape)
    for i in range(len(segway_labelFraction)):
        for j in range(len(chmm_labelFraction)):
            observedFraction_mat[i][j] = overlap_mat[i][j] / total_bp

    obs_exp = np.divide(observedFraction_mat, expectedFraction_mat)

    # get the normalized value
    overlap_mat_colNorm = overlap_mat / np.sum(overlap_mat, axis=0)[np.newaxis, :] # chmm
    overlap_mat_rowNorm = overlap_mat / np.sum(overlap_mat, axis=1)[:, np.newaxis] # segway

    segway_cluster_list_old = list(overlap.index.values)
    
    # load the updated mnemonics
    label_term_mapping = {}
    mnemonics_file = sampleFolderAdd + 'mnemonics_v04.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term

    # change segway_cluster_list_updated based on the mnemonics
    segway_cluster_list = []
    for cluster in segway_cluster_list_old:
        label = cluster.split('_')[0]
        term = label_term_mapping[label]
        segway_cluster_list.append(label + '_' + term)

    # reorder the cluster list
    segway_cluster_list_reordered = []
    for label in segwayLabels:
        for cluster in segway_cluster_list:
            if cluster.split('_')[1] == label:
                segway_cluster_list_reordered.append(cluster)

    # add the ratio to the axis labels
    segwayAxis_list = []
    segEnh_inds = []
    for i, cluster in enumerate(segway_cluster_list):
        segwayAxis_list.append(cluster + '_' + str(round(segway_labelFraction[i], 4)))
        if cluster.split('_')[1] == 'Enhancer':
            segEnh_inds.append(i)

    segwayAxis_list_reordered = []
    for label in segwayLabels:
        for axis in segwayAxis_list:
            if axis.split('_')[1] == label:
                segwayAxis_list_reordered.append(axis)

    chmm_class_list = list(overlap.columns.values)

    chmmAxis_list = []
    for i, chmmclass in enumerate(chmm_class_list):
        chmmAxis_list.append(chmmclass + '_' + str(round(chmm_labelFraction[i], 4)))

    obs_exp_log = np.log10(obs_exp, out=np.zeros_like(obs_exp), where=(obs_exp!=0))
    obs_exp_log = np.where(obs_exp_log < 0, 0, obs_exp_log)

    sumCovVal = np.zeros(len(segwayAxis_list))
    for i in range(len(segwayAxis_list)):
        sumCovVal[i] = np.sum(overlap_mat_colNorm[i,0:6]* obs_exp_log[i, 0:6])

    book = np.argsort(sumCovVal)
    for i in book:
        print('%s  %.4f' %(segwayAxis_list[i], sumCovVal[i]))

    bestMatch[annAccession] = segwayAxis_list[book[-1]]
    bestMatchValue[annAccession] = sumCovVal[book[-1]]
    
    # find the rows with label Enhancer (just that)
    # get the coverage ratio from the cells with log > .3

# OTHERS
########################################


segwayLabels = ['EnhancerLow', 'Enhancer', 'PromoterFlanking', 'Promoter', 'Transcribed', 'CTCF', 'K9K36', 'Bivalent', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent']

allMatch = {}
#selfMatch = {}
#selfMatchValue = {}
selfMatch = {}
for annAccession in annAccessionList[0:30]:

    print(annAccession)

    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')

    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
    # get the name of the segway bed file from annMeta
    originalBedFile = annMeta[annAccession]['bedFile']
    bedAccession = originalBedFile.split('.')[0]

    chmmFile = chmmFile_dict[annAccession]
    if chmmFile == 'none':
        continue
    chmmAccession = re.search('ENCFF.*\.', chmmFile)[0][:-1]
    
    inputFileName = 'overlap_segway_%s_chmm_%s.pkl' %(bedAccession, chmmAccession)
    inputFile = dataFolder + dataSubFolder + annAccession + '/' + inputFileName
    with open(inputFile, 'rb') as f:
        overlap = pickle.load(f)

    # do it for the self file
    overlap_mat = overlap.to_numpy()

    print(inputFileName)

    # total bp count that is covered by both annotations
    total_bp = np.sum(overlap_mat)

    # fraction of base pairs in each label to the total bp count
    chmm_labelFraction = np.sum(overlap_mat, axis = 0) / total_bp
    segway_labelFraction = np.sum(overlap_mat, axis = 1) / total_bp

    # get the observed versus the expected fraction
    expectedFraction_mat = np.zeros(overlap_mat.shape)
    for i in range(len(segway_labelFraction)):
        for j in range(len(chmm_labelFraction)):
            expectedFraction_mat[i][j] = segway_labelFraction[i] * chmm_labelFraction[j]

    observedFraction_mat = np.zeros(overlap_mat.shape)
    for i in range(len(segway_labelFraction)):
        for j in range(len(chmm_labelFraction)):
            observedFraction_mat[i][j] = overlap_mat[i][j] / total_bp

    obs_exp = np.divide(observedFraction_mat, expectedFraction_mat)
    
    # get the normalized value
    overlap_mat_colNorm = overlap_mat / np.sum(overlap_mat, axis=0)[np.newaxis, :] # chmm
    overlap_mat_rowNorm = overlap_mat / np.sum(overlap_mat, axis=1)[:, np.newaxis] # segway

    segway_cluster_list_old = list(overlap.index.values)
    
    # load the updated mnemonics
    label_term_mapping = {}
    mnemonics_file = sampleFolderAdd + 'mnemonics_v04.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term

    # change segway_cluster_list_updated based on the mnemonics
    segway_cluster_list = []
    for cluster in segway_cluster_list_old:
        label = cluster.split('_')[0]
        term = label_term_mapping[label]
        segway_cluster_list.append(label + '_' + term)

    # reorder the cluster list
    segway_cluster_list_reordered = []
    for label in segwayLabels:
        for cluster in segway_cluster_list:
            if cluster.split('_')[1] == label:
                segway_cluster_list_reordered.append(cluster)

    # add the ratio to the axis labels
    segwayAxis_list = []
    segEnh_inds = []
    for i, cluster in enumerate(segway_cluster_list):
        segwayAxis_list.append(cluster + '_' + str(round(segway_labelFraction[i], 4)))
        if cluster.split('_')[1] == 'Enhancer':
            segEnh_inds.append(i)

    segwayAxis_list_reordered = []
    for label in segwayLabels:
        for axis in segwayAxis_list:
            if axis.split('_')[1] == label:
                segwayAxis_list_reordered.append(axis)

    chmm_class_list = list(overlap.columns.values)

    chmmAxis_list = []
    for i, chmmclass in enumerate(chmm_class_list):
        chmmAxis_list.append(chmmclass + '_' + str(round(chmm_labelFraction[i], 4)))

    obs_exp_log = np.log10(obs_exp, out=np.zeros_like(obs_exp), where=(obs_exp!=0))
    obs_exp_log = np.where(obs_exp_log < 0, 0, obs_exp_log)

    sumCovVal = np.zeros(len(segwayAxis_list))
    meanLog = np.zeros(len(segwayAxis_list))
    for i in range(len(segwayAxis_list)):
        sumCovVal[i] = np.sum(overlap_mat_colNorm[i,0:6]* obs_exp_log[i, 0:6])
        meanLog[i] = np.mean(obs_exp_log[i, 0:6])

    totEnhancer = np.sum(overlap_mat[:,0:6], 1)
    chromEnhFracsSelf = totEnhancer /(sum(totEnhancer))

    selfMatch[annAccession] = (segwayAxis_list, sumCovVal, chromEnhFracsSelf, meanLog)

    #book = np.argsort(sumCovVal)
    #for i in book:
        #print('%s  %.4f' %(segwayAxis_list[i], sumCovVal[i]))

    #selfMatch[annAccession] = segwayAxis_list[book[-1]]
    #selfMatchValue[annAccession] = sumCovVal[book[-1]]

    #bestMatch = {}
    #bestMatchVal = {}
    matchList = {}
    matchVal = {}
    chromEnhFracs = {}
    meanLogs = {}
    for nonSelfAccession in annAccessionList:

        if nonSelfAccession == annAccession:
            continue

        if nonSelfAccession == 'ENCSR313QGL' or nonSelfAccession == 'ENCSR592IOP' or nonSelfAccession == 'ENCSR721USS' or nonSelfAccession == 'ENCSR273LUT' or nonSelfAccession == 'ENCSR699DMW':
            continue
        
        chmmFile = chmmFile_dict[nonSelfAccession]

        if chmmFile == 'none':
            continue

        print('here')

        chmmAccession = re.search('ENCFF.*\.', chmmFile)[0][:-1]
    
        inputFileName = 'overlap_segway_%s_chmm_%s_nonself.pkl' %(bedAccession, chmmAccession)
        inputFile = dataFolder + dataSubFolder + annAccession + '/' + inputFileName
        with open(inputFile, 'rb') as f:
            overlap = pickle.load(f)

        overlap_mat = overlap.to_numpy()

        print(inputFileName)

        # total bp count that is covered by both annotations
        total_bp = np.sum(overlap_mat)

        # fraction of base pairs in each label to the total bp count
        chmm_labelFraction = np.sum(overlap_mat, axis = 0) / total_bp
        segway_labelFraction = np.sum(overlap_mat, axis = 1) / total_bp

        # get the observed versus the expected fraction
        expectedFraction_mat = np.zeros(overlap_mat.shape)
        for i in range(len(segway_labelFraction)):
            for j in range(len(chmm_labelFraction)):
                expectedFraction_mat[i][j] = segway_labelFraction[i] * chmm_labelFraction[j]

        observedFraction_mat = np.zeros(overlap_mat.shape)
        for i in range(len(segway_labelFraction)):
            for j in range(len(chmm_labelFraction)):
                observedFraction_mat[i][j] = overlap_mat[i][j] / total_bp

        obs_exp = np.divide(observedFraction_mat, expectedFraction_mat)

        # get the normalized value
        overlap_mat_colNorm = overlap_mat / np.sum(overlap_mat, axis=0)[np.newaxis, :] # chmm
        overlap_mat_rowNorm = overlap_mat / np.sum(overlap_mat, axis=1)[:, np.newaxis] # segway

        segway_cluster_list_old = list(overlap.index.values)
    
        # load the updated mnemonics
        label_term_mapping = {}
        mnemonics_file = sampleFolderAdd + 'mnemonics_v04.txt'
        with open(mnemonics_file, 'r') as mnemonics:
            for line in mnemonics:
                #print(line)
                label = line.strip().split()[0]
                term = line.strip().split()[1]
                label_term_mapping[label] = term

        # change segway_cluster_list_updated based on the mnemonics
        segway_cluster_list = []
        for cluster in segway_cluster_list_old:
            label = cluster.split('_')[0]
            term = label_term_mapping[label]
            segway_cluster_list.append(label + '_' + term)

        # reorder the cluster list
        segway_cluster_list_reordered = []
        for label in segwayLabels:
            for cluster in segway_cluster_list:
                if cluster.split('_')[1] == label:
                    segway_cluster_list_reordered.append(cluster)

        # add the ratio to the axis labels
        segwayAxis_list = []
        segEnh_inds = []
        for i, cluster in enumerate(segway_cluster_list):
            segwayAxis_list.append(cluster + '_' + str(round(segway_labelFraction[i], 4)))
            if cluster.split('_')[1] == 'Enhancer':
                segEnh_inds.append(i)

        segwayAxis_list_reordered = []
        for label in segwayLabels:
            for axis in segwayAxis_list:
                if axis.split('_')[1] == label:
                    segwayAxis_list_reordered.append(axis)

        chmm_class_list = list(overlap.columns.values)

        chmmAxis_list = []
        for i, chmmclass in enumerate(chmm_class_list):
            chmmAxis_list.append(chmmclass + '_' + str(round(chmm_labelFraction[i], 4)))

        obs_exp_log = np.log10(obs_exp, out=np.zeros_like(obs_exp), where=(obs_exp!=0))
        obs_exp_log = np.where(obs_exp_log < 0, 0, obs_exp_log)

        sumCovVal = np.zeros(len(segwayAxis_list))
        meanLog = np.zeros(len(segwayAxis_list))
        for i in range(len(segwayAxis_list)):
            sumCovVal[i] = np.sum(overlap_mat_colNorm[i,0:6]* obs_exp_log[i, 0:6])
            meanLog[i] = np.mean(obs_exp_log[i, 0:6])

        totEnhancer = np.sum(overlap_mat[:,0:6], 1)
        chromEnhFrac = totEnhancer /(sum(totEnhancer))

        #book = np.argsort(sumCovVal)
        #for i in book:
            #print('%s  %.4f' %(segwayAxis_list[i], sumCovVal[i]))


        #bestMatch[nonSelfAccession] = segwayAxis_list[book[-1]]
        #bestMatchVal[nonSelfAccession] = sumCovVal[book[-1]]

        #matchList[nonSelfAccession] = segwayAxis_list
        matchVal[nonSelfAccession] = sumCovVal
        chromEnhFracs[nonSelfAccession] = chromEnhFrac
        meanLogs[nonSelfAccession] = meanLog

    #allMatch[annAccession] = (bestMatch, bestMatchVal)
    allMatch[annAccession] = (segwayAxis_list, matchVal, chromEnhFracs, meanLogs)

        # find the rows with label Enhancer (just that)
    # get the coverage ratio from the cells with log > .3

##### now parsing it all
########################################

allMatch 
selfMatch

selfCov = {}
allCov = {}
for annAccession in annAccessionList[0:30]:

    if not(annAccession in list(selfMatch.keys())):
        continue

    selfMatchThis = selfMatch[annAccession]

    sumCovVal = selfMatchThis[1]
    book = np.argsort(sumCovVal)[::-1]

    topIndsSelf = book[0:3]

    fracs = selfMatchThis[2]
    print(sum(fracs[topIndsSelf]))

    print(selfMatchThis[3][topIndsSelf])

    selfCov[annAccession] = sum(fracs[topIndsSelf])

    labels = selfMatchThis[0]
    for ind in topIndsSelf:
        print(labels[ind])

    allMatchThis = allMatch[annAccession]

    nonKeys = list(allMatchThis[1].keys())
 
    segLabels = allMatchThis[0]

    otherCov = {}
    for nonSelfAccession in nonKeys:

        sumCovVal = allMatchThis[1][nonSelfAccession]
        book = np.argsort(sumCovVal)[::-1]

        topInds = book[0:3]
        fracs = allMatchThis[2][nonSelfAccession]
        print(sum(fracs[topIndsSelf]))
        
        print(sum(fracs[topInds]))

        minMeanLogNonSelf = (allMatchThis[3][nonSelfAccession][topInds])
        
        print(minMeanLogNonSelf)
        print(selfMatchThis[3][topIndsSelf])


        for ind in topInds:
            print(labels[ind])

        print('----------------')
        otherCov[nonSelfAccession] = sum(fracs[topInds])

    allCov[annAccession] = otherCov


annAccession = annAccessionList[6]    
otherCov = list(allCov[annAccession].values())
mainCov = selfCov[annAccession]

sib = np.sort(otherCov/mainCov)

plt.hist(sib)
plt.show()


# TODO: what coverage do I get at that same level of significance (log ratio overlap)
# regarding the coverage, level of significance doe snot matter, as a Quies might have high significance level
# for the other samples, since, it is the other sample!!! BUT, where it is an active lable, we definitley have lower
# significance level for the other thing. That is why it doesn't matter 