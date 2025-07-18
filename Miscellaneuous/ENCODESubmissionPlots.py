# This code is to get the plots for ENCODE submission. There are 5, plots for each sample, the GMTK, the average, the probabilities, the average per term and the gene/transcriptomic plots. The code is similar to those in GMTK_parameter_plot_01.py and annotation_comparison_scaled.py. It also needs QC_transcriptionComparison_util.py
#
# ITEMS IN THE CODE:
# ########################################
# 0. Initials
# 1. The 4 panel plot
# 2. The transcription plot (that should be histograme like)

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

import glob

# submission plot folder:
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/submissionPlots/'

# General data folder
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'

# GTF data structure
fileName = dataFolder + '/geneLists.pkl'
with open(fileName, 'rb') as pickleFile:
    geneListsAndIDs = pickle.load(pickleFile)
    
geneList = geneListsAndIDs[0]
geneIDList = geneListsAndIDs[1]
del geneListsAndIDs

# Segway states:
#segwayStates = ['Enhancer', 'EnhancerLow', 'Promoter', 'PromoterFlanking', 'Transcribed', 'CTCF', 'K9K36', 'Bivalent', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent']

segwayStates = ['Promoter', 'PromoterFlanking', 'Enhancer', 'EnhancerLow', 'Bivalent', 'CTCF', 'Transcribed', 'K9K36', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent']
print(len(segwayStates))
segwayStateCount = len(segwayStates)

trackList = [ 'H3K4me3', 'H3K27ac', 'H3K4me1','H3K36me3', 'H3K27me3', 'H3K9me3', 'CTCF','DNase-seq','ATAC-seq','POLR2A','EP300']

# load the sample list
inputFileName = 'all235Annot_meta_corrected.pkl'
#inputFileName = 'all235Annot_meta.pkl'
inputFile = dataFolder +  inputFileName
with open(inputFile, 'rb') as f:
    allMeta = pickle.load(f)

accessionList = list(allMeta.keys())

########################################
# 1. The 4 panel plot
########################################

# p1: the average values
# p2: the probs
# p3: the GMTK
# p4: the average for the terms

# get the address of the GMTK files for each of the samples
runID_file = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/testBatch105/' + 'runID_accession_map_105run.pkl'
with open(runID_file, 'rb') as pickledFile:
    runID_accession105 = pickle.load(pickledFile)

accession_runID = {}
for runID in list(runID_accession105.keys()):
    ac = runID_accession105[runID]
    accession_runID[ac] = runID

GMTK_files = {}
for accession in accessionList:
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']
    if accession == 'ENCSR121RCG':
        continue

    if '38batch' in annotationFolder:
         GMTK_files[accession] = annotationFolder + 'segOutput/gmtk_parameters/gmtk_parameters.stats.csv'


    if 'May11' in annotationFolder:
        GMTK_files[accession] = annotationFolder + 'call-segtools/gmtk_parameters.stats.csv'

    if 'Batch105' in annotationFolder:
        runID = accession_runID[accession]
        GMTK_files[accession] = dataFolder + 'testBatch105/all_segtools/' + runID +  '/glob-03b7332b8fdb9a1ca33a23093d5878d5' + '/gmtk_parameters.stats.csv'

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


signal_files = {}
for accession in accessionList:
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']
    if accession == 'ENCSR121RCG':
        continue

    if '38batch' in annotationFolder:
        signal_files[accession] = annotationFolder + 'segOutput/call_signal_distribution/signal_distribution.tab'

    if 'May11' in annotationFolder:
        signal_files[accession] = annotationFolder + 'call-segtools/signal_distribution.tab'

    if 'Batch105' in annotationFolder:
        runID = accession_runID[accession]
        signal_files[accession] = dataFolder + 'testBatch105/all_segtools/' + runID +  '/signal_distribution/signal_distribution.tab'

# problem with 205 (no signal file)
plotFileName = {}
for accession in accessionList:

    annotation = allMeta[accession]
        
    annotationFolder = annotation['folder']

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

    # get the Segway cluster names and sorted Segway cluster list (for plot)
    # >>>>>>>>>>>>>>>>>>>>>>>>>>
    # segway cluster names
    segwayClusterNames_list_ordered = [] # name of segway cluster, that is label_state
    segwayClusterNames = {}
    segwayClusterIndex = {} # index of the sorted segway labels
    index = 0
    for state in segwayStates:
        for label in range(labelCount):
            this_state = label_term_mapping[str(label)]
            if this_state == state:
                clusterName = '%d_%s' %(label, state)
                segwayClusterNames_list_ordered.append(clusterName)
                segwayClusterIndex[str(label)] = index
                segwayClusterNames[str(label)] = clusterName
                index+=1

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
            print(line)
            track = track_assay_map[fields[1]]
            track_ind = sampleTrack_list_sorted.index(track) # track list order
            segway_cluster = segwayClusterNames[fields[0]]
            cluster_ind = segwayClusterIndex[fields[0]] # cluster list order
            signal_dist_mat[cluster_ind][track_ind] = round(float(fields[2]), 4)
            
    z_signal_dist_mat = stats.zscore(signal_dist_mat, axis = 0)

    # plot
     # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    from matplotlib import colors
    fig, axs = plt.subplots(2, 2, figsize=(12,8))
    h1 = pd.DataFrame(z_signal_dist_mat, index = segwayClusterNames_list_ordered, columns = sampleTrack_list_sorted)
    
    cmap = sns.diverging_palette(245, 15, s=80, l=40, sep=2, as_cmap=True)
    g1 = sns.heatmap(h1, center =0, cmap=cmap, ax=axs[0,0], linewidths=.1, linecolor='white')
    #plt.show()

    ## For discrete color map
    #cmap = plt.cm.coolwarm
    #norm = colors.BoundaryNorm(np.arange(-2.5, 3, 1), cmap.N)
    #g1 = sns.heatmap(h1, center = 0,cmap=cmap, norm=norm, ax=axs[0,0], linewidths=.1, linecolor='white')
    g1.set_title('mean signal value - zscore')


    # get the probs
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    probs_file = annotationFolder + 'probs_v04.csv'
    h2 = pd.read_csv(probs_file)
    h2.drop(columns=h2.columns[0], inplace=True, axis=1)
    for label in range(labelCount):
        h2.rename(index={int(label) : segwayClusterNames[str(label)]}, inplace=True)

    h2 = h2.reindex(segwayClusterNames_list_ordered)
    h2 = h2[segwayStates]
    cmap = sns.diverging_palette(240, 10, s=80, l=30, as_cmap=True)
    g2 = sns.heatmap(h2, center = 0,cmap=cmap, ax=axs[0,1], yticklabels=True)
    #g1.set_ylabel('this')
    #g1.set_xlabel('that')
    g2.set_title('Interpretation - probabilities')
    #plt.ylabel('segway')
    #plt.ylabel('chmm')
    #plt.subplots_adjust(left=.075, right=.95, top=.9, bottom=.25)
    plt.tight_layout()
    
    # get the GMTK
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    df = pd.read_csv(GMTK_files[accession])

    df = df[df.columns[1:]]

    track_accession = list(df)
    labels = list(df.index)

    # get column list
    track_names = []
    for track in track_accession:
        track_names.append(track_assay_map[track])

    # rename the columns of the dataframe
    df.set_axis(track_names, axis='columns', inplace=True)

    # reorder the columns of the dataframe
    df = df[sampleTrack_list_sorted]

    # rename the rows of dataframe
    labelList_initial = []
    for i in range(labelCount):
        labelList_initial.append(segwayClusterNames[str(i)])
    df.set_axis(labelList_initial, axis='rows', inplace=True)

    # reorder the rows of the dataframe
    df= df.reindex(segwayClusterNames_list_ordered)

    g3 = sns.heatmap(df, center = 0,cmap=cmap, ax=axs[1,0], yticklabels=True)
    #g1.set_ylabel('this')
    #g1.set_xlabel('that')
    g3.set_title('Segway emission probabilities (GMTK)')
    #plt.ylabel('segway')
    #plt.ylabel('chmm')
    #plt.subplots_adjust(left=.075, right=.95, top=.9, bottom=.25)
    plt.tight_layout()
    

    # genome coverage 
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    lfile = length_files[accession]
    barValues = np.zeros(labelCount)
    with open(lfile, 'r') as inputFile:
        header = inputFile.readline()
        header = inputFile.readline()
        for line in inputFile:
            fields = line.strip().split()
            label = fields[0]
            coverage = fields[-1]
            print(label, coverage)
            barValues[segwayClusterIndex[str(label)]] = coverage


    axs[1,1] = plt.bar(range(labelCount), barValues, .9, tick_label=segwayClusterNames_list_ordered, color=(.1,.1,.1,1))
    plt.xticks(rotation=90)
    plt.title('genome coverage')
    plt.tight_layout()
    
    plt.subplots_adjust( top=.9)
    fig.suptitle('Annotation accession: %s' %(accession))


    figFileName = '%s_samplePlot.pdf' %(accession)
    plotFileName[accession] = figFileName
    figFile = '%ssubmissionPlots/%s' %(plotFolder, figFileName)
    print(figFile)
    plt.savefig(figFile)
    plt.close('all')


# get the address of the GMTK files for each of the samples
accession_samplePlot = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/testBatch105/samplePlotMap.pkl'
with open(accession_samplePlot, 'wb') as f:
    pickle.dump(plotFileName, f)

########################################
# 2. The transcriptomic plots
########################################

# 2hrs ish - the plot is running, I am trying one.

# test file for new sample:
annotation = allMeta['ENCSR860AAJ']
geneFile = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/the38batch/ENCSR860AAJ/default_5kg_geneBodyCoverage_newSort.pkl'

countAcc = 0

count = 0
genePlotFileNames = {}
for accession in accessionList[206:]:

    annotation = allMeta[accession]
    #if ((('38batch' in annotation['folder']) or ('May11' in annotation['folder'])) and not(annotation['RNAseqFile'] == 'none')):
    if((annotation['RNAseqFile'] == 'none')):
        count +=1
        print(countAcc)
        annotationFolder = annotation['folder']
        geneFile = annotationFolder + 'default_5kg_geneBodyCoverage_newSort.pkl'

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
        
        # get the sorting index for the labels
        sortIndex = []
        for state in segwayStates:
            for i in range(labelCount):
                thisTerm = label_term_mapping[str(i)]
                if thisTerm == state:
                    sortIndex.append(i)

        with open(geneFile, 'rb') as file:
            geneMat = pickle.load(file)

        thisMat = geneMat

        # get the fraction of coverage for each basepair
        lfile = length_files[accession]
        fractions = np.zeros(labelCount)
        with open(lfile, 'r') as inputFile:
            header = inputFile.readline()
            header = inputFile.readline()
            for line in inputFile:
                fields = line.strip().split()
                label = fields[0]
                coverage = fields[-1]
                print(label, coverage)
                fractions[int(label)] = coverage

        # make it the ratio
        thisMat = thisMat / np.sum(thisMat, axis=0)[np.newaxis,:]
        # versus expected
        thisMat = thisMat / fractions[:, None]
        logMat = np.log10(thisMat, out=np.zeros_like(thisMat), where=(thisMat!=0))

        indexList = np.array(list(range(160)))

        fig, axs = plt.subplots(labelCount, 1, figsize=(5,8))
        xticks = [30, 130]
        xticksLabels = ['TSS', 'TTS']
        for i, label in enumerate(sortIndex):
            print(label)
            print(label_term_mapping[str(label)])
            positiveInds = indexList[logMat[label,:] >= 0]
            negativeInds = indexList[logMat[label,:] < 0]
            posBar = np.copy(logMat[label, :])
            posBar[negativeInds]=0
            axs[i].bar(range(160), posBar, color=[249/256, 80/256, 97/256,1], width=1)
            negBar = np.copy(logMat[label, :])
            negBar[positiveInds]=0
            axs[i].bar(range(160), negBar, color=[50/256,164/256,249/256,1], width=1)
            axs[i].set_ylim((-1.5,1.5))
            ylabel = str(label) + '_' + label_term_mapping[str(label)]
            axs[i].text(60, .65, ylabel, fontsize=8)
            axs[i].set_xticks(xticks)
            axs[i].set_yticks([-1, 1])
            axs[i].set_yticklabels([-1, 1], fontsize=8)
            #axs[i].set_ylabel(ylabel, rotation=0, fontsize=10, labelpad=3, ha='right', va='center')
            
        axs[labelCount-1].set_xticklabels(xticksLabels)
        axs[labelCount-1].set_xlabel('Position relative to gene')
        axs[0].set_title('Enrichment of labels at genes \n (log10(observed/expected))', fontsize = 8)

        #plt.tight_layout()
        plt.subplots_adjust(top=.9)
        fig.suptitle('Annotation accession: %s' %(accession))

        figFileName = '%s_geneBody.pdf' %(accession)
        genePlotFileNames[accession] = figFileName
        figFile = '%ssubmissionPlots/%s' %(plotFolder, figFileName)
        print(figFile)
        plt.savefig(figFile)
        plt.close('all')
        

# get the gene mapping for all the samples, for those without transcriptomic and those with transcriptomic

# Make the matrix several bar plots (negative positive stuff)

# plot plot

count = 0
for accession in accessionList:

    annotation = allMeta[accession]
    #if ((('38batch' in annotation['folder']) or ('May11' in annotation['folder'])) and not(annotation['RNAseqFile'] == 'none')):
    if((annotation['RNAseqFile'] != 'none')):
        count +=1
        annotationFolder = annotation['folder']
        expFile = annotationFolder + 'defaultExp_5kg_expSummary_newSort_Q30.pkl'
        print(expFile)
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

        
        # get the sorting index for the labels
        sortIndex = []
        for state in segwayStates:
            for i in range(labelCount):
                thisTerm = label_term_mapping[str(i)]
                if thisTerm == state:
                    sortIndex.append(i)

        # get the fraction of coverage for each basepair
        lfile = length_files[accession]
        fractions = np.zeros(labelCount)
        with open(lfile, 'r') as inputFile:
            header = inputFile.readline()
            header = inputFile.readline()
            for line in inputFile:
                fields = line.strip().split()
                label = fields[0]
                coverage = fields[-1]
                print(label, coverage)
                fractions[int(label)] = coverage


        fig, axs = plt.subplots(labelCount, 3, figsize=(12,8))
                
        xticks = [30, 130]
        xticksLabels = ['TSS', 'TTS']

        indexList = np.array(list(range(160)))

        thisMat = expMats['clusterMats'][0]

        # make it the ratio
        thisMat = thisMat / np.sum(thisMat, axis=0)[np.newaxis,:]
        # versus expected
        thisMat = thisMat / fractions[:, None]
        logMat = np.log10(thisMat, out=np.zeros_like(thisMat), where=(thisMat!=0))
        for i, label in enumerate(sortIndex):
            print(label)
            print(label_term_mapping[str(label)])
            positiveInds = indexList[logMat[label,:] >= 0]
            negativeInds = indexList[logMat[label,:] < 0]
            posBar = np.copy(logMat[label, :])
            posBar[negativeInds]=0
            axs[i,0].bar(range(160), posBar, color=[249/256, 80/256, 97/256,1], width=1)
            negBar = np.copy(logMat[label, :])
            negBar[positiveInds]=0
            axs[i,0].bar(range(160), negBar, color=[50/256,164/256,249/256,1], width=1)
            axs[i,0].set_ylim((-1.5,1.5))
            ylabel = str(label) + '_' + label_term_mapping[str(label)]
            #axs[i,0].text(60, .55, ylabel, fontsize=8)
            axs[i,0].set_xticks(xticks)
            axs[i,0].set_yticks([-1, 1])
            axs[i,0].set_yticklabels([-1, 1], fontsize=8)
            axs[i,0].set_ylabel(ylabel, rotation=0, fontsize=10, labelpad=3, ha='right', va='center')

        axs[labelCount-1,0].set_xticklabels(xticksLabels)
        axs[labelCount-1,0].set_xlabel('Position relative to gene')
        axs[0,0].set_title('Enrichment of labels at genes with \nzero expression (log10(observed/expected))', fontsize = 9)
        axs[labelCount-1,0].set_xlabel('Position relative to gene')

        thisMat = expMats['clusterMats'][1]

        # make it the ratio
        thisMat = thisMat / np.sum(thisMat, axis=0)[np.newaxis,:]
        # versus expected
        thisMat = thisMat / fractions[:, None]
        logMat = np.log10(thisMat, out=np.zeros_like(thisMat), where=(thisMat!=0))
        for i, label in enumerate(sortIndex):
            print(label)
            print(label_term_mapping[str(label)])
            positiveInds = indexList[logMat[label,:] >= 0]
            negativeInds = indexList[logMat[label,:] < 0]
            posBar = np.copy(logMat[label, :])
            posBar[negativeInds]=0
            axs[i,1].bar(range(160), posBar, color=[249/256, 80/256, 97/256,1], width=1)
            negBar = np.copy(logMat[label, :])
            negBar[positiveInds]=0
            axs[i,1].bar(range(160), negBar, color=[50/256,164/256,249/256,1], width=1)
            axs[i,1].set_ylim((-1.5,1.5))
            ylabel = str(label) + '_' + label_term_mapping[str(label)]
            #axs[i,1].text(60, .55, ylabel, fontsize=8)
            axs[i,1].set_xticks(xticks)
            axs[i,1].set_yticks([-1, 1])
            axs[i,1].set_yticklabels([-1, 1], fontsize=8)

        axs[labelCount-1,1].set_xticklabels(xticksLabels)
        axs[labelCount-1,1].set_xlabel('Position relative to gene')
        axs[0,1].set_title('Enrichment of labels at genes with \nbottom 30% expression (log10(observed/expected))', fontsize = 9)
        axs[labelCount-1,1].set_xlabel('Position relative to gene')


        thisMat = expMats['clusterMats'][2]

        # make it the ratio
        thisMat = thisMat / np.sum(thisMat, axis=0)[np.newaxis,:]
        # versus expected
        thisMat = thisMat / fractions[:, None]
        logMat = np.log10(thisMat, out=np.zeros_like(thisMat), where=(thisMat!=0))
        for i, label in enumerate(sortIndex):
            print(label)
            print(label_term_mapping[str(label)])
            positiveInds = indexList[logMat[label,:] >= 0]
            negativeInds = indexList[logMat[label,:] < 0]
            posBar = np.copy(logMat[label, :])
            posBar[negativeInds]=0
            axs[i,2].bar(range(160), posBar, color=[249/256, 80/256, 97/256,1], width=1)
            negBar = np.copy(logMat[label, :])
            negBar[positiveInds]=0
            axs[i,2].bar(range(160), negBar, color=[50/256,164/256,249/256,1], width=1)
            axs[i,2].set_ylim((-1.5,1.5))
            ylabel = str(label) + '_' + label_term_mapping[str(label)]
            #axs[i,2].text(60, .55, ylabel, fontsize=8)
            axs[i,2].set_xticks(xticks)
            axs[i,2].set_yticks([-1, 1])
            axs[i,2].set_yticklabels([-1, 1], fontsize=8)

        axs[labelCount-1,2].set_xticklabels(xticksLabels)
        axs[labelCount-1,2].set_xlabel('Position relative to gene')
        axs[0,2].set_title('Enrichment of labels at genes with \ntop 70% expression (log10(observed/expected))', fontsize = 9)
        axs[labelCount-1,2].set_xlabel('Position relative to gene')


        plt.subplots_adjust( top=.9)
        fig.suptitle('Annotation accession: %s' %(accession))

        figFileName = '%s_geneBody.pdf' %(accession)
        genePlotFileNames[accession] = figFileName
        figFile = '%ssubmissionPlots/%s' %(plotFolder, figFileName)
        print(figFile)
        plt.savefig(figFile)
        plt.close('all')


        # get the address of the GMTK files for each of the samples
accession_genePlot = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/genePlotMap.pkl'
with open(accession_genePlot, 'wb') as f:
    pickle.dump(genePlotFileNames, f)


# TODO: document the transcriptomic plot (code and stuff)


