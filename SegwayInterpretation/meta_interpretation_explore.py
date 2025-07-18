# in this code I want to go over the classifier reports and improvements
# ########################################
# 0. Initials
# 1. Go through each sample, review classifier info
# 2. Columns: ID, tissue name, track count, ctcf and DNase
#
#
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

# General data folder
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'

plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/testBatch105/'

dataSubFolder = 'testBatch105/fromAPI/'

# GTF data structure
fileName = dataFolder + '/geneLists.pkl'
with open(fileName, 'rb') as pickleFile:
    geneListsAndIDs = pickle.load(pickleFile)
    
geneList = geneListsAndIDs[0]
geneIDList = geneListsAndIDs[1]
del geneListsAndIDs

segwayLabels = ['Quiescent', 'ConstitutiveHet', 'FacultativeHet', 'Transcribed', 'Promoter', 'Enhancer', 'RegPermissive', 'Bivalent', 'LowConfidence']

inputFile = dataFolder + dataSubFolder + 'metaInfo.pkl'
with open(inputFile, 'rb') as pickledFile:
    annMeta = pickle.load(pickledFile)

sample_folder_list = list(annMeta.keys())

#######################################
# 1. Go through each sample, review classifier info
#######################################

outputFile = dataFolder + 'interpretation_summary.txt'
with open(outputFile, 'w') as output:
    output.write('%s\t%s\t%s\t%s\t%s\t%s\n' %('Accession', 'tissue', 'track_count', 'ATAC', 'DNase', 'CTCF'))

    for sampleFolder in sample_folder_list:

        sampleFolderAdd = dataFolder + dataSubFolder + sampleFolder + '/'

        # get the list of input tracks
        mapping_file = sampleFolderAdd + 'trackname_assay.txt'
        signal_file = sampleFolderAdd + 'signal_distribution.tab.txt'

        # read the mapping_file
        track_assay_map = {}
        inputTrack_list = [] # list of the track type for the sample
        ATAC = 0
        DNase = 0
        CTCF = 0
        with open(mapping_file) as inputFile:
            for line in inputFile:
                fields = line.strip().split()
                track_assay_map[fields[0]] = fields[1]
                inputTrack_list.append(fields[1])

                if 'CTCF' in fields[1]:
                    CTCF = 1
                if 'ATAC' in fields[1]:
                    ATAC = 1
                if 'DNase' in fields[1]:
                    DNase = 1

                track_count = len(inputTrack_list)

        output.write('%s\t%s\t%d\t%d\t%d\t%d\n' %(sampleFolder, 'sth', track_count, ATAC, DNase, CTCF))



#######################################
# 2. Do the classifier input heatmap plots
#######################################

inputFile = dataFolder + 'runID_accession_map_105run.pkl'
with open(inputFile, 'rb') as f:
    IDmap = pickle.load(f)

fileList = os.listdir(dataFolder + 'classifier_data/')

orderedCols = [ '(01) H3K9me3', '(02) H3K27me3',  '(03) H3K36me3', '(04) H3K4me3', '(05) H3K27ac', '(06) H3K4me1', "(07) 5' flanking (1000-10000 bp)", "(08) 5' flanking (1-1000 bp)", '(09) initial exon', '(10) initial intron', '(11) internal exons', '(12) internal introns', '(13) terminal exon', '(14) terminal intron',  "(15) 3' flanking (1-1000 bp)", "(16) 3' flanking (1000-10000 bp)", ]

for fileName in fileList:
    id = fileName.split('_')[0]
    accession = IDmap[id]
 
    figureFile = plotFolder + accession + '/classifier_data_thesix.pdf'

    classifier_data_file = dataFolder + 'classifier_data/' + fileName

    with open(classifier_data_file, 'r') as f:
        df = pd.read_csv(f, sep='\t')

    df.drop('orig_label', axis=1, inplace=True)

    #cols = list(df.columns.values)

    df = df[orderedCols]

    book = df[orderedCols[0:6]]

    fig, axs = plt.subplots(1,1, figsize=(6,6))

    #cmap = sns.diverging_palette(240, 10, s=80, l=30, as_cmap=True)
    cmap = sns.diverging_palette(240, 10, s=80, l=30, as_cmap=True)
    #g = sns.heatmap(book, center = 0,cmap=cmap, vmin=-2.5, vmax=3)
    g = sns.heatmap(book, center = 0, cmap=cmap)

    plt.tight_layout()
    plt.savefig(figureFile)
    print(figureFile)
    plt.close('all')

    plt.show()



