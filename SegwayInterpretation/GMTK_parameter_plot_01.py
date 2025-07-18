# this is an update to the GMTK_parameter_plot.py
# 0. Initials
# 

########################################
# 0. Initials
########################################

import util # this is for features_from_segtools_dir
import gzip
import pickle
import pandas as pd

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'

# >>>>>>>>>>>>>>>>>>>> for the 105 batch
dataSubFolder = 'testBatch105/fromAPI/'

inputFileName = 'all_annInfo_list.pkl'
inputFile = dataFolder + dataSubFolder + outputFileName
with open(inputFile, 'rb') as f:
    ann_info_list = pickle.load(f)

runID_file = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/testBatch105/' + 'runID_accession_map_105run.pkl'

with open(runID_file, 'rb') as pickledFile:
    runID_accession = pickle.load(pickledFile)

accession_runID = {}
for runID in list(runID_accession.keys()):
    ac = runID_accession[runID]
    accession_runID[ac] = runID

# <<<<<<<<<<<<<<<<<<<<



# >>>>>>>>>>>>>>>>>>>> for the 112 home run batch
dataSubFolder = 'the112batch/'
inputFile = dataFolder + dataSubFolder + 'accessionList.pkl'
with open(inputFile, 'rb') as f:
    accessionList = pickle.load(f)

# <<<<<<<<<<<<<<<<<<<<

# Segway labels for version 4 classifier
segwayLabels = ['Enhancer', 'EnhancerLow', 'Promoter', 'PromoterFlanking', 'Transcribed', 'CTCF', 'K9K36', 'Bivalent', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent']

track_order = [ 'H3K4me3', 'H3K27ac', 'H3K4me1', 'H3K36me3', 'H3K27me3', 'H3K9me3', 'CTCF', 'DNase-seq', 'ATAC-seq', 'POLR2A', 'EP300']

########################################
# 0. For the 105 runs 
########################################

# NOTE: the Segway labels are for classifier version 4, but this code was for classifer 2, perhaps change the files. 

for index in range(105):
    ann = ann_info_list[index]
    annAccession = ann['accession']
    print(annAccession)

    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
    print(sampleFolderAdd)

    runID = accession_runID[annAccession]

    # get the label from label from mnemonics: a dictionary from labels to terms
    label_term_mapping = {}
    mnemonics_file = sampleFolderAdd + 'mnemonics_v02.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term


    segtools_add = dataFolder + 'testBatch105/all_segtools/' + runID
    print(segtools_add)

    gmtk_file = segtools_add + '/glob-03b7332b8fdb9a1ca33a23093d5878d5' + '/gmtk_parameters.stats.csv'
    df = pd.read_csv(gmtk_file)

    df = df[df.columns[1:]]

    track_accession = list(df)
    labels = list(df.index)

    # get column list
    track_names = []
    for track in track_accession:
        track_names.append(ann['track_assay_map'][track])

    # get row list - this section is obsolete since we got the new interpretation terms that are different
    #ordered_labels = list(df.index)
    #segway_labels = list(ann['segway_anns']['clusters'].keys())
    #for label in segway_labels:
    #    number = int(label.split('_')[0])
    #    ordered_labels[number] = label

    index_list = []
    for i in range(len(df)):
        label_name = str(i) + '_' + label_term_mapping[str(i)]
        index_list.append(label_name)

    df.columns = track_names
    df.index = index_list

    reordered_columns = []
    for track in track_order:
        if track in track_names:
            reordered_columns.append(track)

    reordered_index = []
    for label in segwayLabels:
        for index in index_list:
            if index.split('_')[1] == label:
                reordered_index.append(index)
    
    df = df[reordered_columns]
    df = df.reindex(reordered_index)
    
    kado = df.transpose()

    fig,axs = plt.subplots(1, figsize=(4,3))
    sns.heatmap(kado, cmap='vlag')
    sns.set(font_scale=.7)
    plt.tight_layout()
    plt.title(annAccession)

    figFile = plotFolder + annAccession + '/GMTK_emmission.pdf'
    print(figFile)
    plt.savefig(figFile)
    plt.close('all')


plt.show()

plots_add = plotFolder + annAccession
print(plots_add)


########################################
# 1. For the 112 home runs
########################################

plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/batch112/'

for annAccession in accessionList:
    samplePlotFolder = plotFolder + annAccession
    os.mkdir(samplePlotFolder)

for annAccession in accessionList[2:15]:

    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
    print(sampleFolderAdd)

    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
    print(sampleFolderAdd)

    # get the label from label from mnemonics: a dictionary from labels to terms
    label_term_mapping = {}
    mnemonics_file = sampleFolderAdd + 'mnemonics_v04.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term

    mapping_file = sampleFolderAdd + 'trackname_assay.txt'
    track_assay_map = {}
    with open(mapping_file) as inputFile:
        for line in inputFile:
            fields = line.strip().split()
            track_assay_map[fields[0]] = fields[1]

    gmtk_file = sampleFolderAdd + 'segOutput/gmtk_parameters/gmtk_parameters.stats.csv'
    df = pd.read_csv(gmtk_file)

    df = df[df.columns[1:]]

    track_accession = list(df)
    labels = list(df.index)

    # get column list
    track_names = []
    for track in track_accession:
        track_names.append(track_assay_map[track])

    # get row list - this section is obsolete since we got the new interpretation terms that are different
    #ordered_labels = list(df.index)
    #segway_labels = list(ann['segway_anns']['clusters'].keys())
    #for label in segway_labels:
    #    number = int(label.split('_')[0])
    #    ordered_labels[number] = label

    index_list = []
    for i in range(len(df)):
        label_name = str(i) + '_' + label_term_mapping[str(i)]
        index_list.append(label_name)

    df.columns = track_names
    df.index = index_list

    reordered_columns = []
    for track in track_order:
        if track in track_names:
            reordered_columns.append(track)

    reordered_index = []
    for label in segwayLabels:
        for index in index_list:
            if index.split('_')[1] == label:
                reordered_index.append(index)
    
    df = df[reordered_columns]
    df = df.reindex(reordered_index)
    
    kado = df.transpose()

    fig,axs = plt.subplots(1, figsize=(4,3))
    sns.heatmap(kado, cmap='vlag')
    sns.set(font_scale=.7)
    plt.tight_layout()
    plt.title(annAccession)

    figFile = plotFolder + annAccession + '/GMTK_emmission.pdf'
    print(figFile)
    plt.savefig(figFile)
    plt.close('all')

# The plots look off, labels look almost completely random. I need to investigate the feature Matrix and feature aggregation. It looks like as if the feature matrix is feeded wrongly to the model, but it was fine with the previous batches so I am not sure what could have happened.
# Let's compare a previous data feature matrix with a current data feature Matrix and the resulting labels. Check to make sure all inputsa re right.

# it might be the order to the track assay mapping. Or not, no idea. No, it is not related. 

plt.show()

plots_add = plotFolder + annAccession
print(plots_add)


