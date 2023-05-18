# The new apply function (instead of the one from the Max's code)
# 0. Initials
# 1. For the 105 run batch, with the updated classifier
# 2. for the home run the38batch
# 3. for the Mayrun the92batch

########################################
# 0. Initials
########################################

import util # this is for features_from_segtools_dir
import gzip
import pickle
import pandas as pd

# load the model
runFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/runs/'
model_file = runFolder + "run02/model_278exp_reg0.028_auc0.86_wholeData.pickle.gz"
model_file = runFolder + "run03/model_300_reg.020_auc0.89V04.pickle.gz"
with gzip.open(model_file, "r") as f:
    the_model = pickle.load(f)


# the original 16 features that go to the classifier - just keeping it here
feature_names = ['(09) initial exon', '(01) H3K9me3', '(10) initial intron', '(02) H3K27me3', '(11) internal exons', '(04) H3K4me3', "(16) 3' flanking (1000-10000 bp)", '(12) internal introns', '(03) H3K36me3', '(13) terminal exon', '(06) H3K4me1', '(14) terminal intron', "(07) 5' flanking (1000-10000 bp)", '(05) H3K27ac', "(15) 3' flanking (1-1000 bp)", "(08) 5' flanking (1-1000 bp)"]

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
dataSubFolder = 'testBatch105/fromAPI/'

# for each sample we want the prob and the labels be generated. See how the probs file was made.

inputFileName = 'all_annInfo_list.pkl'
inputFile = dataFolder + dataSubFolder + inputFileName
with open(inputFile, 'rb') as f:
    annInfo_list = pickle.load(f)

########################################
# 1. For the 105 run batch, with the updated classifier
########################################

for ann in annInfo_list:

    index = ann['index']
    print(index)
    sampleFolder =  ann['accession']

    sampleFolderAdd = dataFolder + dataSubFolder + sampleFolder + '/'

    signal_file = sampleFolderAdd + 'signal_distribution.tab.txt'
    feature_file = sampleFolderAdd + 'feature_aggregation.tab.txt'
    mapping_file = sampleFolderAdd + 'trackname_assay.txt'
    track_assay_map = {}
    inputTrack_list = []
    with open(mapping_file) as inputFile:
        for line in inputFile:
            fields = line.strip().split()
            track_assay_map[fields[0]] = fields[1]
            inputTrack_list.append(fields[1])

    ann_features, ann_label_bases, ann_feature_names = features_from_segtools_dir(feature_file, signal_file, track_assay_map)

    df = pd.DataFrame(ann_features)
    dft = df.T
    dftr = dft[feature_names]

    labels = the_model.predict(dftr)
    mnemonics_file = sampleFolderAdd + 'mnemonics_v04.txt'
    data = np.column_stack([dftr.index, labels])
    np.savetxt(mnemonics_file, data, fmt=['%s\t', '%s'])
    # write the mnemonics file
    
    probs = np.round(the_model.predict_proba(dftr), 5)
    probs_file = sampleFolderAdd + 'probs_v04.csv'
    probsdf = pd.DataFrame(probs, index=range(probs.shape[0]), columns=the_model.classes_)
    probsdf.to_csv(probs_file)
    # write the mnemonics file

##### this is just to get the data matrix
from scipy.stats import zscore
term_list = ['Enhancer', 'Enhancer_low', 'Promoter','Promoter_flanking', 'CTCF', 'K9K36', 'Bivalent', 'Transcribed', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent']

term_sum = np.zeros((12, 16))
term_count = np.zeros(12,)

for ann in annInfo_list:

    index = ann['index']
    print(index)
    sampleFolder =  ann['accession']

    sampleFolderAdd = dataFolder + dataSubFolder + sampleFolder + '/'

    signal_file = sampleFolderAdd + 'signal_distribution.tab.txt'
    feature_file = sampleFolderAdd + 'feature_aggregation.tab.txt'
    mapping_file = sampleFolderAdd + 'trackname_assay.txt'
    track_assay_map = {}
    inputTrack_list = []
    with open(mapping_file) as inputFile:
        for line in inputFile:
            fields = line.strip().split()
            track_assay_map[fields[0]] = fields[1]
            inputTrack_list.append(fields[1])
    ann_features, ann_label_bases, ann_feature_names = features_from_segtools_dir(feature_file, signal_file, track_assay_map)

    df = pd.DataFrame(ann_features)
    dft = df.T
    dftr = dft[feature_names_reorder]

    plot_data_z = zscore(dftr, axis = 0)
    plot_data_z_thr = np.where(plot_data_z > 2, 2, plot_data_z)
    plot_data_z_thr = np.where(plot_data_z_thr < -2, -2, plot_data_z_thr)
    label_term_mapping = {}
    mnemonics_file = sampleFolderAdd + 'mnemonics_v02.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term

    for i_label in range(plot_data_z_thr.shape[0]):
        term = label_term_mapping[str(i_label)]
        term_ind = term_list.index(term)
        label_feature = plot_data_z_thr[i_label,]
        term_sum[term_ind,] += label_feature
        term_count[term_ind] += 1
        

summary_data = [term_list, term_sum, term_count]
outputFileName = 'summary_features.pkl'
outputFile = dataFolder + dataSubFolder + outputFileName
with open(outputFile, 'wb') as f:
    pickle.dump(summary_data, f)

with open(outputFile, 'rb') as f:
    summary_data = pickle.load(f)

term_sum = summary_data[1]
term_count = summary_data[2]
term_mean = np.zeros(term_sum.shape)
for i in range(len(term_list)):
    print(i)
    term_mean[i,] = term_sum[i,] / term_count[i]


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.color_palette("vlag", as_cmap=True)
sns.heatmap(term_mean, cmap='coolwarm')
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/testBatch105/'
figFile = plotFolder + 'term_feature_heatmap.pdf'
print(figFile)
plt.savefig(figFile)
plt.close('all')

sns.heatmap(plot_data_z_thr)
plt.show()
plt.close('all')


########################################
# 2. for the home run (the last 38)
########################################

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
dataSubFolder = 'the38batch/'

##### IMPORTANT: classifier training data
feature_names = ['(09) initial exon', '(01) H3K9me3', '(10) initial intron', '(02) H3K27me3', '(11) internal exons', '(04) H3K4me3', "(16) 3' flanking (1000-10000 bp)", '(12) internal introns', '(03) H3K36me3', '(13) terminal exon', '(06) H3K4me1', '(14) terminal intron', "(07) 5' flanking (1000-10000 bp)", '(05) H3K27ac', "(15) 3' flanking (1-1000 bp)", "(08) 5' flanking (1-1000 bp)"]

# load the address file
inputFile = dataFolder + dataSubFolder + 'hg_accessionList.pkl'
with open(inputFile, 'rb') as f:
    accessionList = pickle.load(f)

for i,accession in enumerate(accessionList[15:]):

    index = i
    print(index)
    print(accession)

    sampleFolderAdd = dataFolder + dataSubFolder + accession + '/'

    signal_file = sampleFolderAdd + 'segOutput/call_signal_distribution/signal_distribution.tab'
    feature_file = sampleFolderAdd + 'segOutput/call_feature_aggregation/feature_aggregation.tab'
    mapping_file = sampleFolderAdd + 'trackname_assay.txt'
    track_assay_map = {}
    inputTrack_list = []
    with open(mapping_file) as inputFile:
        for line in inputFile:
            fields = line.strip().split()
            track_assay_map[fields[0]] = fields[1]
            inputTrack_list.append(fields[1])

    ann_features, ann_label_bases, ann_feature_names = features_from_segtools_dir(feature_file, signal_file, track_assay_map)

    df = pd.DataFrame(ann_features)
    dft = df.T
    dftr = dft[feature_names]
    rowCount = dftr.shape[0]
    inds = np.linspace(0,rowCount-1, rowCount).astype(int)
    indList = [str(x) for x in inds]
    dftr = dftr.reindex(indList)

    labels = the_model.predict(dftr)
    mnemonics_file = sampleFolderAdd + 'mnemonics_v04.txt'
    data = np.column_stack([dftr.index, labels])
    np.savetxt(mnemonics_file, data, fmt=['%s\t', '%s'])
    # write the mnemonics file
    
    probs = np.round(the_model.predict_proba(dftr), 5)
    probs_file = sampleFolderAdd + 'probs_v04.csv'
    probsdf = pd.DataFrame(probs, index=range(probs.shape[0]), columns=the_model.classes_)
    probsdf.to_csv(probs_file)
    # write the mnemonics file

########################################
# 3. for the May run
########################################

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
dataSubFolder = 'testBatch_May112022/'

##### IMPORTANT: classifier training data
feature_names = ['(09) initial exon', '(01) H3K9me3', '(10) initial intron', '(02) H3K27me3', '(11) internal exons', '(04) H3K4me3', "(16) 3' flanking (1000-10000 bp)", '(12) internal introns', '(03) H3K36me3', '(13) terminal exon', '(06) H3K4me1', '(14) terminal intron', "(07) 5' flanking (1000-10000 bp)", '(05) H3K27ac', "(15) 3' flanking (1-1000 bp)", "(08) 5' flanking (1-1000 bp)"]

runIDFile = 'runID_list.pkl'

inputFile = dataFolder + dataSubFolder + runIDFile
with open(inputFile, 'rb') as f:
    runIDs = pickle.load(f)

for i,runID in enumerate(runIDs[46:]):

    index = i
    print(index)

    #if runID == '5857d68b-e559-4776-9c12-a6e10aea7f76': # this runID doesn't have the signal_dist file
    #    continue

    sampleFolderAdd = dataFolder + dataSubFolder + runID + '/'

    signal_file = sampleFolderAdd + 'call-segtools/signal_distribution.tab'
    feature_file = sampleFolderAdd + 'call-segtools/feature_aggregation.tab'
    mapping_file = sampleFolderAdd + 'trackname_assay.txt'
    track_assay_map = {}
    inputTrack_list = []
    with open(mapping_file) as inputFile:
        for line in inputFile:
            fields = line.strip().split()
            track_assay_map[fields[0]] = fields[1]
            inputTrack_list.append(fields[1])

    ann_features, ann_label_bases, ann_feature_names = features_from_segtools_dir(feature_file, signal_file, track_assay_map)

    
    df = pd.DataFrame(ann_features)
    dft = df.T
    dftr = dft[feature_names]

    labels = the_model.predict(dftr)
    mnemonics_file = sampleFolderAdd + 'mnemonics_v04.txt'
    data = np.column_stack([dftr.index, labels])
    np.savetxt(mnemonics_file, data, fmt=['%s\t', '%s'])
    # write the mnemonics file
    
    probs = np.round(the_model.predict_proba(dftr), 5)
    probs_file = sampleFolderAdd + 'probs_v04.csv'
    probsdf = pd.DataFrame(probs, index=range(probs.shape[0]), columns=the_model.classes_)
    probsdf.to_csv(probs_file)
    # write the mnemonics file

# writing all the mnemonic files into one file
ID_mnemonics = {}
for runID in runIDs:

    print(runID)

    if runID == '5857d68b-e559-4776-9c12-a6e10aea7f76': # this runID doesn't have the signal_dist file
        continue

    sampleFolderAdd = dataFolder + dataSubFolder + runID + '/'
    mnemfile = sampleFolderAdd + 'mnemonics_v03.txt'
    data = pd.read_csv(mnemfile, '\t', header=None)
    book = data.to_dict()
    ID_mnemonics[runID] = book[1]


outputFile = dataFolder + dataSubFolder + 'classifier_output.pkl'
with open(outputFile, 'wb') as output:
    pickle.dump(ID_mnemonics, output)


 # the original 16 features that go to the classifier - just keeping it here
feature_names = ['(09) initial exon', '(01) H3K9me3', '(10) initial intron', '(02) H3K27me3', '(11) internal exons', '(04) H3K4me3', "(16) 3' flanking (1000-10000 bp)", '(12) internal introns', '(03) H3K36me3', '(13) terminal exon', '(06) H3K4me1', '(14) terminal intron', "(07) 5' flanking (1000-10000 bp)", '(05) H3K27ac', "(15) 3' flanking (1-1000 bp)", "(08) 5' flanking (1-1000 bp)"]

feature_names_reorder = ['(09) initial exon', "(08) 5' flanking (1-1000 bp)", '(10) initial intron','(11) internal exons', '(12) internal introns', '(13) terminal exon', '(14) terminal intron', "(15) 3' flanking (1-1000 bp)",  "(07) 5' flanking (1000-10000 bp)", "(16) 3' flanking (1000-10000 bp)", '(04) H3K4me3','(05) H3K27ac',  '(06) H3K4me1','(03) H3K36me3', '(02) H3K27me3', '(01) H3K9me3']


# draft    
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # get signal_dist file
    signal_file = '/Users/marjanfarahbod/Downloads/signal_distribution.tab'
    # get feature_agg file
    feature_file = '/Users/marjanfarahbod/Downloads/feature_aggregation.tab'
    # get trackname_assay file
    mapping_file = '/Users/marjanfarahbod/Downloads/trackname_assay_meh.txt'
    track_assay_map = {}
    inputTrack_list = []
    with open(mapping_file) as inputFile:
        for line in inputFile:
            fields = line.strip().split()
            track_assay_map[fields[0]] = fields[1]
            inputTrack_list.append(fields[1])

        ann_features, ann_label_bases, ann_feature_names = util.features_from_segtools_dir(feature_file, signal_file, track_assay_map)

        df = pd.DataFrame(ann_features)
        dft = df.T
        dftr = dft[feature_names]

        labels = the_model.predict(dftr)
        # save the label file

