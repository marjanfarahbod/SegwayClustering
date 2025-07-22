#
# 1. Get the old classifier training data 
## 1.0 Get the old labels
## 1.1 Get the old featrues and fetch their labels based on the label_mappings
## 1.2 Get the classifier old data (from saved updated file)
# 2. Get data from the 105 sample runs, getting the feature matrix.
# 2.0 sample check plot (old vs new classifier data check)
# 3. Training the classifier with the old data
# 4. Extend the training data with the new data
# 4.1 leave one out cross validation
# 4.2 train the model with the whole training data and save the model
# 5. Apply the model on the 105 samples
# 5.0 Get the new labels
# 5.1 The clustering heatmap plot
# 5.2 The probabilities heatmap plot
# 100. new classifier data for Abe
# OBSOLETE

########################################
# 0. Initiation and set up
########################################

import sys
import os
import argparse
import subprocess
import gzip
import math
import random
from path import Path
#from pathlib import Path
from collections import defaultdict
#import bedtools
import shutil
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from util import *

from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.tree import DecisionTreeClassifier, export_graphviz
from sklearn.ensemble import RandomForestClassifier
import numpy
import pandas as pd
import pickle
from copy import deepcopy
import pandas as pd

import psutil
import resource


dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/testBatch105/'
runFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/runs/'

dataSubFolder = 'testBatch105/fromAPI/'

# get the file: classifier_data.tab

classifier_tab_fname = dataFolder + 'classifier_data.tab' # this is the old classifier data

#########################################
# 1. Get the old classifier training data
#########################################

## 1.0 Get the old labels
#########################################
'''
In this part, only the mapping of the labels are made and the 
paring of the features and labels remains for the next part
'''

'''
I will have the new data in this one as well
'''
# these are the old labels (from 2019)
label_mapping = dataFolder + 'label_mappings.txt'
label_mappings = {}
labels_set = set()
orig_labels_set = set()
with open(label_mapping, "r") as f:
    for line in f:
        line = line.split()
        if line[0] == "concatenation_key":
            continue
        concatenation_key = line[0]
        orig_label = line[1]
        label = line[2]
        if not concatenation_key in label_mappings:
            label_mappings[concatenation_key] = {}
        label_mappings[concatenation_key][orig_label] = label
        labels_set.add(label)
        orig_labels_set.add(orig_label)

all_ref_bio_labels = set.union(*map(set, map(lambda x: x.values(), label_mappings.values())))

# to get the first N biolabels
# these are the new labels (from 2022)
label_mapping = runFolder + 'run02/label_mappings_trainSelect.csv'
N = 295 # the old sample biolables
label_mappings = {}
new_label_archive = {}
#orig_labels_set = set()
with open(label_mapping, "r") as f:
    for i in range(N):
        line = next(f).strip().split(',')
        print(line)
        #line = next(f).strip().split('\t')
        if line[0] == "concatenation_key":
            continue
        concatenation_key = line[0]
        orig_label = line[1]
        label = line[2]
        new_label = line[-1]
        if not concatenation_key in label_mappings:
            label_mappings[concatenation_key] = {}
            new_label_archive[concatenation_key] = {}
            
        label_mappings[concatenation_key][orig_label] = label
        new_label_archive[concatenation_key][orig_label] = new_label
        #orig_labels_set.add(orig_label)

all_ref_bio_labels = set.union(*map(set, map(lambda x: x.values(), label_mappings.values())))

## 1.1 Get the old featrues and fetch their labels based on the label_mappings
#######################################################

'''
We need to run this part to get the example_bio_labels
'''

classifier_data_frame = pd.read_csv(classifier_tab_fname, sep="\t")

# extracting the classifier label from the file
example_bio_labels = [None for i in range(classifier_data_frame.shape[0])]
for i in range(classifier_data_frame.shape[0]):
    concatenation_key = classifier_data_frame.concatenation_key[i]
    orig_label = classifier_data_frame.orig_label[i]
    if orig_label in label_mappings[concatenation_key]:
        example_bio_labels[i] = label_mappings[concatenation_key][orig_label]
    else:
        logger.warning("No label mapping for {concatenation_key} {orig_label}".format(**locals()))
        example_bio_labels[i] = "??"

labels_set = set(example_bio_labels)
example_bio_labels = numpy.array(example_bio_labels) # classifier label
feature_names = numpy.array(classifier_data_frame.drop("orig_label", 1).drop("concatenation_key",1).columns)
example_orig_labels = numpy.array(classifier_data_frame.orig_label) # we don't need to use this

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> the new count
print(labels_set)
segwayLabels = list(labels_set)

label_training_counts_223 = np.zeros([len(labels_set),1])
for label in model_labels_extended:#example_bio_labels:
    ind = segwayLabels.index(label)
    label_training_counts_223[ind] += 1
    
print(label_training_counts_223)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> the old count
# old labels segwayLabels = ['Quiescent', 'ConstitutiveHet', 'FacultativeHet', 'Transcribed', 'Promoter', 'Enhancer', 'RegPermissive', 'Bivalent', 'LowConfidence']
# note: the 'LowConfidence' up there, is called 'Unclassified' in the labels
label_training_counts_223 = np.zeros([9,1])
for label in example_bio_labels:
    if label == 'Unclassified':
        label_training_counts_223[8] +=1
    else: 
        ind = segwayLabels.index(label)
        label_training_counts_223[ind] += 1

print(classifier_data_frame.mean(axis = 0))

################# Hard Coding Feature Order ####################
feature_names = ['(09) initial exon', '(01) H3K9me3', '(10) initial intron', '(02) H3K27me3', '(11) internal exons', '(04) H3K4me3', "(16) 3' flanking (1000-10000 bp)", '(12) internal introns', '(03) H3K36me3', '(13) terminal exon', '(06) H3K4me1', '(14) terminal intron', "(07) 5' flanking (1000-10000 bp)", '(05) H3K27ac', "(15) 3' flanking (1-1000 bp)", "(08) 5' flanking (1-1000 bp)"]

# reordering the features 
example_features = features_frame_to_matrix(classifier_data_frame, feature_names)
example_bio_labels = numpy.array(example_bio_labels) # classifier label

# Saving the original training data
model_labels = example_bio_labels
model_features = example_features

sampleNames = classifier_data_frame['concatenation_key'] # this for keep track of train samples

## 1.2 Getting the classifier old data (from saved updated file)
##################################################

original_training_data_file = runFolder + 'run02/trainingData_223_2022labels_00.pkl'
original_training_data = {}
original_training_data['model_labels'] = model_labels
original_training_data['model_features'] = model_features
with open(original_training_data_file, 'wb') as f:
    pickle.dump(original_training_data, f)

with open(original_training_data_file, 'rb') as f:
    original_training_data = pickle.load(f)

model_labels = original_training_data['model_labels']
model_features = original_training_data['model_features']
 
for i,label in enumerate(model_labels):
    if label == 'Enhancer_low':
        model_labels[i] = 'EnhancerLow'
    if label == 'Promoter_flanking':
        model_labels[i] = 'PromoterFlanking'

for i, label in enumerate(model_labels):
    if label == 'Unclassified':
        print(label)
        print(i)

del sampleNames[1] # only unclassified was index 1

# deleting the label 'Unclassified'
model_labels_del = np.delete(model_labels, 1)
model_features_del = np.delete(model_features, 1, 0)

model_labels = model_labels_del
model_features = model_features_del

#########################################
# 2. Get data from the 105 sample runs, getting the feature matrix.
#########################################

inputFileName = 'all_annInfo_list.pkl'
inputFile = dataFolder + dataSubFolder + inputFileName
with open(inputFile, 'rb') as f:
    ann_info_list = pickle.load(f)

allFeatureAgg_mat = np.zeros((105*16, 16)) # this one fills with index
agg_cluster_list = []
afa_index = 0
for ann in ann_info_list:

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

    for cluster in ann_features:
        ann_feature_names = list(ann_features[cluster].keys())
        for ann_feature in ann_feature_names:
            # feature = ann_feature[5:] # stripping the first 4 characters
            value = ann_features[cluster][ann_feature]
            if ann_feature in feature_names:
                feature_index = feature_names.index(ann_feature)
                allFeatureAgg_mat[afa_index, feature_index] = value
        index_cluster = '%d___%s' %(index, cluster)
        agg_cluster_list.append(index_cluster)
        afa_index +=1

allFeatureAgg_mat = allFeatureAgg_mat[0:afa_index,] # this has all the features to the order of feature_names (like the classifier data)

# 2.0 sample check plot (old vs new classifier data check)
##################################################

# ** plotting a set of samples just to see how different it is from the old classifier data
df = pd.DataFrame(allFeatureAgg_mat[0:100,])
classifier_data_frame # the other one

fig, axs = plt.subplots(1, 2, figsize=[8, 10])

# plotting the transcript data
cmap = sns.diverging_palette(240, 10, s=80, l=30, as_cmap=True)
sns.heatmap(df,  ax=axs[0])
df2 = classifier_data_frame
df3 = df2.drop(['orig_label', 'concatenation_key'], axis = 1)
sns.heatmap(df3, ax=axs[1])

plt.show()
plt.close('all')

##################################################
# 3. Training the classifier with the old data
##################################################

def make_model(reg): 
    return RandomForestClassifier(min_samples_leaf=reg, criterion="entropy")

#reg = 1e-2
reg = .067 # it seems that reg 15 or 15/223 is the best one. I will balance the samples and train the model then.

model = make_model(reg)
model.fit(model_features, model_labels)

# save the model
'''
model_file = runFolder + "model_223exp_reg15_auc0.75.pickle.gz"
with gzip.open(model_file, "w") as f:
    pickle.dump(model, f)
'''

from sklearn.metrics import precision_recall_fscore_support as score
import numpy as np
precision, recall, fscore, support = score(model_labels, model.predict(model_features), labels=['Quiescent', 'ConstitutiveHet', 'FacultativeHet', 'Transcribed', 'Promoter', 'Enhancer', 'RegPermissive', 'Bivalent', 'Unclassified'])
print('precision: {}'.format(precision))
print('recall: {}'.format(recall))
print('fscore: {}'.format(fscore))
from sklearn.metrics import accuracy_score
Accuracy = accuracy_score(model_labels, model.predict(model_features))
print("Accuracy: {}".format(Accuracy))
print(list(model.predict(model_features)))
print(list(model_labels))

# save the model
model_file = runFolder + "run02/model_223exp_reg15_auc0.78.pickle.gz"
with gzip.open(model_file, "w") as f:
    pickle.dump(model, f)

# save the confusion matrix
from sklearn.metrics import confusion_matrix
confusion_mat = confusion_matrix(model_labels, model.predict(model_features))
print(confusion_mat)
import matplotlib.pyplot as plt
from sklearn.metrics import plot_confusion_matrix
plot_confusion_matrix(model, model_features, model_labels, xticks_rotation='vertical')
plt.grid(False)
plt.tight_layout()
#plt.show()
pltFile = runFolder + 'run01/model_223exp_reg15_auc0.78.pdf'
plt.savefig(pltFile)

##################################################
# 4. Extend the training data the training data with the new data
##################################################

# read the extended line from label_mappings file:
label_mapping = runFolder + 'run02/label_mappings_trainSelect.csv' # lines = f.readlines()[296:360]
label_mapping = runFolder + 'run03/label_mappings_Jan31_2023.csv'  # lines = f.readlines()[296:383]
IDlist_labels = {}
extendedSampleNames = []
with open(label_mapping, "r") as f:
    lines = f.readlines()[296:385]
    for line in lines:
        if '?' in line:
            continue
        else:
            elements = line.strip().split(',')
            clusterID = elements[0].split('_')[1] + '___' + elements[1]
            label = elements[2]
            IDlist_labels[clusterID] = label
            extendedSampleNames.append(elements[0])
    

# extend the classifier
model_features_extended = model_features
model_labels_extended = list(model_labels)

for item in IDlist_labels.keys():
    print(item)
    theIndex = agg_cluster_list.index(item)
    print(theIndex)
    this_features = allFeatureAgg_mat[theIndex,]
    this_features_reshape = np.reshape(this_features, (1,16))
    model_features_extended = np.append(model_features_extended, this_features_reshape, axis=0)
    model_labels_extended.append(IDlist_labels[item])
    
model_labels_extended_array = numpy.array(model_labels_extended)

train = {'features': model_features_extended, 'labels':model_labels_extended}

training_data_file = runFolder + "run03/training_data.pkl"
with gzip.open(training_data_file, "wb") as f:
    pickle.dump(train, f)

supp_table05 = dataFolder + 'training_data.tsv'

with open(supp_table05, 'w') as f:
    lineString = ''
    for i in range(15):
        lineString = lineString + feature_names[i] + '\t'

    
    lineString = lineString + feature_names[15] + '\tlabel' + '\n'
    f.write(lineString)
    
    for i in range(300):
        lineString = ''
        for j in range(15):
            lineString = lineString + str(model_features_extended[i,j]) + '\t'

        lineString = lineString + str(model_features_extended[i, 15]) + '\t' + model_labels_extended[i] + '\n'
        f.write(lineString)

# 4.0.1 saving the training sampleID and training labels 
########################################
# TODO: save in one file 
sampleNames_list = list(sampleNames)
extendedSampleNames
allSampleNames = sampleNames_list + extendedSampleNames
model_labels_extended



# TODO: make a plot or something about count of the labels etc - and table of label by samples 

    
# 4.1 leave one out cross validation
########################################

reg = .020
test_labels = []
for i in range(len(model_labels_extended_array)):
    print(i)
    train_features = np.delete(model_features_extended, i, axis = 0)
    train_labels = np.delete(model_labels_extended_array, i, axis = 0)
    test_features = model_features_extended[i,]

    model = make_model(reg)
    model.fit(train_features, train_labels)
    label = model.predict(test_features.reshape(1,-1))
    test_labels.append(label)

segwayLabels = ['Bivalent', 'ConstitutiveHet', 'EnhancerLow', 'FacultativeHet', 'Transcribed', 'Unclassified', 'CTCF', 'Enhancer', 'Promoter', 'K9K36', 'Quiescent', 'PromoterFlanking']    
#segwayLabels_order = [7,2 ,8, 11, 4, 9, 6, 0, 1, 3, 10]
segwayLabels_order = [8, 11, 7, 2, 0, 6, 4, 9, 3, 1, 10]
segwayLabels_reorder = [segwayLabels[i] for i in segwayLabels_order]
Accuracy = accuracy_score(model_labels_extended, test_labels)
print(Accuracy)
from sklearn.metrics import confusion_matrix
cmat = confusion_matrix(model_labels_extended, test_labels, labels=segwayLabels_reorder)

from sklearn.metrics import ConfusionMatrixDisplay
disp = ConfusionMatrixDisplay(cmat, display_labels = segwayLabels_reorder)
disp.plot()
plt.rcParams.update({'font.size': 10})
plt.xticks(fontsize=10, rotation=90)
plt.yticks(fontsize=10)
plt.xlabel('predicted', fontsize=12)
plt.xlabel('predicted', fontsize=12)
plt.ylabel('true', fontsize=12)
plt.grid(False)
plt.tight_layout()
plt.show()
plot_file = runFolder + 'run03/confusion_mat_loo300_auc.72_reg.2.pdf' 
plt.savefig(plot_file)
plt.close('all')

##################################################
# 4.2 train the model with the whole training data and save the model
##################################################
reg = 0.020
model = make_model(reg)
model.fit(model_features_extended, model_labels_extended_array)
test_labels = model.predict(model_features_extended)
Accuracy = accuracy_score(model_labels_extended, test_labels)
print(Accuracy)

model_file = runFolder + "run03/model_300_reg.020_auc0.89.pickleV04.gz"
with gzip.open(model_file, "wb") as f:
    pickle.dump(model, f)


model_file = runFolder + "run03/model_300_reg.020_auc0.89V04.pickle.gz"
with gzip.open(model_file, "rb") as f:
    model = pickle.load(f)


#########################################
# 5. Apply the model on the 105 samples 
######################################### 

# loading the clusterer indices for the large heatmap plot
inputFileName = 'clusterings_dendro_indices.pkl'
inputFile = dataFolder + dataSubFolder + inputFileName
with open(inputFile, 'rb') as f:
    clusterings = pickle.load(f)

# the data we need for clustering
inputFileName = 'allData_mat_105run.pkl'
inputFile = dataFolder + dataSubFolder + inputFileName
with open(inputFile, 'rb') as f:
    data_105 = pickle.load(f)

# test print
print(data_105['cluster_list'][0:10])
print(data_105.keys())
print(clusterings.keys())

# load the model
# inputFile = runFolder + 'model_296exp_reg0.067_auc0.77on.32test.pickle'
model_file = runFolder + "run03/model_298_reg.020_auc0.90.pickle.gz"
model_file = runFolder + "run03/model_300_reg.020_auc0.89.pickleV04.gz"
with gzip.open(model_file, "rb") as f:
    the_model = pickle.load(f)

# 5.0 Get new labels
########################################
allFeatureAgg_mat  # features - see section 2 of this code
all_labels = the_model.predict(allFeatureAgg_mat) # classifier output
probs = model.predict_proba(allFeatureAgg_mat)
agg_cluster_list # cluster id

# 5.1 The clustering heatmap plot
########################################


'''
# save the labels
classifier_output = {'cluster_list': agg_cluster_list, 'interpretation_term': all_labels}
classifier_output_file = runFolder + "model_296exp_reg0.067_auc0.77on.32test_classifier_labels.pickle"
with open(classifier_output_file, "wb") as f:
    pickle.dump(classifier_output, f)
'''

# this is for the classifier data - the data that we used in this code. Not the data that we brought in for plotting. For the new plot, we only need the new assigned predicted interpretation terms, instead of the old ones.
feature_order = ['(09) initial exon',"(08) 5' flanking (1-1000 bp)", '(10) initial intron', '(11) internal exons','(01) H3K9me3',  '(02) H3K27me3', '(04) H3K4me3', "(16) 3' flanking (1000-10000 bp)", '(12) internal introns', '(03) H3K36me3', '(13) terminal exon', '(06) H3K4me1', '(14) terminal intron', "(07) 5' flanking (1000-10000 bp)", '(05) H3K27ac', "(15) 3' flanking (1-1000 bp)", ]

plot_data = data_105['whole_mat'][clusterings['ward_16'],0:16]
# plot_data = data_105['whole_mat'][0:17, 17:20] # FOR THE POSTER
feature_list_original = data_105['feature_list']
feature_list_reorder = ['initial exon', "5' flanking (1-1000 bp)",'initial intron', 'internal exons', 'internal introns',  'terminal exon', 'terminal intron', "3' flanking (1-1000 bp)", "5' flanking (1000-10000 bp)", "3' flanking (1000-10000 bp)", 'H3K4me3', 'H3K27ac', 'H3K4me1', 'H3K36me3', 'H3K27me3',  'H3K9me3']

book = data_105['cluster_list']

plot_original_labels = [book[i] for i in clusterings['ward_16']] # list of original labels
agg_cluster_list # list of cluster IDs
all_labels # list of predictions, it has the same order of the agg_cluster_list

'''
colormapping = {'Quiescent':[255,255,255], 'Promoter':[255,0,0], 'RegPermissive':[255,255,0], 'LowConfidence':[128,128,128], 'FacultativeHet':[128,0,128], 'Enhancer':[255,195,77], 'Bivalent':[200,200,100], 'Transcribed':[0,128,0], 'ConstitutiveHet':[137,236,218]}
'''

colormapping = {'Quiescent':[255,255,255], 'Promoter':[255,0,0], 'PromoterFlanking':[255,69, 0],'CTCF':[138,236,208], 'K9K36':[102,205,170], 'FacultativeHet':[128,0,128], 'Enhancer':[255,195,77], 'EnhancerLow':[255, 255, 0], 'Bivalent':[200,200,100], 'Transcribed':[0,128,0], 'ConstitutiveHet':[138,145,208], 'badSample1': [0,0,0], 'badSample2': [70,70,70]}

# list of bad samples based on the mean_low from prob_filter.py, 105 run samples
badsampleList1 = [ 86,  20,  22,  69,  57,   9,  92,  37,  54,  39, 103,  82,  23, 27,  25]
#badsampleList2 = [35, 40, 97, 15, 28]

# getting the predicted labels and color for the heatmap plot
plot_predicted_labels = []
plot_colors = []
plot_labels = []
p_index_list = []
for o_label in plot_original_labels:
    tokens = o_label.split('_')
    if int(tokens[0]) in badsampleList1:
        term = tokens[-1]
        id = tokens[0] + '___' + tokens[3]
        predict_label_index = agg_cluster_list.index(id)
        p_label = all_labels[predict_label_index]
        p_index_list.append(predict_label_index)
        plot_labels.append(p_label)
        if p_label == 'Unclassified':
            p_label = 'LowConfidence'
        new_label = id + '_' + p_label
        plot_predicted_labels.append(new_label)
        print('bs1')
        color = np.array(colormapping['badSample1'])/255
        plot_colors.append(color)
        continue

    #if int(tokens[0]) in badsampleList2:
    #    print('bs2')
    #    color = np.array(colormapping['badSample2'])/255
    #    plot_colors.append(color)
    #    continue
        
    term = tokens[-1]
    id = tokens[0] + '___' + tokens[3]
    predict_label_index = agg_cluster_list.index(id)
    p_label = all_labels[predict_label_index]
    p_index_list.append(predict_label_index)
    plot_labels.append(p_label)
    if p_label == 'Unclassified':
        p_label = 'LowConfidence'
    new_label = id + '_' + p_label
    plot_predicted_labels.append(new_label)
    color = np.array(colormapping[p_label])/255
    plot_colors.append(color)


# The following normalizaions are just for the heatmap plot
from scipy.stats import zscore
plot_data_z = zscore(plot_data, axis = 0)
plot_data_z_thr = np.where(plot_data_z > 4, 4, plot_data_z)

df = pd.DataFrame(plot_data_z_thr, index = plot_predicted_labels, columns = data_105['feature_list'])
df = df.apply(zscore, axis=0)
df = df[feature_list_reorder]

# the black lines are not zero, they are just very close to zero
sib = sns.clustermap(df, figsize=(6, 150), row_colors= plot_colors, row_cluster=False, col_cluster=False, dendrogram_ratio = [.25,.01], cbar_pos=(.02, .985,.03,.01))
sns.set(font_scale = .4)

figFile = '%sclustering_example_model_trained_ward_16_zcol_thr__v04.pdf' %(plotFolder)
print(figFile)
plt.savefig(figFile)
plt.close('all')

# just testing
for label in segwayLabels:
    print(plot_labels.count(label))
 
plot_original_labels_noID = []
for label in plot_original_labels:
    stripped = label.split('_')[-1]
    plot_original_labels_noID.append(stripped)
    
for label in segwayLabels:
    print(plot_original_labels_noID.count(label))


# doing the plot for the poster for a conf
plot_original_labels = book[0:100] # FOR THE POSTER
df = pd.DataFrame(plot_data_z_thr, columns = data_105['feature_list']) # FOR POSTER
df = pd.DataFrame[plot_data] # FOR POSTER
df = df.apply(zscore, axis=0)
df = df[feature_list_reorder]

sib = sns.clustermap(df, figsize=(6, 15), row_cluster=False, col_cluster=False, dendrogram_ratio = [.25,.01], cbar_pos=(.02, .985,.03,.01)) # FOR POSTER

figFile = figFile = '/Users/marjanfarahbod/Documents/talks/RSGDREAM2022/fig5.pdf'
print(figFile)
plt.savefig(figFile)
plt.close('all')


# 5.2 The probabilities heatmap plot
########################################
plot_probs = probs[p_index_list,0:-1]
df = pd.DataFrame(plot_probs, index=plot_predicted_labels, columns=model.classes_[0:-1])

fig, axs = plt.subplots(1, 1, figsize=[4, 150])

sns.heatmap(df)
plt.tight_layout()

figFile = runFolder + 'run02/probsplot.pdf'
print(figFile)
plt.savefig(figFile)

########################################
# 100. new classifier data for Abe
########################################

newLabel_data ={}
newLabel_data['labels'] = all_labels
newLabel_data['cluster_list'] = agg_cluster_list

ann_info = data_105['key']
index_accession = {}
for item in ann_info:
    index_accession[item['index']] = item['accession']
    
newLabel_data['index_accession']  = index_accession
#newLabel_data['annotation_info'] = data_105['key']

label_info = runFolder + "run02/labelDataForAbe.pickle"
with open(label_info, 'wb') as f:
    pickle.dump(newLabel_data, f)


##################################################
# OBSOLETE
##################################################

# getting the confusion matrix myself
confusion_mat = np.zeros([9,9], dtype=np.intc)
for i, label in enumerate(plot_original_labels_noID):
    x = segwayLabels.index(label)
    new_label = plot_labels[i]
    if new_label == 'Unclassified':
        y = 8
    else:
        y = segwayLabels.index(new_label)
    confusion_mat[x, y]+=1
        

# 4.2 3 fold cross validation
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# fix the cross_validation sets
cross_count = 3
index_sets = []
total_sampleCount = len(model_labels_extended_array)
select_count = int((1/cross_count)*total_sampleCount)
select_index = set(range(total_sampleCount))
for i in range(cross_count):
    print(i)
    one_set = set(random.sample(select_index, select_count))
    index_sets.append(one_set)
    select_index = select_index - one_set


# for a model, do the set
model_list = []
for i in range(cross_count):

    test_index = index_sets[i]
    train_index = set()
    for j in set(range(cross_count))-{i}: train_index = train_index.union(index_sets[j])

    
    train_index = list(train_index)
    test_index = list(test_index)

    train_features = model_features_extended[train_index,]
    test_features = model_features_extended[test_index,]
    train_labels = model_labels_extended_array[train_index]
    test_labels = model_labels_extended_array[test_index]

    reg = .02
    model = make_model(reg)
    model.fit(train_features, train_labels)

    Accuracy = accuracy_score(train_labels, model.predict(train_features))
    print(Accuracy)

    plot_confusion_matrix(model, train_features, train_labels, xticks_rotation='vertical')
    plt.grid(False)
    plt.tight_layout()
    plot_file = runFolder + 'run02/confusion_mat_train_3fold__%d.pdf' %(i)
    plt.savefig(plot_file)
    plt.close('all')

    Accuracy = accuracy_score(test_labels, model.predict(test_features))
    print(Accuracy)

    plot_confusion_matrix(model, test_features, test_labels, xticks_rotation='vertical')
    plt.grid(False)
    plt.tight_layout()
    plot_file = runFolder + 'run02/confusion_mat_test_3fold__%d.pdf' %(i)
    plt.savefig(plot_file)
    plt.close('all')
    
    model_list.append(model)

# doing the whole data train
train_index = range(total_sampleCount)
train_features = model_features_extended[train_index,]
train_labels = model_labels_extended_array[train_index]

reg = .027
model = make_model(reg)
model.fit(train_features, train_labels)
Accuracy = accuracy_score(train_labels, model.predict(train_features))
print(Accuracy)
plot_confusion_matrix(model, train_features, train_labels, xticks_rotation='vertical')
plt.grid(False)
plt.tight_layout()
plot_file = runFolder + 'run02/model_278exp_reg.027_auc.84_allData.pdf' 
plt.savefig(plot_file)
plt.show()

'''
What I did: I did the model with the 3 fold cross validation. With this setting, I got the accuracy .90+ for  train set, and ~70 for the test set. When I trained with the whole data, I got the accuracy around .84. My conclusion for now is that the variation is not enough in the cross validaion set, and the model got overfit. While with the complete set, there is less overfitting. I got confusion matrices and the rest of the stuff to make my report. The good thing is that model seems to be able to pick on the CTCF, k9k36 and all the newly added groups. at this point I should examine the big plot 

'''

