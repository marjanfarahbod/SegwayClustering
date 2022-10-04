#### THIS CODE IS HEAVILY IN PROGRESS ###
### NEEDS CLEAN UP AND DOCUMENTATION ####
# get the util file - done
# load the training data - done
# get the code for setting up the classifier and training it - done
# did the data use any normalization - no, but the mean varies a lot - can we do some form of normalization? - done. the mean is alright.
# get the new data. Let's reach all of them to 45, for the first round the unclassified is fine, or let them have wrong label.
# train the classifier and get the model - ...
# plot the two datasets and pick what you want to do for normalization
#
# add new samples
# train the new model

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
import pandas
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

classifier_tab_fname = dataFolder + 'classifier_data.tab'
label_mapping = dataFolder + 'label_mappings.txt'


#########################################
# Get bio label conversions
#########################################

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

#######################################################
# Convert classifier data to numpy matrix
#######################################################
# checkpoint: are the labels transferred correctly? 

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


example_bio_labels = numpy.array(example_bio_labels) # classifier label
feature_names = numpy.array(classifier_data_frame.drop("orig_label", 1).drop("concatenation_key",1).columns)
example_orig_labels = numpy.array(classifier_data_frame.orig_label) # we don't need to use this

segwayLabels = ['Quiescent', 'ConstitutiveHet', 'FacultativeHet', 'Transcribed', 'Promoter', 'Enhancer', 'RegPermissive', 'Bivalent', 'LowConfidence']

# note: the 'LowConfidence' up there, is called 'Unclassified' in the labels
label_training_counts_223 = np.zeros([9,1])
for label in example_bio_labels:
    if label == 'Unclassified':
        label_training_counts_223[8] +=1
    else:
        ind = segwayLabels.index(label)
        label_training_counts_223[ind] += 1

print(label_training_counts_223)

print(classifier_data_frame.mean(axis = 0))

################# Hard Coding Feature Order ####################
feature_names = ['(09) initial exon', '(01) H3K9me3', '(10) initial intron', '(02) H3K27me3', '(11) internal exons', '(04) H3K4me3', "(16) 3' flanking (1000-10000 bp)", '(12) internal introns', '(03) H3K36me3', '(13) terminal exon', '(06) H3K4me1', '(14) terminal intron', "(07) 5' flanking (1000-10000 bp)", '(05) H3K27ac', "(15) 3' flanking (1-1000 bp)", "(08) 5' flanking (1-1000 bp)"]

example_features = features_frame_to_matrix(classifier_data_frame, feature_names)


# getting data from the 105 sample runs
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

# plotting a set of samples just to see how different it is from the old classifier data
df = pd.DataFrame(allFeatureAgg_mat[0:100,])
classifier_data_frame # the other one

    fig, axs = plt.subplots(1, 2, figsize=[8, 10])

    # plotting the transcript data
    cmap = sns.diverging_palette(240, 10, s=80, l=30, as_cmap=True)
    sns.heatmap(df,  ax=axs[0])
    #sns.heatmap(df, center = 0, cmap=cmap, vmin=-2, vmax=2, ax=axs[0], cbar_kws=dict(use_gridspec=False, location="top"))

    # todo: plot the rest of the features
    # get the features:
    # get the index
    df2 = classifier_data_frame
    df3 = df2.drop(['orig_label', 'concatenation_key'], axis = 1)
    sns.heatmap(df3, ax=axs[1])
    #sns.heatmap(df, ax=axs[1], cbar_kws=dict(use_gridspec=False, location="top"))

    plt.show()

    plt.close('all')


# train the classifier- test it and plot the new labels. for the old label.

feature_names = ['(09) initial exon', '(01) H3K9me3', '(10) initial intron', '(02) H3K27me3', '(11) internal exons', '(04) H3K4me3', "(16) 3' flanking (1000-10000 bp)", '(12) internal introns', '(03) H3K36me3', '(13) terminal exon', '(06) H3K4me1', '(14) terminal intron', "(07) 5' flanking (1000-10000 bp)", '(05) H3K27ac', "(15) 3' flanking (1-1000 bp)", "(08) 5' flanking (1-1000 bp)"]

# converting the classifier_data_frame to the new feature order (we only do this to keep the order the same)
example_features = features_frame_to_matrix(classifier_data_frame, feature_names)

model_labels = example_bio_labels
model_features = example_features

def make_model(reg): 
    return RandomForestClassifier(min_samples_leaf=reg, criterion="entropy")

#reg = 1e-2
reg = .067 # it seems that reg 15 or 15/223 is the best one. I will balance the samples and train the model then. 
model = make_model(reg)
model.fit(model_features, model_labels)


model_file = runFolder + "model_223exp_reg15_auc0.75.pickle.gz"
with gzip.open(model_file, "w") as f:
    pickle.dump(model, f)


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
from sklearn.metrics import confusion_matrix
confusion_mat = confusion_matrix(model_labels, model.predict(model_features))
print(confusion_mat)
import matplotlib.pyplot as plt
from sklearn.metrics import plot_confusion_matrix
plot_confusion_matrix(model, model_features, model_labels)
plt.show()
plt.savefig("Confusion.pdf")
    
# example assigned term
# example identified term

# TODO: add rows to the training data, train the new classifier

# let's get everything to 43
# enhancer
e_labels = ['48___14', '102___5', '41___13', '38___10', '19___10', '3___10', '87___0', '60___9', '78___7', '43___12', '101___5', '39___14', '53___2', '43___12', ]
# promoter
p_labels = ['48___10', '101___1', '102___8', '87___0', '39___15', '61___6', '79___10', '84___15', '104___2', '85___3', '47___13', '5___12', '77___3', '79___0', '54___6', '92___15', '57___8', '93___8', '20___12', '41___3']
# bivalent
b_labels = ['26___10', '67___7', '29___3', '77___6', '17___4', '99___12', '93___14', '14___2', '31___0', ]
# quiescent
q_labels = ['17___0', '82___14', '28___0', '57___5', '37___2', '39___9', '18___4', '31___14', '58___15', '93___7', '86___10', '38___11', '103___9', '82___1']
# constit het
c_labels = ['19___12', '63___10', '98___12', '5___0', '31___4', '71___4', '52___8', '97___14', '102___0', '77___2', '79___9', '104___10', '26___0', '70___4', '80___5', '58___2']
# facultative het
f_labels = []

added_examples = {'e':e_labels, 'p':p_labels, 'b':b_labels, 'q':q_labels, 'c':c_labels, 'f':f_labels}
added_examples_label = {'e':'Enhancer', 'p':'Promoter', 'b':'Bivalent', 'q':'Quiescent', 'c':'ConstitutiveHet', 'f':'FacultativeHet'}

#model_labels = example_bio_labels
#model_features = example_features

model_features_extended = model_features
model_labels_extended = list(model_labels)

#index_cluster = '%d___%s' %(index, cluster)
#agg_cluster_list.append(index_cluster)

allFeatureAgg_mat = allFeatureAgg_mat[0:afa_index,] # this has all the features to the order of feature_names (like

for item in added_examples:
    print(item)
    for label in added_examples[item]:
        print(label)
        theIndex = agg_cluster_list.index(label)
        this_features = allFeatureAgg_mat[theIndex,]
        this_features_reshape = np.reshape(this_features, (1,16))
        model_features_extended = np.append(model_features_extended, this_features_reshape, axis=0)
        model_labels_extended.append(added_examples_label[item])

# just look at the plot - see how samples are built
df = pd.DataFrame(model_features_extended)
fig, axs = plt.subplots(1, 2, figsize=[8, 10])
sns.heatmap(df,  ax=axs[0])
plt.show()

# ge the labels back to numpy array

model_labels_extended_array = numpy.array(model_labels_extended)

# fit and train the model
total_sampleCount = len(model_labels_extended_array)
train_index = set(random.sample(range(total_sampleCount), 200))
test_index = set(range(total_sampleCount)) - train_index

train_index = list(train_index)
test_index = list(test_index)

train_features = model_features_extended[train_index,]
test_features = model_features_extended[test_index,]
train_labels = model_labels_extended_array[train_index]
test_labels = model_labels_extended_array[test_index]

#reg = 1e-2
reg = .067 # it seems that reg 15 or 15/223 is the best one. I will balance the samples and train the model then. 
model = make_model(reg)
model.fit(train_features, train_labels)


model_file = runFolder + "model_296exp_reg0.067_auc0.77on.32test.pickle.gz"
with gzip.open(model_file, "w") as f:
    pickle.dump(model, f)


from sklearn.metrics import precision_recall_fscore_support as score
import numpy as np
precision, recall, fscore, support = score(model_labels, model.predict(model_features), labels=['Quiescent', 'ConstitutiveHet', 'FacultativeHet', 'Transcribed', 'Promoter', 'Enhancer', 'RegPermissive', 'Bivalent', 'Unclassified'])
print('precision: {}'.format(precision))
print('recall: {}'.format(recall))
print('fscore: {}'.format(fscore))
from sklearn.metrics import accuracy_score
Accuracy = accuracy_score(train_labels, model.predict(train_features))
print("Accuracy: {}".format(Accuracy))
Accuracy = accuracy_score(test_labels, model.predict(test_features))
print("Accuracy: {}".format(Accuracy))
print(list(model.predict(model_features)))
print(list(model_labels))
from sklearn.metrics import confusion_matrix
confusion_mat = confusion_matrix(train_labels, model.predict(train_features))
confusion_mat = confusion_matrix(test_labels, model.predict(test_features))
print(confusion_mat)
import matplotlib.pyplot as plt
from sklearn.metrics import plot_confusion_matrix
plot_confusion_matrix(model, train_features, train_labels)
plot_confusion_matrix(model, test_features, test_labels)
plt.show()
figFile = runFolder + 'model_296exp_reg0.067_auc0.77on.32test_trainConfustion.pdf'
figFile = runFolder + 'model_296exp_reg0.067_auc0.77on.32test_testConfustion.pdf'
plt.savefig(figFile)


####### predict the model for 

# todo: run it on all the data based on the classification. Let's just see what we get. 



