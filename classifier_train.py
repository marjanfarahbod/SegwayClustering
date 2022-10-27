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

# TODO: plot the confusion matrix - DONE
# TODO: get the probs from the classifier
# TODO: print the text file

# 0. Initiation and set up
# 1. Get the old classifier training data 
## 1.0 Get the old  label_mappings (this is also where we get the updated labels, but from the old samples)
## 1.1 Get the old featrues and fetch their labels based on the label_mappings
# 2. Getting data from the 105 sample runs, getting the feature matrix.
# ** plotting a set of samples just to see how different it is from the old classifier data
# 3. Training the classifier
# 4. Extend the training data
# 5. Apply the model on the whole data and get the large matrix with colors and all
# DRAFT

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

# 1.0 Getting the labels (this is also where we get the updated labels, but from the old samples)
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
#orig_labels_set = set()
with open(label_mapping, "r") as f:
    for i in range(N):
        line = next(f).strip().split(',')
        if line[0] == "concatenation_key":
            continue
        concatenation_key = line[0]
        orig_label = line[1]
        label = line[2]
        if not concatenation_key in label_mappings:
            label_mappings[concatenation_key] = {}
        label_mappings[concatenation_key][orig_label] = label
        #orig_labels_set.add(orig_label)


all_ref_bio_labels = set.union(*map(set, map(lambda x: x.values(), label_mappings.values())))

# 1.1 Getting the features
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

#########################################
# 2. Getting data from the 105 sample runs, getting the feature matrix.
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

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

########################################
# 3. Training the classifier 
########################################

# 3.0 Training with the old training data
########################################
 
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

########################################
# 4. Extend the training data - report the phase 2 of the classifier
########################################
# OUT OF CODE TODO: add samples to the training data - DONE

# read the extended line from label_mappings file:
label_mapping = runFolder + 'run02/label_mappings_trainSelect.csv'
IDlist_labels = {}
with open(label_mapping, "r") as f:
    lines = f.readlines()[296:360]
    for line in lines:
        if '?' in line:
            continue
        else:
            elements = line.strip().split(',')
            clusterID = elements[0].split('_')[1] + '___' + elements[1]
            label = elements[2]
            IDlist_labels[clusterID] = label
    

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


# DRAFT
########################################
# 1. load the original traning data - add their labels

# 2. load the list of extended - add the labels

# 3. document the new training data

# 4. document the train process with the parameters

# 5. document the model

Accuracy = accuracy_score(test_labels, model.predict(test_features))
print(Accuracy)
plot_confusion_matrix(model, test_features, test_labels, xticks_rotation='vertical')
plt.grid(False)
plt.tight_layout()

plt.show()

model_list = []
model = model_list[9]
for i in range(50):

    print(i)
    train_percentage = 2/3
    total_sampleCount = len(model_labels_extended_array)
    train_index = set(random.sample(range(total_sampleCount), int(train_percentage*total_sampleCount)))
    test_index = set(range(total_sampleCount)) - train_index

    train_index = list(train_index)
    test_index = list(test_index)

    train_features = model_features_extended[train_index,]
    test_features = model_features_extended[test_index,]
    train_labels = model_labels_extended_array[train_index]
    test_labels = model_labels_extended_array[test_index]

    #train_features = model_features_extended
    #train_labels = model_labels_extended_array

    reg = .02
    #model = make_model(reg)
    #model.fit(train_features, train_labels)
    #model_list.append(model)
    Accuracy = accuracy_score(train_labels, model.predict(train_features))
    print(Accuracy)
    trainA = Accuracy


    Accuracy = accuracy_score(test_labels, model.predict(test_features))
    print(Accuracy)
    testA = Accuracy

    print(trainA - testA)
    #model_list.append(model)
    


# let's get everything to 43
# enhancer
e_labels = ['48___14', '102___5', '41___13', '38___10', '19___10', '3___10', '87___0', '60___9', '78___7', '43___12', '101___5', '39___14', '53___2', '43___12', ]
# promoter
p_labels = ['48___10', '101___1', '102___8',  '39___15', '61___6', '79___10', '84___15', '104___2', '85___3', '47___13', '5___12', '77___3', '79___0', '54___6', '92___15', '57___8', '93___8', '20___12', '41___3']
# bivalent fn = 77-4, 14-3: many things that are bivalent are expressed, we can go either way. 
b_labels = ['26___10', '67___7', '29___3', '77___6', '17___4', '99___12', '93___14', '14___2', '31___0', ]
# quiescent
q_labels = ['17___0', '82___14', '28___0', '57___5', '37___2', '39___9', '18___4', '31___14', '58___15', '93___7', '86___10', '38___11', '103___9', '82___1']
# constit het
c_labels = ['19___12', '63___10', '98___12', '5___0', '31___4', '71___4', '52___8', '97___14', '102___0', '77___2', '79___9', '104___10', '26___0', '70___4', '80___5', '58___2']
# facultative het
f_labels = []
#
# broadk4me1
weak_enhancer = [] # shorter region, lower signal
k4me1_labels = [] # lower signal, larger region
promoter_flanking = []# it tends to be some k4me1, could be short.  a little lower strength than promoters - it is derived by segmentation process picking it as something
# how to organize them = make a folder , make an html. 
k36_k9 = []
# this is the two signal high 
k_labels = []

# TODO: pick the regpermissive, longer and weaker signal, and k-labels


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


total_sampleCount = len(model_labels_extended_array)
train_index = set(random.sample(range(total_sampleCount), 200))
test_index = set(range(total_sampleCount)) - train_index

train_index = list(train_index)
test_index = list(test_index)

train_features = model_features_extended[train_index,]
test_features = model_features_extended[test_index,]
train_labels = model_labels_extended_array[train_index]
test_labels = model_labels_extended_array[test_index]

#train_features = model_features_extended
#train_labels = model_labels_extended_array


#reg = 1e-2
reg = .055 # it seems that reg 15 or 15/223 is the best one. I will balance the samples and train the model then.
reg = .02
model = make_model(reg)
model.fit(train_features, train_labels)
Accuracy = accuracy_score(train_labels, model.predict(train_features))
plot_confusion_matrix(model, train_features, train_labels, xticks_rotation='vertical')
plt.grid(False)
plt.tight_layout()

plt.show()


Accuracy = accuracy_score(test_labels, model.predict(test_features))
plot_confusion_matrix(model, test_features, test_labels, xticks_rotation='vertical')
plt.grid(False)
plt.tight_layout()

plt.show()


model_features
model_labels
model.fit(model_features, model_labels)
Accuracy = accuracy_score(model_labels, model.predict(model_features))

Accuracy = accuracy_score(test_labels, model.predict(test_features))
print("Accuracy: {}".format(Accuracy))
plot_confusion_matrix(model, model_features, model_labels, xticks_rotation='vertical')

plot_confusion_matrix(model, test_features, test_labels, xticks_rotation='vertical')
plt.grid(False)
plt.tight_layout()

plt.show()
figFile = runFolder + 'model_223_reg0.067_auc0.74.pdf'
plt.savefig(figFile)

model_file = runFolder + "run02/model_278exp_reg0.028_auc0.86_wholeData.pickle.gz"
with gzip.open(model_file, "w") as f:
    pickle.dump(model, f)


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
pfrom sklearn.metrics import plot_confusion_matrix
plot_confusion_matrix(model, train_features, train_labels)
plot_confusion_matrix(model, test_features, test_labels)
plt.show()
figFile = runFolder + 'model_296exp_reg0.067_auc0.77on.32test_trainConfustion.pdf'
figFile = runFolder + 'model_296exp_reg0.067_auc0.77on.32test_testConfustion.pdf'
plt.savefig(figFile)

####### predict the model for

#########################################
# Predict the model and examine that
######################################### 


# stats on the previous labels

inputFileName = 'clusterings_dendro_indices.pkl'
inputFile = dataFolder + dataSubFolder + inputFileName
with open(inputFile, 'rb') as f:
    clusterings = pickle.load(f)

inputFileName = 'allData_mat_105run.pkl'
inputFile = dataFolder + dataSubFolder + inputFileName
with open(inputFile, 'rb') as f:
    data_105 = pickle.load(f)

print(data_105['cluster_list'][0:10])

# load the model
inputFile = runFolder + 'model_296exp_reg0.067_auc0.77on.32test.pickle'
with open(inputFile, 'rb') as f:
    the_model = pickle.load(f)

allFeatureAgg_mat  # features
all_labels = the_model.predict(allFeatureAgg_mat) # classifier output
agg_cluster_list # cluster id

classifier_output = {'cluster_list': agg_cluster_list, 'interpretation_term': all_labels}
classifier_output_file = runFolder + "model_296exp_reg0.067_auc0.77on.32test_classifier_labels.pickle"
with open(classifier_output_file, "wb") as f:
    pickle.dump(classifier_output, f)


# here, we only need the all_labels and the clusterIDs, we have the actual values normalized for the plot from the data_105.

print(data_105.keys())
print(clusterings.keys())

# this is for the classifier data - the data that we used in this code. Not the data that we brought in for plotting. For the new plot, we only need the new assigned predicted interpretation terms, instead of the old ones.
feature_order = ['(09) initial exon',"(08) 5' flanking (1-1000 bp)", '(10) initial intron', '(11) internal exons','(01) H3K9me3',  '(02) H3K27me3', '(04) H3K4me3', "(16) 3' flanking (1000-10000 bp)", '(12) internal introns', '(03) H3K36me3', '(13) terminal exon', '(06) H3K4me1', '(14) terminal intron', "(07) 5' flanking (1000-10000 bp)", '(05) H3K27ac', "(15) 3' flanking (1-1000 bp)", ]


plot_data = data_105['whole_mat'][clusterings['ward_16'],0:16]
feature_list_original = data_105['feature_list']
feature_list_reorder = ['initial exon', "5' flanking (1-1000 bp)",'initial intron', 'internal exons', 'internal introns',  'terminal exon', 'terminal intron', "3' flanking (1-1000 bp)", "5' flanking (1000-10000 bp)", "3' flanking (1000-10000 bp)", 'H3K4me3', 'H3K27ac', 'H3K4me1', 'H3K36me3', 'H3K27me3',  'H3K9me3']

book = data_105['cluster_list']
plot_original_labels = [book[i] for i in clusterings['ward_16']] # list of original labels
agg_cluster_list # list of cluster IDs
all_labels # list of predictions, it has the same order of the agg_cluster_list

colormapping = {'Quiescent':[255,255,255], 'Promoter':[255,0,0], 'RegPermissive':[255,255,0], 'LowConfidence':[128,128,128], 'FacultativeHet':[128,0,128], 'Enhancer':[255,195,77], 'Bivalent':[200,200,100], 'Transcribed':[0,128,0], 'ConstitutiveHet':[137,236,218]}

plot_predicted_labels = []
plot_colors = []
plot_labels = []
for o_label in plot_original_labels:
    tokens = o_label.split('_')
    term = tokens[-1]
    id = tokens[0] + '___' + tokens[3]
    predict_label_index = agg_cluster_list.index(id)
    p_label = all_labels[predict_label_index]
    plot_labels.append(p_label)
    if p_label == 'Unclassified':
        p_label = 'LowConfidence'
    new_label = id + '_' + p_label
    plot_predicted_labels.append(new_label)
    color = np.array(colormapping[p_label])/255
    plot_colors.append(color)

df = pd.DataFrame(data_105['color'][clusterings['ward_16'][700:1000],])
sns.heatmap(df)
plt.show()

for label in segwayLabels:
    print(plot_labels.count(label))
 
plot_original_labels_noID = []
for label in plot_original_labels:
    stripped = label.split('_')[-1]
    plot_original_labels_noID.append(stripped)
    
for label in segwayLabels:
    print(plot_original_labels_noID.count(label))

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
        

#from sklearn.metrics import confusion_matrix
#confusion_mat = confusion_matrix(plot_labels, plot_original_labels_noID)
fig, axs = plt.subplots(1,1, figsize = (6,6))

df = pd.DataFrame(confusion_mat, index = segwayLabels, columns = segwayLabels)
sns.heatmap(df, annot=True, fmt='d')
sns.set(font_scale=1)
plt.xlabel('classifier 01')
plt.ylabel('classifier 00')
plt.title('classifier label overlap')
plt.tight_layout()
figFile = plotFolder + 'classifier_label_overlap_00_01.pdf'
print(figFile)
plt.savefig(figFile)
plt.close('all')

plt.show()

### TODO: print the labels in a file
# plot_predicted_labels
# labels selected for training
# plot_original_labels

classifier_label_file = dataFolder + dataSubFolder + 'label_file.txt'
with open(classifier_label_file, 'w') as output:
    header = 'classifier_00_labels\tclassifier_01_train_select\tclassifier_01_labels\n'
    output.write(header)
    for i, olabel in enumerate(plot_original_labels):
        line =  '%s\t%s\t%s\n' %(olabel, ' ', plot_predicted_labels[i])
        output.write(line)


from scipy.stats import zscore

plot_data_z = stats.zscore(plot_data, axis = 0)
plot_data_z_thr = np.where(plot_data_z > 4, 4, plot_data_z)

#myMat = plot_data[0:20, 0:16]
#df = pd.DataFrame(plot_data, index = plot_original_labels[0:20], columns = data_105['feature_list'])
df = pd.DataFrame(plot_data_z_thr, index = plot_predicted_labels, columns = data_105['feature_list'])
df = df.apply(zscore, axis=0)
df = df[feature_list_reorder]

# the black lines are not zero, they are just very close to zero
sib = sns.clustermap(df, figsize=(6, 150), row_colors= plot_colors, row_cluster=False, col_cluster=False, dendrogram_ratio = [.25,.01], cbar_pos=(.02, .985,.03,.01))
sns.set(font_scale = .4)


figFile = '%sclustering_example_model_trained_ward_16_zcol_thr.pdf' %(plotFolder)
print(figFile)
plt.savefig(figFile)
plt.close('all')

# reordering dataframe columns: 
sns.heatmap(df)
plt.show()
# the black lines are not zero, they are just very close to zero


# getting count of diferent labels, from the original to the new one. 



# todo: run it on all the data based on the classification. Let's just see what we get. 


# todo: make a tab delimitered file, one column the old labels, one column, new labels. The other column 
