# At this stage, we are not sure what classifier is doing and with a classifier we have the limitation of training and training data. I want to know how do the samples look with regards to the values that go into the classifier. I also want to compare them to the transcript data report.

# TODO:
# 1. Do the clustering of the Segway labels using the Classifier input features, all of them.
# 2. Add the Expression data to the plots, how does it look? How do you expect it to look?
# 3. Add the remaining track data (those not used in the classification) to the plot, how do they look?

# The outcome:
# 1. To identify label groupings and examine their transcript data
# 2. To pick training samples for the classifier for the under-represented labels

# code sections
# 1. Get the meta info for the annotations
# 2. Get the track values for labels
# 3. Getting the feature aggregation data - I am using the previous code
# 4. Get the transcript data for labels
# 5. Put together all the matrices
# 6. Make the clustering plots
# 7. Get the transcript plots based on the clusterings
# 8. Examinatation of individual samples/labels
# 9. Notes and drafts

# 1. Get the meta info for the annotations
########################################

'''
For a batch of annotations, we need these sets of meta info:
1. A dictionary of class "annotation" - annMeta, with ENCODE accession as key
2. A dictionary of tissue info with name and uberon tissue id - tissue_info, with ENCODE accession as key
3. List of ENCODE accessions for samples: sampleFolder_list
4. List of histone tracks that went into the classifier
5. List of Segway interpretation terms/biological labels
'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
from util import *

# for all the samples, extract their value and cluster label and the sample ID

# mapping Segway accession to a number - since the accession is too big or long, this I can do with that accession tissue thing.

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
dataSubFolder = 'testBatch105/fromAPI/'

inputFile = dataFolder + dataSubFolder + 'metaInfo.pkl'
with open(inputFile, 'rb') as f:
    annMeta = pickle.load(f)

inputFile = dataFolder + dataSubFolder + 'biosample_tissue_info.pkl'
with open(inputFile, 'rb') as f:
    tissue_info = pickle.load( f)

sampleFolder_list = list(tissue_info.keys())

# list of tracks that are used for classification
classifier_tracks = ['H3K36me3', 'H3K27me3', 'H3K4me3', 'H3K9me3', 'H3K27ac', 'H3K4me1']

# list of segway labels to the order that I like
segwayLabels = ['Quiescent', 'ConstitutiveHet', 'FacultativeHet', 'Transcribed', 'Promoter', 'Enhancer', 'RegPermissive', 'Bivalent', 'LowConfidence']
# build the whole mat:

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
# 1.1 Gathering data for each annotation 
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
ann_info_list = []
ann_info = {}
c = 0 
for sampleFolder in sampleFolder_list:
    print(sampleFolder)

    info = {} # keeps the info for the annotations - it also has an index which is just 0-104 for the 105 runs
    info['index'] = c
    c = c+1
    # load the overlap dataframe
    sampleFolderAdd = dataFolder + dataSubFolder + sampleFolder + '/'
    mapping_file = sampleFolderAdd + 'trackname_assay.txt'
    signal_file = sampleFolderAdd + 'signal_distribution.tab.txt'

    originalBedFile = annMeta[sampleFolder]['bedFile']
    originalBed_accession = originalBedFile.split('.')[0]

    segway_ann_sum_file = dataFolder + dataSubFolder + sampleFolder + '/' + originalBed_accession + '_annotationSummary.pkl'

    # get the segway annotation summary file
    with open(segway_ann_sum_file, 'rb') as pickledFile:
        segway_anns = pickle.load(pickledFile)

    info['segway_anns'] = segway_anns
    
    segway_cluster_list = list(segway_anns['clusters'].keys())
    # sorting the cluster labels for segway
    sortedClusterList = []
    for label in segwayLabels:
        #print(label)
        for item in segway_cluster_list:
            #print(item)
            if item.split('_')[1] == label:
                sortedClusterList.append(item)

    segway_cluster_list = sortedClusterList

    # read the mapping file
    track_assay_map = {}
    inputTrack_list = []
    with open(mapping_file) as inputFile:
        print('I am here')
        for line in inputFile:
            fields = line.strip().split()
            track_assay_map[fields[0]] = fields[1]
            if (fields[1] in inputTrack_list):
                print('here we are')
            inputTrack_list.append(fields[1])

    info['track_assay_map'] = track_assay_map
    info['inputTrack_list'] = inputTrack_list

    # getting the segway input track count
    segway_track_count = 0
    with open(signal_file, 'r') as inputFile:
        header = inputFile.readline().split()[0]
        previousCluster_int = inputFile.readline().split()[0]
        cluster_int = previousCluster_int
        while previousCluster_int == cluster_int:
            previousCluster_int = cluster_int
            cluster_int = inputFile.readline().split()[0]
            segway_track_count += 1

    print(segway_track_count)
    print(len(inputTrack_list))
    if not(segway_track_count==len(inputTrack_list)):
        print('not equal')

    # the above two are the same for all samples, meaning that tracks in singal_distribution and track_assay_map are the same

    info['segway_track_count'] = segway_track_count

    ann_info[sampleFolder] = info
    info['accession'] = sampleFolder
    ann_info_list.append(info)

outputFileName = 'all_annInfo.pkl'
outputFile = dataFolder + dataSubFolder + outputFileName
with open(outputFile, 'wb') as f:
    pickle.dump(ann_info, f)

inputFileName = 'all_annInfo.pkl'
inputFile = dataFolder + dataSubFolder + inputFileName
with open(inputFile, 'rb') as f:
    ann_info = pickle.load(f)

outputFileName = 'all_annInfo_list.pkl'
outputFile = dataFolder + dataSubFolder + outputFileName
with open(outputFile, 'wb') as f:
    pickle.dump(ann_info_list, f)

inputFileName = 'all_annInfo_list.pkl'
inputFile = dataFolder + dataSubFolder + inputFileName
with open(inputFile, 'rb') as f:
    ann_info_list = pickle.load(f)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
# 1.2 Some stat from tracks per sample
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
# get the total different tracks in the data: 11
# get the track count distriubtion for samples
# get the count of samples for each track

inputFileName = 'all_annInfo_list.pkl'
inputFile = dataFolder + dataSubFolder + outputFileName
with open(inputFile, 'rb') as f:
    ann_info_list = pickle.load(f)

unique_track_list = []
for ann in ann_info_list:
    ann_track_list = ann['track_assay_map'].values()
    for track in ann_track_list:
        if not(track in unique_track_list):
            unique_track_list.append(track)


track_counts = np.zeros(105)
for i, ann in enumerate(ann_info_list):
    track_counts[i] = ann['segway_track_count']

sample_count_per_track = np.zeros(len(unique_track_list))
for ann in ann_info_list:
    ann_track_list = ann['track_assay_map'].values()
    #book = np.zeros((len(unique_track_list)))
    for track in ann_track_list:
        i = unique_track_list.index(track)
        #book[i] +=1
        sample_count_per_track[i] += 1

    print(book)

fix, axs = plt.subplots(1, figsize=(6,4))
plt.grid()
plt.bar(range(11), sample_count_per_track)
plt.ylim([15,110])
plt.xticks(range(11))
axs.set_xticklabels(unique_track_list, rotation=90)
plt.tight_layout()
plt.show()
plotFolder_add = dataFolder + dataSubFolder
figFile = plotFolder_add + 'sample_count_for_tracks.pdf'
print(figFile)
plt.savefig(figFile)
plt.close('all')

########################################
# 2. Get the track values for labels
########################################

# build the whole matrix of cluster/track values, with the index label of the clusters
# first get the segway tracks, then CTCF, DNase, ATAC-seq, EPS is last

segway_track_ordered = classifier_tracks
segway_track_ordered.append('CTCF')
segway_track_ordered.append('DNase-seq')
segway_track_ordered.append('ATAC-seq')
segway_track_ordered.append('POLR2A')
segway_track_ordered.append('EP300')

allSignalDist_mat = np.zeros((105*16, len(segway_track_ordered))) * np.nan # this one fills with index
cluster_classCode = np.zeros(105*16)
cluster_color = np.zeros((105*16, 3))
asd_index = 0 # allSignalDist_mat row index
allClusters_list = [] # this one just append
print(allSignalDist_mat.shape)

for ann in ann_info_list:

    index = ann['index']
    print(index)
    sampleFolder =  ann['accession']

    sampleFolderAdd = dataFolder + dataSubFolder + sampleFolder + '/'
    signal_file = sampleFolderAdd + 'signal_distribution.tab.txt'

    # getting the segway input track count
    segway_track_count = len(ann['track_assay_map'])


    segway_cluster_list = list(ann['segway_anns']['clusters'].keys())
    sortedClusterList = []
    for label in segwayLabels:
        #print(label)
        for item in segway_cluster_list:
            #print(item)
            if item.split('_')[1] == label:
                sortedClusterList.append(item)

    segway_cluster_list = sortedClusterList

    if len(segway_cluster_list)<16:
        print('here')
        print(len(segway_cluster_list))

    
    cluster_class_map = {}
    for cluster_label in segway_cluster_list:
        int_label = cluster_label.split('_')[0]
        cluster_class_map[int_label] = cluster_label


    track_assay_map = ann['track_assay_map']
    signal_dist_mat = np.zeros((len(segway_cluster_list), len(segway_track_ordered))) * np.nan
    with open(signal_file, 'r') as inputFile:
        header = inputFile.readline()
        for line in inputFile:
            fields = line.strip().split()
            track = track_assay_map[fields[1]]
            track_ind = segway_track_ordered.index(track) # track list order
            segway_cluster = cluster_class_map[fields[0]]
            cluster_ind = segway_cluster_list.index(segway_cluster) # cluster list order
            signal_dist_mat[cluster_ind][track_ind] = round(float(fields[2]), 4)

            
    # do it or not?
    # z_signal_dist_mat = stats.zscore(signal_dist_mat, axis = 0)

    # fill up the whole cluster list and the signal matrix
    for cluster in segway_cluster_list:
        index_cluster = '%d___%s' %(index, cluster)
        allClusters_list.append(index_cluster)
        #class_code = segwayLabels.index(cluster.split('_')[1])
        #cluster_classCode[asd_index] = class_code
        #color = list(map(int, ann['segway_anns']['clusters'][cluster].color.split(',')))
        #cluster_color[asd_index,:] = np.array(color)
        

    for i in range(signal_dist_mat.shape[0]):
        print(i)

        #if math.isnan(signal_dist_mat[i, 0]):
        #   print('ind: %d' %(asd_index))
        allSignalDist_mat[asd_index,:] = signal_dist_mat[i,:]
        cluster = segway_cluster_list[i]
        class_code = segwayLabels.index(cluster.split('_')[1])
        cluster_classCode[asd_index] = class_code
        color = list(map(int, ann['segway_anns']['clusters'][cluster].color.split(',')))
        cluster_color[asd_index,:] = np.array(color)
        asd_index +=1

allSignalDist_mat = allSignalDist_mat[0:asd_index,:]

#########################################
# 3. Getting the feature aggregation data - I am using the previous code
#########################################

# the original 16 features that go to the classifier - just keeping it here
feature_names = ['(09) initial exon', '(01) H3K9me3', '(10) initial intron', '(02) H3K27me3', '(11) internal exons', '(04) H3K4me3', "(16) 3' flanking (1000-10000 bp)", '(12) internal introns', '(03) H3K36me3', '(13) terminal exon', '(06) H3K4me1', '(14) terminal intron', "(07) 5' flanking (1000-10000 bp)", '(05) H3K27ac', "(15) 3' flanking (1-1000 bp)", "(08) 5' flanking (1-1000 bp)"]

# getting the classifier features from the previous code:

import util

allFeatureAgg_mat = np.zeros((105*16, 10)) # this one fills with index
featureAgg_names = ["5' flanking (1-1000 bp)", "5' flanking (1000-10000 bp)", 'initial exon', 'initial intron', 'internal exons', 'internal introns', 'terminal exon', 'terminal intron', "3' flanking (1-1000 bp)", "3' flanking (1000-10000 bp)"]

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
            feature = ann_feature[5:] # stripping the first 4 characters
            value = ann_features[cluster][ann_feature]
            if feature in featureAgg_names:
                feature_index = featureAgg_names.index(feature)
                allFeatureAgg_mat[afa_index, feature_index] = value
        index_cluster = '%d___%s' %(index, cluster)
        agg_cluster_list.append(index_cluster)
        afa_index +=1 

allFeatureAgg_mat = allFeatureAgg_mat[0:afa_index,]

########################################
# 4. Get the transcript data for labels
########################################

# we have agg_cluster_list (from the allaggregation feature matrix) and allClusters_list (from the all signal distribution matrix). Now let's gather the transcribed matrix.

allTranscribe_mat = np.zeros((105*16, 160*2)) # this one fills with index
transcribed_cluster_list = []
transcribed_cluster_index = 0
transcribed_accession_list = []
transcribed_index_list = []

for ann in ann_info_list:

    # get the index
    annAccession = ann['accession']
    index = ann['index']
    
    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
    print(sampleFolderAdd)
    
    expFile = 'defaultExp_5kg_expSummary.pkl'
    inputFile = sampleFolderAdd + expFile

    try:
        with open(inputFile, 'rb') as f:
            expMat = pickle.load(f)
            transcribed_accession_list.append(annAccession)
            transcribed_index_list.append(index)
    except FileNotFoundError:
        print('no transcript data for this file')
        continue

    # getting the list of cluster list for the expression matrix
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


    clusterCount = len(clusterList)
    fractions = np.zeros([clusterCount])
    total_bp = 0

    for i, cluster in enumerate(sortedClusterList):
        print(clusters[cluster].bp_count)
        fractions[i] = clusters[cluster].bp_count
        total_bp += clusters[cluster].bp_count

    fractions = fractions / total_bp
    
    zeroMat = expMat['clusterMats'][0]
    # make it the ratio
    zeroMat = zeroMat / np.sum(zeroMat, axis=0)[np.newaxis,:]
    # versus expected
    zeroMat = zeroMat / fractions[:, None]
    zeroMat = np.log10(zeroMat, out=np.zeros_like(zeroMat), where=(zeroMat!=0))

    highMat = expMat['clusterMats'][2]
    # make it the ratio
    highMat = highMat / np.sum(highMat, axis=0)[np.newaxis,:]
    # versus expected
    highMat = highMat / fractions[:, None]
    highMat = np.log10(highMat, out=np.zeros_like(highMat), where=(highMat!=0))

        
    allTranscribe_mat[transcribed_cluster_index:transcribed_cluster_index+len(sortedClusterList),0:160] = zeroMat
    allTranscribe_mat[transcribed_cluster_index:transcribed_cluster_index+len(sortedClusterList),160:] = highMat
    
    transcribed_cluster_index += len(sortedClusterList)

    for cluster in sortedClusterList:
        index_cluster = '%d___%s' %(index, cluster)
        transcribed_cluster_list.append(index_cluster)

allTranscribe_mat = allTranscribe_mat[0:transcribed_cluster_index,]

#########################################
# 5. Putting  together all the matrices
#########################################

# let's get samples with transcription data only, and investigate them only
allTranscribe_mat # matrix of transcribed values
transcribed_cluster_list #

allFeatureAgg_mat
agg_cluster_list
# transfer the allFeatureAgg_mat to the range of allSignalDist_mat
omin = np.min(allFeatureAgg_mat)
omax = np.max(allFeatureAgg_mat)
tmin = np.nanmin(allSignalDist_mat)
tmax = np.nanmax(allSignalDist_mat)

t_allFeatureAgg_mat = (((allFeatureAgg_mat - omin)/(omax - omin)) * (tmax-tmin)) + tmin

allSignalDist_mat
allClusters_list

wholeMat_transcribed = np.zeros((len(transcribed_cluster_list), 16+(2*160)))
whole_cluster_class_color = np.zeros(cluster_color.shape)
extra_features_mat = np.zeros((len(transcribed_cluster_list), 5))

for cluster in transcribed_cluster_list:
    cluster_sig_index = allClusters_list.index(cluster)
    #print(cluster)
    #print(cluster_sig_index)
    #print(cluster_color[cluster_sig_index])
    
    tempSplit = cluster.split('_')
    temp_agg = tempSplit[0] + '___' + tempSplit[3]
    cluster_agg_index = agg_cluster_list.index(temp_agg)

    cluster_trans_index = transcribed_cluster_list.index(cluster)
    extra_features_mat[cluster_trans_index, :] = allSignalDist_mat[cluster_sig_index, 6:]
    
    whole_cluster_class_color[cluster_trans_index,] = cluster_color[cluster_sig_index,]

    wholeMat_transcribed[cluster_trans_index, 0:6] = allSignalDist_mat[cluster_sig_index, 0:6]
    wholeMat_transcribed[cluster_trans_index, 6:16] = t_allFeatureAgg_mat[cluster_agg_index, 0:10]
    wholeMat_transcribed[cluster_trans_index, 16:] = allTranscribe_mat[cluster_trans_index,]

# these three:
transcribed_cluster_list # list of clusters
wholeMat_transcribed # list of values for each cluster
whole_cluster_class_color # list of colors

########################################
# 6. Make the clustering plots
########################################

import seaborn as sns

plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/testBatch105/'

feature_list = segway_track_ordered[0:6]
for agg_feature in featureAgg_names:
    feature_list.append(agg_feature)

feature_count = [6,16] # either 6 or 16
linkage_method = ['average', 'ward'] # either 'average' or 'ward'
fig_width = [2,7] # either 6 or 12, based on the feature_count
pcc = wholeMat_transcribed.shape[0] #150 # plot cluster count
colors = whole_cluster_class_color[0:pcc,] / 255
 
fw = 6
dendro_indices = {}
for fc in feature_count:
    print('feature count: %d' %(fc))
    for method in linkage_method:
        print('method: %s' %(method))

        myMat = wholeMat_transcribed[0:pcc, 0:fc]
        df = pd.DataFrame(myMat, index = transcribed_cluster_list[0:pcc], columns = feature_list[0:fc])
        if fc == 16:
            fw = fig_width[1]

        # the black lines are not zero, they are just very close to zero
        sib = sns.clustermap(df, figsize=(fw, 150), row_colors= colors, method = method, dendrogram_ratio = [.25,.01], cbar_pos=(.02, .985,.03,.01))
        sns.set(font_scale = .4)

        clustering_name = '%s_%d' %(method, fc)
        dendro_indices[clustering_name] = sib.dendrogram_row.reordered_ind

        figFile = '%sclustering_%s.pdf' %(plotFolder, clustering_name)
        print(figFile)
        plt.savefig(figFile)
        plt.close('all')


outputFileName = 'clusterings_dendro_indices.pkl'
outputFile = dataFolder + dataSubFolder + outputFileName
with open(outputFile, 'wb') as f:
    pickle.dump(dendro_indices, f)

########################################
# 7. Get the transcript plots based on the clusterings
########################################

clusterings_list = list(dendro_indices.keys())[1:]

for clustering in clusterings_list:

    clustering_name = clustering

    indices = dendro_indices[clustering_name]

    myMat = wholeMat_transcribed[indices, 17:]
    label_list = [transcribed_cluster_list[i] for i in indices]

    fig, axs = plt.subplots(1, 2, figsize=[8, 150])

    # plotting the transcript data
    df = pd.DataFrame(myMat, index = label_list)
    cmap = sns.diverging_palette(240, 10, s=80, l=30, as_cmap=True)
    sns.heatmap(df, center = 0, cmap=cmap, vmin=-2, vmax=2, ax=axs[0])
    #sns.heatmap(df, center = 0, cmap=cmap, vmin=-2, vmax=2, ax=axs[0], cbar_kws=dict(use_gridspec=False, location="top"))

    # todo: plot the rest of the features
    # get the features:
    extra_features_mat
    # get the index
    myMat = extra_features_mat[indices,:]
    df = pd.DataFrame(myMat, index = label_list, columns = segway_track_ordered[6:])
    sns.heatmap(df, ax=axs[1])
    #sns.heatmap(df, ax=axs[1], cbar_kws=dict(use_gridspec=False, location="top"))

    plt.tight_layout()

    figFile = plotFolder + clustering_name + '.pdf'
    print(figFile)
    plt.savefig(figFile)

plt.savefig(figFile, bbox_inches='tight')

plt.close('all')

plt.show()

########################################
# 8. Examinatation of individual samples/labels
########################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> working line

# a function for, given the index, give me the length of the segment

index = 98
index = 79
index = 56
index = 77

ann = ann_info_list[index]
annAccession = ann['accession']
print(annAccession)

total_bp = 0
labels_list = list(ann['segway_anns']['clusters'].keys())
for lab in labels_list:
    total_bp += ann['segway_anns']['clusters'][lab].bp_count

# test code, for getting the reorder 
label_digit = str(2)
for i in range(16):
    label_digit = str(11)
    for lab in labels_list:
        if label_digit in lab:
            print(labels_list.index(lab))
            print(lab)


p_inds = [12, 13]
e_inds = [10, 8]
we_inds = [9] 
t_inds = [11, 4, 3]
q_inds = [2, 7, 0, 1 , 5]
b_inds = []
f_inds = [6]
c_inds = [14]
n_inds = [15]
all_inds = {'p': [10], 
            'e': [12],
            'we': [],
            't' : [2,9,3],
            'c' : [6],
            'f' : [5],
            'b' : [14],
            'q' : [4,7,0,11,8,1],
            'n' : [8]}

new_labels = ['promoter', 'enhancer', 'weak-enhancer-signal', 'transcribed', 'cons-het', 'fac-het', 'bivalent', 'quiescent', 'not identified']

relabel = {}
relabel['labels_list'] = labels_list
relabel['labels_inds'] = all_inds

# save the file
    
relabel_file = dataFolder + dataSubFolder + annAccession + '/relabel.pkl'
print(relabel_file)
with open(relabel_file, 'wb') as f:
    pickle.dump(relabel, f)

# do the bar plot - list what went in it
relabel_bp_ratios = np.zeros([len(all_inds.keys())])
relabel_segwayLabels_list = []
for i,item in enumerate(list(all_inds.keys())):
    relabel_inds = all_inds[item]
    relabel_bp_ratio = 0
    segwayLabels = ''
    for ind in relabel_inds:
        segway_label = labels_list[ind]
        segwayLabels = segwayLabels + segway_label + ' | '
        relabel_bp_ratio += ann['segway_anns']['clusters'][segway_label].bp_count / total_bp
    relabel_bp_ratios[i] = relabel_bp_ratio
    segwayLabels = segwayLabels[:-2] + ' ::: ' + new_labels[i]
    relabel_segwayLabels_list.append(segwayLabels)

fig, ax = plt.subplots(1,1, figsize=[7,3])
#ax.barh(range(0,len(all_inds)), relabel_bp_ratios, .8)

ax.tick_params(axis='both', which='major', labelsize=8)
plt.barh(relabel_segwayLabels_list, relabel_bp_ratios, .7, color = 'black')
plt.title('genome coverage', fontsize=10)
plt.tight_layout()
plotFolder_add = plotFolder + annAccession + '/'
figFile = plotFolder_add + 'genome_ratio_reorder_dist.pdf'
plt.savefig(figFile)

plt.show()
plt.close('all')

print(ann['accession'])
print(label)
print(round(ann['segway_anns']['clusters'][label].bp_count/sib, 3))

#########################################
# 9. Notes and drafts
#########################################

# exploring plots
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Questions:
# the train sample question - how many training actual value do we have? 
# normalization for the train and test data: what do we use for feature and signal distribution
# Why do we have zero Quiescent values
# summarize the transcript data?

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> draft line

        clustering_name = '%s_%d' %(method, fc)
        dendro_indices[clustering_name] = sib.dendrogram_row.reordered_ind

        figFile = '%sclustering_%s.pdf' %(plotFolder, clustering_name)
        print(figFile)
        plt.savefig(figFile)
        plt.close('all')


# doing the clustering
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> example clustering
import seaborn as sns

feature_list = segway_track_ordered[0:6]
for agg_feature in featureAgg_names:
    feature_list.append(agg_feature)
    
fig, axs = plt.subplots(1, 3, figsize =(12, 40), gridspec_kw={'width_ratios':[1,1,1]})
pcc = wholeMat_transcribed.shape[0] #150 # plot cluster count
myMat = wholeMat_transcribed[0:pcc, 0:6]
#myccc = cluster_classCode[0:pcc]
colors = whole_cluster_class_color[0:pcc,] / 255
df = pd.DataFrame(myMat, index = transcribed_cluster_list[0:pcc], columns = feature_list[0:6])
df = pd.DataFrame(myMat, index = allClusters_list[0:pcc], columns = segway_track_ordered[0:6])
df_trans = pd.DataFrame(myMat.transpose(), columns = transcribed_cluster_list[0:pcc], index = feature_list)
df_trans = pd.DataFrame(myMat.transpose(), columns = allClusters_list[0:pcc], index = segway_track_ordered[0:6])

#sns.heatmap(myMat)
#plt.show()

# the black lines are not zero, they are just very close to zero
sib = sns.clustermap(df, figsize=(6, 150), row_colors= colors, method = 'average', dendrogram_ratio = [.25,.01], cbar_pos=(.02, .985,.03,.01))
sns.set(font_scale = .3)
sib = sns.clustermap(df, row_colors= colors, method = 'ward', ax = axs[0])
sib = sns.clustermap(df_trans, figsize=(40,5), col_colors= colors, method = 'ward' )
sib = sns.clustermap(df_trans, figsize=(40,5), col_colors= colors, method = 'average')
plt.show()
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/testBatch105/'
figFile = plotFolder + 'clustering_histone_average.pdf'
print(figFile)
plt.savefig(figFile)
plt.close('all')


    #summary_zero_mat = np.zeros((thisMat.shape[0], 16))
    #for i in range(len(sortedClusterList)):
    #    for j in range(16):
    #        summary_zero_mat[i,j] = np.mean(thisMat[i, j*10:((j*10)+10)])

    plotMat = allTranscribe_mat[0:50,:]
    plotMat = summary_zero_mat
    h1 = pd.DataFrame(plotMat, index = transcribed_cluster_list[0:50])

    h1 = pd.DataFrame(np.log10(plotMat, out=np.zeros_like(plotMat), where=(plotMat!=0)), index = transcribed_cluster_list[0:50])
    cmap = sns.diverging_palette(240, 10, s=80, l=30, as_cmap=True)
    # g1 = sns.heatmap(h1, center=0, cmap=cmap)
    g1 = sns.heatmap(h1, center=0, cmap=cmap, vmin=-2, vmax=2)
    g1 = sns.heatmap(h1, center=0, cmap=cmap, vmin=-2, vmax=2, ax=axs[0])
    plt.show()

    transcribe_cluster_list.append(index_cluster)


    info = {} # keeps the info for the annotations - it also has an index which is just 0-104 for the 105 runs
    info['index'] = c
    c = c+1
    # load the overlap dataframe
    sampleFolderAdd = dataFolder + dataSubFolder + sampleFolder + '/'
    mapping_file = sampleFolderAdd + 'trackname_assay.txt'
    signal_file = sampleFolderAdd + 'signal_distribution.tab.txt'

    originalBedFile = annMeta[sampleFolder]['bedFile']
    originalBed_accession = originalBedFile.split('.')[0]

    segway_ann_sum_file = dataFolder + dataSubFolder + sampleFolder + '/' + originalBed_accession + '_annotationSummary.pkl'

    # get the segway annotation summary file
    with open(segway_ann_sum_file, 'rb') as pickledFile:
        segway_anns = pickle.load(pickledFile)

    info['segway_anns'] = segway_anns
    
    segway_cluster_list = list(segway_anns['clusters'].keys())
    # sorting the cluster labels for segway
    sortedClusterList = []
    for label in segwayLabels:
        #print(label)
        for item in segway_cluster_list:
            #print(item)
            if item.split('_')[1] == label:
                sortedClusterList.append(item)

    segway_cluster_list = sortedClusterList

    # read the mapping file
    track_assay_map = {}
    inputTrack_list = []
    with open(mapping_file) as inputFile:
        for line in inputFile:
            fields = line.strip().split()
            track_assay_map[fields[0]] = fields[1]
            inputTrack_list.append(fields[1])

    info['track_assay_map'] = track_assay_map
    info['inputTrack_list'] = inputTrack_list

    # getting the segway input track count
    segway_track_count = 0
    with open(signal_file, 'r') as inputFile:
        header = inputFile.readline().split()[0]
        previousCluster_int = inputFile.readline().split()[0]
        cluster_int = previousCluster_int
        while previousCluster_int == cluster_int:
            previousCluster_int = cluster_int
            cluster_int = inputFile.readline().split()[0]
            segway_track_count += 1

    info['segway_track_count'] = segway_track_count

    ann_info[sampleFolder] = info




    ann_clusters = ann['segway_anns']['clusters']
    for label in segwayLabels:
        #print(label)
        for cluster in ann_clusters:
            #print(item)
            if ann_clusters[cluster].biolabel == label:
                index_cluster = '%d___%s' %(index, cluster)
                allClusters_list.append(index_cluster)

    segway_cluster_list = sortedClusterList









