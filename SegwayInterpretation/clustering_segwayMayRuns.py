# To load the data, make the clustering matrix, do the clustering.

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
dataSubFolder = 'testBatch_May112022/'
dataSubFolder = 'the38batch/'

##### IMPORTANT: classifier training data
feature_names = ['(09) initial exon', '(01) H3K9me3', '(10) initial intron', '(02) H3K27me3', '(11) internal exons', '(04) H3K4me3', "(16) 3' flanking (1000-10000 bp)", '(12) internal introns', '(03) H3K36me3', '(13) terminal exon', '(06) H3K4me1', '(14) terminal intron', "(07) 5' flanking (1000-10000 bp)", '(05) H3K27ac', "(15) 3' flanking (1-1000 bp)", "(08) 5' flanking (1-1000 bp)"]

# this is for plot only, just renamed
feature_names_mod = ['initial exon', 'H3K9me3', 'initial intron', 'H3K27me3', 'internal exons', 'H3K4me3', "3' flanking (1000-10000 bp)", 'internal introns', 'H3K36me3', 'terminal exon', 'H3K4me1', 'terminal intron', "5' flanking (1000-10000 bp)", 'H3K27ac', "3' flanking (1-1000 bp)", "5' flanking (1-1000 bp)"]

# This is for plotting only - reordered
feature_list_reorder = ['initial exon', "5' flanking (1-1000 bp)",'initial intron', 'internal exons', 'internal introns',  'terminal exon', 'terminal intron', "3' flanking (1-1000 bp)", "5' flanking (1000-10000 bp)", "3' flanking (1000-10000 bp)", 'H3K4me3', 'H3K27ac', 'H3K4me1', 'H3K36me3', 'H3K27me3',  'H3K9me3']


# get the list of folders / May92
runIDFile = 'runID_list.pkl'
inputFile = dataFolder + dataSubFolder + runIDFile
with open(inputFile, 'rb') as f:
    runIDs = pickle.load(f)

# get the list of folders / 38
inputFile = dataFolder + dataSubFolder + 'hg_accessionList.pkl'
with open(inputFile, 'rb') as f:
    accessionList = pickle.load(f)


classifier_tracks = ['H3K36me3', 'H3K27me3', 'H3K4me3', 'H3K9me3', 'H3K27ac', 'H3K4me1']
segway_track_ordered = classifier_tracks
segway_track_ordered.append('CTCF')
segway_track_ordered.append('DNase-seq')
segway_track_ordered.append('ATAC-seq')
segway_track_ordered.append('POLR2A')
segway_track_ordered.append('EP300')

runID = '5857d68b-e559-4776-9c12-a6e10aea7f76' # this one has a problem in the May batch for some reason.
removeInd = runIDs.index(runID)
runIDs.remove(runID)

allFeatureAgg_mat = np.zeros((92*16, 16)) # this one fills with index
agg_cluster_list = []
afa_index = 0
for i,runID in enumerate(runIDs):

    index = i
    print(index)

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

df = pd.DataFrame(allFeatureAgg_mat, index = agg_cluster_list, columns = feature_names_mod)

sib = sns.clustermap(df, figsize=(6, 150), method = 'ward', dendrogram_ratio = [.25,.01], cbar_pos=(.02, .985,.03,.01))
sns.set(font_scale = .4)
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/testBatch92/'
figFile = '%sclustergram_ward_May92runs_V04.pdf' %(plotFolder)
print(figFile)
plt.savefig(figFile)
plt.close('all')

dendro_indices = sib.dendrogram_row.reordered_ind
dendro_feature_indices = sib.dendrogram_col.reordered_ind

clustering_data = {}
clustering_data['indices'] = dendro_indices
clustering_data['labels'] = agg_cluster_list
clustering_data['sampleIDs'] = IDs
clustering_data['featureInds'] = dendro_feature_indices

outputFileName = 'clusterings_dendro_indices_92runs.pkl'
outputFile = dataFolder + dataSubFolder + outputFileName
with open(outputFile, 'wb') as f:
    pickle.dump(clustering_data, f)

runFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/runs/'
model_file = runFolder + "run03/model_300_reg.020_auc0.89V04.pickle.gz"
with gzip.open(model_file, "rb") as f:
    the_model = pickle.load(f)

allFeatureAgg_mat  # features - see section 2 of this code
all_labels = the_model.predict(allFeatureAgg_mat) # classifier output
probs = the_model.predict_proba(allFeatureAgg_mat)
agg_cluster_list # cluster id
plot_original_labels = [agg_cluster_list[i] for i in dendro_indices] # list of original labels

colormapping = {'Quiescent':[255,255,255], 'Promoter':[255,0,0], 'PromoterFlanking':[255,69, 0],'CTCF':[138,236,208], 'K9K36':[102,205,170], 'FacultativeHet':[128,0,128], 'Enhancer':[255,195,77], 'EnhancerLow':[255, 255, 0], 'Bivalent':[200,200,100], 'Transcribed':[0,128,0], 'ConstitutiveHet':[138,145,208], 'badSample1': [0,0,0], 'badSample2': [70,70,70]}

# getting the predicted labels and color for the heatmap plot
plot_predicted_labels = []
plot_colors = []
plot_labels = []
p_index_list = []
for o_label in plot_original_labels:
    tokens = o_label.split('_')
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


from scipy.stats import zscore
plot_data_z = zscore(allFeatureAgg_mat[dendro_indices,], axis = 0)
plot_data_z_thr = np.where(plot_data_z > 4, 4, plot_data_z)

df = pd.DataFrame(plot_data_z_thr, index = plot_predicted_labels, columns = feature_names_mod)
df = df.apply(zscore, axis=0)
df = df[feature_list_reorder]

# the black lines are not zero, they are just very close to zero
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/testBatch105/'
sib = sns.clustermap(df, figsize=(6, 150), row_colors= plot_colors, row_cluster=False, col_cluster=False, dendrogram_ratio = [.25,.01], cbar_pos=(.02, .985,.03,.01))
sns.set(font_scale = .4)
figFile = '%sclustergram_ward_May92runs_thr_zscore_classifierLabels_testrun_V04.pdf' %(plotFolder)
print(figFile)
plt.savefig(figFile)
plt.close('all')

sib = sns.clustermap(df, figsize=(fw, 150), row_colors= colors, method = method, dendrogram_ratio = [.25,.01], cbar_pos=(.02, .985,.03,.01))
sns.set(font_scale = .4)

# to the order of the classifier
feature_names = ['(09) initial exon', '(01) H3K9me3', '(10) initial intron', '(02) H3K27me3', '(11) internal exons', '(04) H3K4me3', "(16) 3' flanking (1000-10000 bp)", '(12) internal introns', '(03) H3K36me3', '(13) terminal exon', '(06) H3K4me1', '(14) terminal intron', "(07) 5' flanking (1000-10000 bp)", '(05) H3K27ac', "(15) 3' flanking (1-1000 bp)", "(08) 5' flanking (1-1000 bp)"]


########################################
allFeatureAgg_mat = np.zeros((38*16, 16)) # this one fills with index
agg_cluster_list = []
afa_index = 0
for i,accession in enumerate(accessionList):

    index = i
    print(index)

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

    for cluster in ann_features:
        print(cluster)
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

df = pd.DataFrame(allFeatureAgg_mat, index = agg_cluster_list, columns = feature_names_mod)

sib = sns.clustermap(df, figsize=(6, 150), method = 'ward', dendrogram_ratio = [.25,.01], cbar_pos=(.02, .985,.03,.01))
sns.set(font_scale = .4)
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/testBatch92/'
figFile = '%sclustergram_ward_May92runs_V04.pdf' %(plotFolder)
print(figFile)
plt.savefig(figFile)
plt.close('all')

dendro_indices = sib.dendrogram_row.reordered_ind
dendro_feature_indices = sib.dendrogram_col.reordered_ind

clustering_data = {}
clustering_data['indices'] = dendro_indices
clustering_data['labels'] = agg_cluster_list
clustering_data['sampleIDs'] = IDs
clustering_data['featureInds'] = dendro_feature_indices

outputFileName = 'clusterings_dendro_indices_92runs.pkl'
outputFile = dataFolder + dataSubFolder + outputFileName
with open(outputFile, 'wb') as f:
    pickle.dump(clustering_data, f)

runFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/runs/'
model_file = runFolder + "run03/model_300_reg.020_auc0.89V04.pickle.gz"
with gzip.open(model_file, "rb") as f:
    the_model = pickle.load(f)

allFeatureAgg_mat  # features - see section 2 of this code
all_labels = the_model.predict(allFeatureAgg_mat) # classifier output
probs = the_model.predict_proba(allFeatureAgg_mat)
agg_cluster_list # cluster id
plot_original_labels = [agg_cluster_list[i] for i in dendro_indices] # list of original labels

colormapping = {'Quiescent':[255,255,255], 'Promoter':[255,0,0], 'PromoterFlanking':[255,69, 0],'CTCF':[138,236,208], 'K9K36':[102,205,170], 'FacultativeHet':[128,0,128], 'Enhancer':[255,195,77], 'EnhancerLow':[255, 255, 0], 'Bivalent':[200,200,100], 'Transcribed':[0,128,0], 'ConstitutiveHet':[138,145,208], 'badSample1': [0,0,0], 'badSample2': [70,70,70]}

# getting the predicted labels and color for the heatmap plot
plot_predicted_labels = []
plot_colors = []
plot_labels = []
p_index_list = []
for o_label in plot_original_labels:
    tokens = o_label.split('_')
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


from scipy.stats import zscore
plot_data_z = zscore(allFeatureAgg_mat[dendro_indices,], axis = 0)
plot_data_z_thr = np.where(plot_data_z > 4, 4, plot_data_z)

df = pd.DataFrame(plot_data_z_thr, index = plot_predicted_labels, columns = feature_names_mod)
df = df.apply(zscore, axis=0)
df = df[feature_list_reorder]

# the black lines are not zero, they are just very close to zero
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/the38batch/'
sib = sns.clustermap(df, figsize=(6, 70), row_colors= plot_colors, row_cluster=False, col_cluster=False, dendrogram_ratio = [.25,.01], cbar_pos=(.02, .985,.03,.01))
sns.set(font_scale = .4)
figFile = '%sclustergram_ward_the38batch_thr_zscore_classifierLabels_testrun_V04.pdf' %(plotFolder)
print(figFile)
plt.savefig(figFile)
plt.close('all')

# Use the model to get the labels, add figures and all.


