# TODO: Get the classifier long heatmap for all the samples

# get the feature mat for each of the samples. put them together. (this, get from the new_apply_01.py)
# keep a track of interpretaion terms, also their clustering label 
# normalize the whole matrix like the way you did before.(this, get from classifier_train.py)
# sort based on the label 
# use some parts in the paperFigure01.py to extract the new figures

    # signal file
    # feature file
    # mnemonics file
    # mapping file

from util import *

segwayStates = ['Promoter', 'PromoterFlanking', 'Enhancer', 'EnhancerLow', 'Bivalent', 'CTCF', 'Transcribed', 'K9K36', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent']

feature_names = ['(09) initial exon', '(01) H3K9me3', '(10) initial intron', '(02) H3K27me3', '(11) internal exons', '(04) H3K4me3', "(16) 3' flanking (1000-10000 bp)", '(12) internal introns', '(03) H3K36me3', '(13) terminal exon', '(06) H3K4me1', '(14) terminal intron', "(07) 5' flanking (1000-10000 bp)", '(05) H3K27ac', "(15) 3' flanking (1-1000 bp)", "(08) 5' flanking (1-1000 bp)"]

feature_names_plotOrder = ['(09) initial exon',  "(08) 5' flanking (1-1000 bp)",'(10) initial intron', '(11) internal exons','(12) internal introns', '(13) terminal exon',  '(14) terminal intron', "(15) 3' flanking (1-1000 bp)" ,"(07) 5' flanking (1000-10000 bp)","(16) 3' flanking (1000-10000 bp)",'(04) H3K4me3','(05) H3K27ac','(06) H3K4me1', '(03) H3K36me3','(02) H3K27me3','(01) H3K9me3']

feature_names_plotOrder_renamed = ['initial exon',  "5' flanking (1-1000 bp)",'initial intron', 'internal exons','internal introns', 'terminal exon',  'terminal intron', "3' flanking (1-1000 bp)" , "5' flanking (1000-10000 bp)","3' flanking (1000-10000 bp)",'H3K4me3','H3K27ac','H3K4me1', 'H3K36me3','H3K27me3','H3K9me3']

feature_names_dict = {}
for i, oft in enumerate(feature_names_plotOrder):
    feature_names_dict[oft] = feature_names_plotOrder_renamed[i]


# load the meta file 
inputFileName = 'all235Annot_meta_corrected.pkl'
#inputFileName = 'all235Annot_meta.pkl'
inputFile = dataFolder +  inputFileName
with open(inputFile, 'rb') as f:
    allMeta = pickle.load(f)

accessionList = list(allMeta.keys())

runID_file = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/testBatch105/' + 'runID_accession_map_105run.pkl'
with open(runID_file, 'rb') as pickledFile:
    runID_accession105 = pickle.load(pickledFile)

accession_runID = {}
for runID in list(runID_accession105.keys()):
    ac = runID_accession105[runID]
    accession_runID[ac] = runID

signal_files = {}
for accession in accessionList:
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']
    if accession == 'ENCSR424JDX':
        continue

    if '38batch' in annotationFolder:
        signal_files[accession] = annotationFolder + 'segOutput/call_signal_distribution/signal_distribution.tab'

    if 'May11' in annotationFolder:
        signal_files[accession] = annotationFolder + 'call-segtools/signal_distribution.tab'

    if 'Batch105' in annotationFolder:
        runID = accession_runID[accession]
        signal_files[accession] = dataFolder + 'testBatch105/all_segtools/' + runID +  '/signal_distribution/signal_distribution.tab'

        
feature_files = {}
for accession in accessionList:
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']
    if accession == 'ENCSR424JDX':
        continue

    if '38batch' in annotationFolder:
        feature_files[accession] = annotationFolder + 'segOutput/call_feature_aggregation/feature_aggregation.tab'

    if 'May11' in annotationFolder:
        feature_files[accession] = annotationFolder + 'call-segtools/feature_aggregation.tab'

    if 'Batch105' in annotationFolder:
        runID = accession_runID[accession]
        feature_files[accession] = dataFolder + 'testBatch105/all_segtools/' + runID +  '/feature_aggregation/feature_aggregation.tab'

accessionList = list(allMeta.keys())
wholeMat = np.zeros((16*234, 16))
sindex = 0
accessionIndex = 0

termIndexList = []
labelList = []
accessionListFill = []
for accession in accessionList:
    annotation = allMeta[accession]

    if accession == 'ENCSR424JDX':
        print(accessionList.index(accession))
        continue

    annotationFolder = annotation['folder']
    print(annotationFolder)
    
    mnemFile = annotationFolder + 'mnemonics_v04.txt'
        
    label_term_mapping = {}
    with open(mnemFile, 'r') as mnemonics:
        for line in mnemonics:
            print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term

    labelCount = len(label_term_mapping)
    signalFile = signal_files[accession]
    featureFile = feature_files[accession]

    mapping_file = annotationFolder + 'trackname_assay.txt'
    track_assay_map = {}
    inputTrack_list = []
    with open(mapping_file) as inputFile:
        for line in inputFile:
            fields = line.strip().split()
            track_assay_map[fields[0]] = fields[1]
            inputTrack_list.append(fields[1])

    ann_features, ann_label_bases, ann_feature_names = features_from_segtools_dir(featureFile, signalFile, track_assay_map)

    df = pd.DataFrame(ann_features)
    dft = df.T
    dftr = dft[feature_names]

    # the index might not be the same as the label
    rowCount = dftr.shape[0]
    inds = np.linspace(0,rowCount-1, rowCount).astype(int)
    indList = [str(x) for x in inds]
    dftr = dftr.reindex(indList)

    dftr = dftr[feature_names_plotOrder]
    dftr = dftr.rename(columns=feature_names_dict)
    
    # put it inot a matrix
    wholeMat[sindex:sindex + labelCount, :] = dftr.to_numpy()
    sindex = sindex + labelCount

    for l in range(labelCount):
        term = label_term_mapping[str(l)]
        termIndex = segwayStates.index(term)
        termIndexList.append(termIndex)
        labelList.append(l)
        accessionListFill.append(accession)

    accessionIndex +=1
    print(accessionIndex)
    # put it into the final matrix list

wholeMat = wholeMat[0:sindex,:]
allFeatureMatData = {}
allFeatureMatData['mat'] = wholeMat
allFeatureMatData['colLabel'] = feature_names_plotOrder_renamed
allFeatureMatData['labelList'] = labelList
allFeatureMatData['accessionList'] = accessionListFill
allFeatureMatData['termIndex'] = termIndexList
allFeatureMatData['segwayStates'] = segwayStates

outputFile = dataFolder + 'allsamples_featureMat.pkl'
with open(outputFile, 'wb') as f:
    pickle.dump(allFeatureMatData, f)

inputFile = dataFolder + 'allsamples_featureMat.pkl'
with open(inputFile, 'rb') as f:
    allFeatureMatData = pickle.load(f)

################
# 2. get the same thing, sorted for each sample
################

wholeMat = np.zeros((16*15, 16))
sindex = 0
accessionIndex = 0

termIndexList = []
labelList = []
accessionListFill = []
for accession in accessionList[0:15]:
    annotation = allMeta[accession]

    if accession == 'ENCSR424JDX':
        print(accessionList.index(accession))
        continue

    annotationFolder = annotation['folder']
    print(annotationFolder)
    
    mnemFile = annotationFolder + 'mnemonics_v04.txt'
        
    label_term_mapping = {}
    with open(mnemFile, 'r') as mnemonics:
        for line in mnemonics:
            print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term

    labelCount = len(label_term_mapping)
    signalFile = signal_files[accession]
    featureFile = feature_files[accession]

    mapping_file = annotationFolder + 'trackname_assay.txt'
    track_assay_map = {}
    inputTrack_list = []
    with open(mapping_file) as inputFile:
        for line in inputFile:
            fields = line.strip().split()
            track_assay_map[fields[0]] = fields[1]
            inputTrack_list.append(fields[1])

    ann_features, ann_label_bases, ann_feature_names = features_from_segtools_dir(featureFile, signalFile, track_assay_map)

    df = pd.DataFrame(ann_features)
    dft = df.T
    dftr = dft[feature_names]

    # the index might not be the same as the label
    rowCount = dftr.shape[0]
    inds = np.linspace(0,rowCount-1, rowCount).astype(int)
    indList = [str(x) for x in inds]
    dftr = dftr.reindex(indList)

    dftr = dftr[feature_names_plotOrder]
    dftr = dftr.rename(columns=feature_names_dict)

    termList = []
    for l in range(labelCount):
        term = label_term_mapping[str(l)]
        termIndex = segwayStates.index(term)
        termList.append(termIndex)

        
    myMat = dftr.to_numpy()
    termList_arr = np.array(termList)
    sortedTermInds_inds = np.argsort(termList_arr)
    
    sortedMat = np.zeros(myMat.shape)
    for i in range(labelCount):
        sortedMat[i,:] = myMat[sortedTermInds_inds[i],:]
        termIndexList.append(termList_arr[sortedTermInds_inds[i]])
    
    # put it inot a matrix
    wholeMat[sindex:sindex + labelCount, :] = sortedMat
    sindex = sindex + labelCount

    print(accessionIndex)
    # put it into the final matrix list

wholeMat = wholeMat[0:sindex,:]
allFeatureMatData = {}
allFeatureMatData['mat'] = wholeMat
allFeatureMatData['colLabel'] = feature_names_plotOrder_renamed
allFeatureMatData['accessionList'] = '15 first samples from the accessionList'
allFeatureMatData['termIndex'] = termIndexList
allFeatureMatData['segwayStates'] = segwayStates

outputFile = dataFolder + 'plot15samples_termSortedFeatureMat.pkl'
with open(outputFile, 'wb') as f:
    pickle.dump(allFeatureMatData, f)


inputFile = dataFolder + 'plot15samples_termSortedFeatureMat.pkl'
with open(inputFile, 'rb') as f:
    allFeaturesMetaData = pickle.load(f)


