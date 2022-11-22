

import util # this is for features_from_segtools_dir
import gzip
import pickle
import pandas as pd

# load the model
runFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/runs/'
model_file = runFolder + "run02/model_278exp_reg0.028_auc0.86_wholeData.pickle.gz"
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
    mnemonics_file = sampleFolderAdd + 'mnemonics_v02.txt'
    data = np.column_stack([range(len(labels)), labels])
    np.savetxt(mnemonics_file, data, fmt=['%s\t', '%s'])
    # write the mnemonics file
    
    probs = np.round(the_model.predict_proba(dftr), 5)
    probs_file = sampleFolderAdd + 'probs_v02.csv'
    probsdf = pd.DataFrame(probs, index=range(probs.shape[0]), columns=the_model.classes_)
    probsdf.to_csv(probs_file)
    # write the mnemonics file 


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

