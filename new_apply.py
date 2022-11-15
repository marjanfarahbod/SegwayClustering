
# for each of your samples
import util

# the original 16 features that go to the classifier - just keeping it here
feature_names = ['(09) initial exon', '(01) H3K9me3', '(10) initial intron', '(02) H3K27me3', '(11) internal exons', '(04) H3K4me3', "(16) 3' flanking (1000-10000 bp)", '(12) internal introns', '(03) H3K36me3', '(13) terminal exon', '(06) H3K4me1', '(14) terminal intron', "(07) 5' flanking (1000-10000 bp)", '(05) H3K27ac', "(15) 3' flanking (1-1000 bp)", "(08) 5' flanking (1-1000 bp)"]

signal_file = '/Users/marjanfarahbod/Downloads/signal_distribution.tab'
feature_file = '/Users/marjanfarahbod/Downloads/feature_aggregation.tab'
mapping_file = '/Users/marjanfarahbod/Downloads/trackname_assay_meh.txt'
track_assay_map = {}
inputTrack_list = []
with open(mapping_file) as inputFile:
    for line in inputFile:
        fields = line.strip().split()
        track_assay_map[fields[0]] = fields[1]
        inputTrack_list.append(fields[1])

ann_features, ann_label_bases, ann_feature_names = util.features_from_segtools_dir(feature_file, signal_file, track_assay_map)

import pandas as pd
df = pd.DataFrame(ann_features)
dft = df.T
dftr = dft[feature_names]

# load the model
import gzip
import pickle
runFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/runs/'
model_file = runFolder + "run02/model_278exp_reg0.028_auc0.86_wholeData.pickle.gz"
with gzip.open(model_file, "r") as f:
    the_model = pickle.load(f)

labels = the_model.predict(dftr)



