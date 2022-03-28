# Scaling up the transcription data comparison

########################################
# 0. Initials 
########################################

import linecache
import pickle
import re
import numpy as np
# TODO: not sure if it works
from QC_transcriptionComparison_util import Gene, Exon, Annotation, AnnotationClass

# General data folder
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'

dataSubFolder = 'testBatch105/fromAPI/'

# GTF data structure
fileName = dataFolder + '/geneLists.pkl'
with open(fileName, 'rb') as pickleFile:
    geneListsAndIDs = pickle.load(pickleFile)
    
geneList = geneListsAndIDs[0]
geneIDList = geneListsAndIDs[1]
del geneListsAndIDs

segwayLabels = ['Quiescent', 'ConstitutiveHet', 'FacultativeHet', 'Transcribed', 'Promoter', 'Enhancer', 'RegPermissive', 'Bivalent', 'LowConfidence']

########################################
# 1. Specifics
########################################

# list of annotation folders
sampleList = os.listdir(dataFolder + dataSubFolder)

########################################
# 2. Pre processings 
########################################

# 2.0 annotation folder prep

# NOTE: we have a class for annotations, we are not using it at the moment

annMeta = {}
c = 1
for annFolder in sampleList:
    print(c)
    c += 1

    sampleFolderAdd = dataFolder + dataSubFolder + annFolder + '/'
    bedFile = 'none'
    expressionFile = 'none'

    files = os.listdir(dataFolder + dataSubFolder + annFolder)

    for file in files:
        if file.endswith('.bed'):
            bedFile = file
        
        if file.endswith('bed.gz'):
            fileAdd = sampleFolderAdd + file
            os.system('gunzip %s' %(fileAdd))
            bedFile = file[0:-3]
        
        if file.startswith('preferred_default'):
            expressionFile = file
            
    annMeta[annFolder] = {"bedFile": bedFile, "expressionFile": expressionFile}

# 2.1. preping expression data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

annAccessionList = list(annMeta.keys())
# annAccession = annAccessionList[3]

for annAccession in annAccessionList:

    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
    print(sampleFolderAdd)

    expFileName = annMeta[annAccession]['expressionFile']
    if (expFileName != 'none'):

#        mistakeFile = sampleFolderAdd + 'geneExp_dict_ENCFF776YSH.pkl'
#        print(mistakeFile)
#        os.remove(mistakeFile)
#        print('doneRemove')

        print('doing the expression')
        expFile = sampleFolderAdd + expFileName
        print(expFile)
        expression = {}
        with open(expFile, 'r') as file:
            line = file.readline() # there is one header line

            for line in file:
                fields = line.strip().split()
                geneID = fields[0]
                transcriptID = fields[1]

                if geneID in expression:
                    expression[geneID] += np.log10(float(fields[5]) + 1)
                else:
                    expression[geneID] = np.log10(float(fields[5]) + 1)

        # saving expression dict            
        expAccession = re.split('_|\.', expFileName)[2]
        outputFile = sampleFolderAdd + 'geneExp_dict_' + expAccession + '.pkl'
        print('printing this file %s' %(outputFile))
        
        with open(outputFile, 'wb') as f:
            pickle.dump(expression, f)


# loading expression dict
inputFile = sampleFolderAdd + 'geneExp_dict_' + expAccession + '.pkl'
with open(inputFile, 'rb') as pickledFile:
    expression = pickle.load(pickledFile)


# 2.2. preping annotation data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# the function gets annotation bed file and documents the report on it for length and etc. for clusters, we also get the length distribution which is useful for classifier performance correction

# fetch the length distribution data from length_distribution

## >>>>>>>>>>>>>>>>>>>>>>>> code for one file
annAccession = annAccessionList[0]
sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
bedFileName = annMeta[annAccession]['bedFile']

bedFileAdd = sampleFolderAdd + bedFileName

annotation_generalInfo_clusters(bedFileAdd)

inputFile = sampleFolderAdd + 'ENCFF172LXH_annotationSummary.pkl'
with open(inputFile, 'rb') as pickledFile:
    annSum = pickle.load(pickledFile)

# the function
def annotation_generalInfo_clusters(bedFileAdd):

    splitFileName = bedFileAdd.split('/')
    bedFileName = splitFileName[-1].split('.')[0]
    index = bedFileAdd.index(bedFileName)
    outputFolder = bedFileAdd[0:index]

    previous_class = ''

    clusters = {}
    classes = {}
    c = 0
    with open(bedFileAdd, 'r') as annotations:

        # annotations have no header
        # header = annotations.readline()

#        f = open(bedFileAdd, 'r')
#        line = f.readline()
        
        for line in annotations:
            c += 1
            fields = line.strip().split()

            # doing the clusters first
            if fields[3] in clusters.keys():

                clusters[fields[3]].bp_count += int(fields[2]) - int(fields[1])
                clusters[fields[3]].region_count += 1
                
            else:
                
                called = fields[3]
                biolabel = fields[3].split('_')[1]
                cluster = fields[3].split('_')[0]
                color = fields[8]
                bp_count = int(fields[2]) - int(fields[1])
                region_count = 1
                region_dist = 1
                clusters[fields[3]] = Annotation(called, biolabel, cluster, color, bp_count, region_count, region_dist)


            # doing the class
            if previous_class == fields[3].split('_')[1]:
                classes[previous_class].bp_count +=  int(fields[2]) - int(fields[1])
            else:
                current_class = fields[3].split('_')[1]
                if current_class in classes.keys():
                    classes[current_class].bp_count +=  int(fields[2]) - int(fields[1])
                    classes[current_class].region_count += 1
                else:
                    clusterList = [] # not filling it now, it can be filled later using annotations 
                    biolabel = current_class
                    color = fields[8]
                    bp_count = int(fields[2]) - int(fields[1])
                    region_count = 1
                    region_dist = 1
                    classes[biolabel] = AnnotationClass(biolabel, clusterList, color, bp_count, region_count, region_dist)

            previous_class = fields[3].split('_')[1]

            
    # filling up the cluster distribution from distribution file
    
    segmentSizesFile = outputFolder + 'segment_sizes.tab.txt'
    with open(segmentSizesFile, 'r') as inputFile:
        lines = inputFile.readlines()
    
    clusterList = list(clusters.keys())
    for cluster in clusterList:
        cluster_number = int(cluster.split('_')[0])
        fields = lines[cluster_number + 2].strip().split('\t')
        info = {'num.segs': int(fields[1]), 'mean.len': float(fields[2]), 'median.len': float(fields[3]), 'stdev.len': float(fields[4]), 'num.bp': float(fields[5]), 'frac.bp': float(fields[6])}
        clusters[cluster].region_dist = info

    annotationSummary = {"classes": classes, "clusters": clusters}
    outputFile = outputFolder + bedFileName + '_annotationSummary.pkl'
    with open(outputFile, 'wb') as f:
        pickle.dump(annotationSummary, f)

    print('annotation summary saved in %s' %(outputFile))

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# code for all files:

for annAccession in annAccessionList:

    print(annAccession)
    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
    bedFileName = annMeta[annAccession]['bedFile']

    bedFileAdd = sampleFolderAdd + bedFileName

    annotation_generalInfo_clusters(bedFileAdd)


