# 1. clustering of Enhancer data


#########################################
# 0. Initials
#########################################
from numpy import random
import pickle
import re
import numpy as np
import pandas as pd
import subprocess as sp
import seaborn as sns
import matplotlib.pyplot as plt
import linecache
import subprocess
from QC_transcriptionComparison_util import Gene, Exon, Annotation, AnnotationClass

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/'

# load the meta file 
inputFileName = 'all235Annot_meta_corrected.pkl'
#inputFileName = 'all235Annot_meta.pkl'
inputFile = dataFolder +  inputFileName
with open(inputFile, 'rb') as f:
    allMeta = pickle.load(f)

accessionList = list(allMeta.keys())
accessionList.remove(accessionList[205])
# Segway states:
segwayStates = ['Promoter', 'PromoterFlanking', 'Enhancer', 'EnhancerLow', 'Bivalent', 'CTCF', 'Transcribed', 'K9K36', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent']


#########################################
# 1. Get chr19 files
#########################################

#coverage = np.zeros((max_region_count, int(region_length/10)+40))
# create the en file
for accession in accessionList:

    if accession == 'ENCSR424JDX':
        print(accessionList.index(accession))
        continue
    
    a = accessionList.index(accession)
 
    print(a, accession)
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']
    annFile = annotation['bedFile']

    print(annFile)

    if 'Batch105' in annFile:
        continue

 
    if annFile.endswith('gz'):
        os.system('gunzip %s' %(annFile))
        annFile = annFile[0:-3]
    else:
        os.system('gunzip %s.gz' %(annFile))


    # load the mnomics, and extract the enhancer labels for mnemonics
    label_term_mapping = {}
    mnemonics_file = annotationFolder + 'mnemonics_v04.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term

    enLabels = []
    for label in label_term_mapping.keys():
        if 'Enhancer' == label_term_mapping[label]:
            enLabels.append(label)

    one_digit_label = []
    two_digit_label = []
    for label in enLabels:
        if len(label) == 2:
            two_digit_label.append(label)
        else:
            one_digit_label.append(label)

    # making the label string for the label
    label_string = ''
    if len(one_digit_label)>0:
        label_string = label_string + '['
        for label in one_digit_label:
            label_string = label_string + label

        label_string = label_string + ']'
        if len(two_digit_label)>0:
            label_string = label_string + '|'


    if len(two_digit_label)>0:
        for label in two_digit_label:
            label_string = label_string + '(' + label + ')' + '|'
        label_string = label_string[0:-1]

    
    # pattern is different for 105 batch
    # command = "grep -E 'chr19.*\\t(%s)_\' %s" %(label_string, annFile)
    command = "grep -E 'chr19.*\\t(%s)\' %s" %(label_string, annFile)
    print(command)
    
    out = annotationFolder + 'chr19_enhOnly.bed'
    f = open(out, 'w')
    subprocess.run(command, shell=True, stdout=f)

    #os.system('gzip %s' %(annFile))

    print('chr19 enhancer regions selected')

    # if it is the first file:



# chr19 divided into 300,000 , 59000000
print((59000000/300000)*10)
 
filCount = 14 # choice
bsStep = 180000 # choice
regionCount = int(59000000/bsStep) + 1
regionVLength = int(bsStep/300) # it should have no remainders

embMat = np.zeros((234, (regionCount*filCount)))

# create the 5*1000 filters
filterMat = np.zeros((regionVLength, filCount))
for i in range(filCount):
    filterMat[:,i] = np.random.choice([0,1], size=[regionVLength], p=[2./3, 1./3])

kado = np.corrcoef(filterMat, rowvar=False)
sns.heatmap(kado)
plt.show()

for accession in accessionList:

    if accession == 'ENCSR424JDX':
        print(accessionList.index(accession))
        continue
    
    a = accessionList.index(accession)
 
    print(a, accession)
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']

    enhFile = annotationFolder + 'chr19_enhOnly.bed'

    fileMat = np.zeros((regionCount, regionVLength))

    label_term_mapping = {}
    mnemonics_file = annotationFolder + 'mnemonics_v04.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term

    enLabels = []
    for label in label_term_mapping.keys():
        if 'Enhancer' == label_term_mapping[label]:
        #if 'Enhancer' in label_term_mapping[label]:
            enLabels.append(label)
    
    with open(enhFile, 'r') as efile:
        line = efile.readline()

        for line in efile:

            # load the mnomics, and extract the enhancer labels for mnemonics

            #print(line)
            fields = line.split()
            seg_start = int(fields[1])
            seg_end = int(fields[2])
            if not('Batch105' in annotationFolder) and not(fields[3] in enLabels):
                continue

            vindex = int(seg_start/bsStep) # which vector it is
            seg_start_index = int(int((seg_start%bsStep)/300)) # which starting point
            offset = 3 - int(((seg_start%bsStep)%300)/100)# remainder of the index from the starting point of the segment
            #print(vindex, seg_start_index, offset)
            seg_length_100 = int((seg_end - seg_start)/100) # how many vector units to fill

            if seg_length_100 <= offset:
                fileMat[vindex, seg_start_index] += seg_length_100
                continue

            if seg_length_100 > offset:
                fileMat[vindex, seg_start_index] += offset
                seg_start_index +=1
                seg_start = seg_start + (offset*100)

                seg_length_300 = int((seg_end - seg_start)/300) # how many vector units to fill

                end_offset = ((seg_end - seg_start)%300)/100

                seg_end_index = seg_start_index + seg_length_300
                
                if seg_end_index < regionVLength:
                    fileMat[vindex, seg_start_index:seg_end_index] = 3
                    fileMat[vindex, seg_end_index] = end_offset
                else:
                    while seg_end_index >= regionVLength:
                        end = regionVLength-1
                        fileMat[vindex, seg_start_index:regionVLength] = 3
                        #fileMat[vindex+1, 0:end] = 3
                        seg_end_index -= regionVLength
                        seg_start_index = 0
                        vindex+=1

                    fileMat[vindex, seg_start_index:seg_end_index] = 3
                    #fileMat[vindex+1, 0:seg_end_index-1000] = 3
                    fileMat[vindex, seg_end_index] = end_offset

        print(np.sum(fileMat), np.max(fileMat))
        book = np.matmul(fileMat, filterMat)
        embMat[a, :] = book.reshape(1, (regionCount*filCount))

        
df = pd.DataFrame(embMat)
book = np.corrcoef(embMat+.1)
book = np.nan_to_num(book)
bookdf = pd.DataFrame(book)

sib = sns.clustermap(df, method = 'ward', col_cluster=False)
rowInds = sib.dendrogram_row.reordered_ind

figFile = plotFolder + 'randomClustering_03.pdf'
plt.savefig(figFile)
plt.close('all')

sibbook = sns.clustermap(bookdf, method = 'ward')
plt.show()
bookrowInds = sibbook.dendrogram_row.reordered_ind

for i,sind in enumerate(bookrowInds):
    accession = accessionList[sind]
    print(i, sind,allMeta[accession]['tissueInfo'][0])


results = {}
results['embMat'] = embMat
results['corr'] = bookdf
results['embMatCInds'] = rowInds
results['corrCInds'] = bookrowInds

file = dataFolder + 'sampleEnhancerClustering.pkl'
with open(file, 'wb') as f:
    pickle.dump(results, f)

with open(file, 'rb') as f:
    results = pickle.load(f)


 
