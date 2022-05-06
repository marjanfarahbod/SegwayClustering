# here to compare the annotations for two files.

# TODO: write the annotation summary function for chmm - done
# TODO: write the annotation overlap - done
# TODO: in QC_transcriptionComparison_util.py, fix the annotation_generalInfo_clusters() for calling, if needed, the main one is currently in the QC_transcriptionComparison_02.py. 
# TODO: do the same for ccre

########################################
# 0. Initials 
########################################

import linecache
import pickle
import re
import numpy as np
import pandas as pd
import subprocess as sp
import seaborn as sns
import matplotlib.pyplot as plt
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
'''
Active TSS
Flanking TSS
Flanking TSS upstream
Flanking TSS downstream
Strong transcription
Weak transcription
Genic enhancer 1
Genic enhancer 2
Active enhancer 1
Active enhancer 2
Weak enhancer
ZNF genes & repeats
Heterochromatin
Bivalent/Poised TSS
Bivalent enhancer
Repressed PolyComb
Weak Repressed PolyComb
Quiescent/Low
'''

# TODO: fix it

# TODO load annMeta from the folder
book =  'EnhA1 EnhA2 EnhBiv EnhG1 EnhG2 EnhWk Het Quies ReprPC ReprPCWk TssA TssBiv TssFlnk TssFlnkD TssFlnkU Tx TxWk ZNF/Rpts'
chromLabels = book.split()

chmmFile = dataFolder + dataSubFolder + 'ENCSR121RCG/ChromHMM_SegwayDonor_age44_ENCFF261PTV.bed'
segwayFile = dataFolder + dataSubFolder + 'ENCSR121RCG/ENCFF172LXH_filteredSorted.bed'
ccreFile = dataFolder + dataSubFolder + 'ENCSR121RCG/ccre_age1_ENCFF174POB.bed'

ccreLabels =

# count the basepairs in chromhmm

bpSum_chmm = 0 # 2,820,428,400
with open(chmmFile, 'r') as file:
    for line in file:
        fields = line.strip().split()
        bpSum += int(fields[2]) - int(fields[1])

bpSum_segway = 0 # 3,088,458,224
with open(segwayFile, 'r') as file:
    for line in file:
        fields = line.strip().split()
        bpSum_segway += int(fields[2]) - int(fields[1])

bpSum_ccre = 0  # 253,321,371
with open(ccreFile, 'r') as file:
    for line in file:
        fields = line.strip().split()
        bpSum_ccre += int(fields[2]) - int(fields[1])

# we go by linecache

chroms = 'chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8  chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY'

chrom_list = chroms.split()

# linecache.clearcache() # TODO: use this at some point

# we by 100bp, with a chance of skipping for both chromhmm and segway - both files should be sorted TODO: sort the CCRE

# the output matrix will be overlap

# we walk on Segway 100 to 100, mark the labels.
sampleFolderAdd = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/testBatch105/fromAPI/ENCSR121RCG/'

sortedCmmFile = sampleFolderAdd + 'sorted_ChromHMM_SegwayDonor_age44_ENCFF261PTV.bed'
os.system('sort -V -k1,1 -k2,2 %s > %s' %(chmmFile, sortedCmmFile))

# use annotation_generalInfo_classes_chmm(bedFileAdd) from QC_transcriptionComparison_util.py

chmm_ann_summ_file = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/testBatch105/fromAPI/ENCSR121RCG/sorted_ChromHMM_SegwayDonor_age44_ENCFF261PTV_annotationSummary.pkl'

# the function for pickle the annotation summary for chmm is in QC_transcriptionComparison_util.py

segway_ann_summ_file = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/testBatch105/fromAPI/ENCSR121RCG/ENCFF172LXH_annotationSummary.pkl'

# load both class lists

with open(chmm_ann_summ_file, 'rb') as pickledFile:
    chmm_anns = pickle.load(pickledFile)

with open(segway_ann_summ_file, 'rb') as pickledFile:
    segway_anns = pickle.load(pickledFile)

segway_cluster_list = list(segway_anns['clusters'].keys())

# sorting the cluster labels for segway
print(segwayLabels)
sortedClusterList = []
for label in segwayLabels:
    #print(label)
    for item in segway_cluster_list:
        #print(item)
        if item.split('_')[1] == label:
            sortedClusterList.append(item)

segway_cluster_list = sortedClusterList

chmm_class_list = list(chmm_anns['classes'].keys())
# todo: sort chmm, need to figure out how for now

# define the matrix for counting based on the label summary info for both of them
segway_cluster_count = len(segway_cluster_list)
chmm_class_count = len(chmm_class_list)

overlap_mat = np.zeros((segway_cluster_count, chmm_class_count))

# read the lines

segwayLineInd = 1
segLine = linecache.getline(segwayFile, segwayLineInd)
sfields = segLine.split()
schr = sfields[0]
sstart = int(sfields[1])
send = int(sfields[2])

chmmLineInd = 1
chmmLine = linecache.getline(sortedCmmFile, chmmLineInd)
cfields = chmmLine.split()
cchr = cfields[0]
cstart = int(cfields[1])
cend = int(cfields[2])

# while not at the end 
while cchr != 'chr2' and schr != 'chr2': # while not at the end

    #print(chmmLineInd)

    if cchr == schr: # if we are in the same chromosome
        current_chr = schr
        overlap = min(send, cend) - max(sstart, cstart)
        if overlap > 0: # if there is an overlap
            sindex = segway_cluster_list.index(sfields[3])
            cindex = chmm_class_list.index(cfields[3])
            overlap_mat[sindex][cindex] += overlap # we counted the overlap

            # the one with smaller end should move forward, or both if equal ends
            if send < cend:
                segwayLineInd += 1
                segLine = linecache.getline(segwayFile, segwayLineInd)
                sfields = segLine.split()
                schr = sfields[0]
                sstart = int(sfields[1])
                send = int(sfields[2])

            else:
                if cend < send:
                    chmmLineInd += 1
                    chmmLine = linecache.getline(sortedCmmFile, chmmLineInd)
                    cfields = chmmLine.split()
                    cchr = cfields[0]
                    cstart = int(cfields[1])
                    cend = int(cfields[2])

                else: # cend == send:
                    chmmLineInd += 1
                    chmmLine = linecache.getline(sortedCmmFile, chmmLineInd)
                    cfields = chmmLine.split()
                    cchr = cfields[0]
                    cstart = int(cfields[1])
                    cend = int(cfields[2])

                    segwayLineInd += 1
                    segLine = linecache.getline(segwayFile, segwayLineInd)
                    sfields = segLine.split()
                    schr = sfields[0]
                    sstart = int(sfields[1])
                    send = int(sfields[2])

        else: # else : if overlap <= zero
            while cstart > send and schr == cchr: # if at the beginning of chromosome, chmm is missing bps - I know segway won't be missing basepairs : if overlap < 0
                segwayLineInd += 1
                segLine = linecache.getline(segwayFile, segwayLineInd)
                sfields = segLine.split()
                schr = sfields[0]
                sstart = int(sfields[1])
                send = int(sfields[2])

            if overlap == 0:
                chmmLineInd += 1
                chmmLine = linecache.getline(sortedCmmFile, chmmLineInd)
                cfields = chmmLine.split()
                cchr = cfields[0]
                cstart = int(cfields[1])
                cend = int(cfields[2])
                
                segwayLineInd += 1
                segLine = linecache.getline(segwayFile, segwayLineInd)
                sfields = segLine.split()
                schr = sfields[0]
                sstart = int(sfields[1])
                send = int(sfields[2])
                

    else: # one of them has moved forward, so the next should read until it reaches it 
        while schr == current_chr:
            segwayLineInd += 1
            segLine = linecache.getline(segwayFile, segwayLineInd)
            sfields = segLine.split()
            schr = sfields[0]
            sstart = int(sfields[1])
            send = int(sfields[2])

        while schr == current_chr:
            chmmLineInd += 1
            chmmLine = linecache.getline(sortedCmmFile, chmmLineInd)
            cfields = chmmLine.split()
            cchr = cfields[0]
            cstart = int(cfields[1])
            cend = int(cfields[2])

# the heatmap plot (will have the count of basepairs basically)

overlap_mat_colNorm = overlap_mat / np.sum(overlap_mat, axis=0)[np.newaxis, :]
overlap_mat_rowNorm = overlap_mat / np.sum(overlap_mat, axis=1)[:, np.newaxis]

segway_cluster_list
chmm_class_list

book = pd.DataFrame(overlap_mat_rowNorm, index = segway_cluster_list, columns = chmm_class_list)
book = pd.DataFrame(overlap_mat_colNorm, index = segway_cluster_list, columns = chmm_class_list)
sns.heatmap(book)
plt.show()

figFile = dataFolder + dataSubFolder + annAccession + '/exp2class_heatmap.pdf'
plt.savefig(figFile)
plt.close('all')

# scaling - for all files, but only run for one chromosome

# fix the chmm file

# run the files

# run the ccre files : I want to know where the items fall

# get all the plots and investigate 
        
# read the lines until on of them reaches chr22 (for now I am skipping 22 as well), if chrs are the same, if there is an overlap, do the overlap game.
## if chrs are not the same, the one that has the same chr should move forward until chrs are the same

# do the code for the middle overlaps

# when we are in the middle overlap, the chrs are the same.









