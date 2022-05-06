# scaling up the code annotation_comparison.py for current folders
#
# ITEMS IN THE CODE:
# ########################################
# 0. Initials
# 1. Unzip and sort the chmm files for each sample
# 2. Get the summary info for each chrom file (but we dont need it now)

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

book =  'EnhA1 EnhA2 EnhBiv EnhG1 EnhG2 EnhWk Het Quies ReprPC ReprPCWk TssA TssBiv TssFlnk TssFlnkD TssFlnkU Tx TxWk ZNF/Rpts'
chromLabels = book.split()

# get the list of subfolders

# list of annotation folders

inputFile = dataFolder + dataSubFolder + 'metaInfo.pkl'
with open(inputFile, 'rb') as pickledFile:
    annMeta = pickle.load(pickledFile)

sample_folder_list = list(annMeta.keys())

### TODO: fix chromhmm

chroms = 'chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8  chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY'

chrom_list = chroms.split()


########################################
# 1. Unzip and sort the chmm files for each sample
########################################


# a dictionary of sampleFolder and sorted chmm file
chmmFile_dict = {}

for sampleFolder in sample_folder_list:
    sampleFolderAdd = dataFolder + dataSubFolder + sampleFolder + '/'

    '''
    get list of files in the folder, pick one that starts with ChromHmm_segwayDonor - 
    catch exception if it is missing

    '''
#    fileList = os.listdir(sampleFolderAdd + 'ChromHMM_SegwayDonor_*.bed')
#   import glob
    fileList = list(glob.iglob(sampleFolderAdd + 'ChromHMM_SegwayDonor_*', recursive = True))

    # if we have chromhmm, we get the first file, if we don't, we get any chmm file
    if len(fileList) > 0:
        chmmFile = fileList[0]
        
        # unzip the file
        if chmmFile.endswith('.gz'):
            os.system('gunzip %s' %(chmmFile))
            chmmFile = chmmFile[:-3]
            
        # modify the file
        fileName = fileList[0].split('/')[-1]
        sortedCmmFile = sampleFolderAdd + 'sorted_' + fileName
        os.system('sort -V -k1,1 -k2,2 %s > %s' %(chmmFile, sortedCmmFile))
        chmmFile_dict[sampleFolder] = sortedCmmFile
        
    else:
        fileList = list(glob.iglob(sampleFolderAdd + 'ChromHMM_age*', recursive = True))
        if len(fileList) > 0:
            chmmFile = fileList[0]
            
            # unzip the file
            if chmmFile.endswith('.gz'):
                os.system('gunzip %s' %(chmmFile))
                chmmFile = chmmFile[:-3]
                
            # modify the file
            fileName = fileList[0].split('/')[-1]
            sortedCmmFile = sampleFolderAdd + 'sorted_' + fileName
            os.system('sort -V -k1,1 -k2,2 %s > %s' %(chmmFile, sortedCmmFile))
            chmmFile_dict[sampleFolder] = sortedCmmFile
        else:
            chmmFile_dict[sampleFolder] = 'none'

outputFile = dataFolder + dataSubFolder + 'chmmFile_list_dict.pkl'
with open(outputFile, 'wb') as f:
    pickle.dump(chmmFile_dict, f)

#########################################    
# 2. Get the summary info for each chrom file (but we dont need it now)
#########################################

from QC_transcriptionComparison_util import *
# TODO: fix a problem: when calling the function annotation_generalInfo_classes_chmm from QC_ranscriptionComparison_util, it gives an error that name 'pickle' is not defined. I don't know how to fix this. I tried importing the pickle in the util file or even within the function, it doesn't work. for now, I just run the function first (instead of calling it from here and expecting it to read it from the file, I just compile it once and then run the following code. It works, but I need to fix the problem) - I BELIEVE it should be enough to import it in the QC_transcriptionComparison_util, and then import this. 

sampleFolder_list = list(chmmFile_dict.keys())

count = 0
for sample in sampleFolder_list:
    print(count)
    count += 1
    if chmmFile_dict[sample] != 'none':
        annFile = chmmFile_dict[sample]
        annotation_generalInfo_classes_chmm(annFile)

#########################################    
# 3. do the comparison, save the matrix and plot
#########################################

# for segway, I need to sort the clusters for each file since the cluster id (the number) is assigned randomly and clusters are identified by this number but they also have a class label. But for chromhmm that is not needed





#########################################

# TODO : first I do segway label correction and then I do the overlap thing. For the Segway label correction I also do it with the ccre 

######
# for each sample folder, fetch the chmm, fetch the segway, fix the chmm



