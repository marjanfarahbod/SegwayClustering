# for comparing ccre and chmm and segway annotations
#
# ITEMS IN THE CODE:
# ########################################
# 0. Initials
# 1. Unzip and sort the ccre
# 2. Get the summary info for each ccre file
# 3. do the comparison, save the matrix 
# 4. get the plots - the two heatmaps for now

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
# 1. Unzip and sort the ccrefiles
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
        fileName = chmmFile.split('/')[-1]
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
            fileName = chmmFile.split('/')[-1]
            sortedCmmFile = sampleFolderAdd + 'sorted_' + fileName
            os.system('sort -V -k1,1 -k2,2 %s > %s' %(chmmFile, sortedCmmFile))
            chmmFile_dict[sampleFolder] = sortedCmmFile
        else:
            chmmFile_dict[sampleFolder] = 'none'

####################################################
# just correcting the chmmFile list
####################################################
inputFile = dataFolder + dataSubFolder + 'chmmFile_list_dict_obsolete.pkl'
with open(inputFile, 'rb') as f:
    chmmFile_dict_obsolete = pickle.load(f)
            
chmmFile_dict_corrected = {}
count = 0
for sampleFolder in sampleFolder_list: # for each sample


    #sampleFolder = sampleFolder_list[count]
    print(sampleFolder)

    chmmFile = chmmFile_dict_obsolete[sampleFolder]
    print(chmmFile)
    
    if chmmFile.endswith('.gz'):
        print('it ended with .gz')
        chmmFileCorrected = chmmFile[:-3]

        print(chmmFileCorrected)

        #command = 'mv %s %s' %(chmmFile, chmmFileCorrected)
        #os.system(command)

        chmmFile_dict_corrected[sampleFolder] = chmmFileCorrected

    else:
        chmmFile_dict_corrected[sampleFolder] = chmmFile

    print(count)
    count +=1
    
chmmFile_dict = chmmFile_dict_corrected
outputFile = dataFolder + dataSubFolder + 'chmmFile_list_dict.pkl'
with open(outputFile, 'wb') as f:
    pickle.dump(chmmFile_dict, f)

#########################################    
# 2. Get the summary info for each chrom file (but we dont need it now)
#########################################

inputFile = dataFolder + dataSubFolder + 'chmmFile_list_dict.pkl'
with open(inputFile, 'rb') as f:
    chmmFile_dict = pickle.load(f)

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

# Note: for segway, I need to sort the clusters for each file since the cluster id (the number) is assigned randomly and clusters are identified by this number but they also have a class label. But for chromhmm that is not needed

# load the list of chromhmm classes
book =  'EnhA1 EnhA2 EnhBiv EnhG1 EnhG2 EnhWk Het Quies ReprPC ReprPCWk TssA TssBiv TssFlnk TssFlnkD TssFlnkU Tx TxWk ZNF/Rpts'
chmm_class_list = book.split()

inputFile = dataFolder + dataSubFolder + 'metaInfo.pkl'
with open(inputFile, 'rb') as f:
    annMeta = pickle.load(f)

count = 0
for sampleFolder in sampleFolder_list: # for each sample

    print(count)
    count += 1

    '''
    1. get the segway stuff: cluster list, annotation file name

    '''
    
    # get the name of the segway bed file from annMeta
    originalBedFile = annMeta[sampleFolder]['bedFile']
    originalBed_accession = originalBedFile.split('.')[0]

    segway_ann_sum_file = dataFolder + dataSubFolder + sampleFolder + '/' + originalBed_accession + '_annotationSummary.pkl'

    # get the segway annotation summary file
    with open(segway_ann_sum_file, 'rb') as pickledFile:
        segway_anns = pickle.load(pickledFile)

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

    # get the segway annotation file
    segwayFile = dataFolder + dataSubFolder + sampleFolder + '/' + originalBed_accession + '_filteredSorted.bed'

    '''
    2. the chmm file

    '''
    chmmFile = chmmFile_dict[sampleFolder]

    if chmmFile == 'none':
        continue

    #if chmmFile.endswith('.gz'):
    #    chmmFile = chmmFile[:-3]

    '''
    3. we have both chmmFile and segwayFile, so we can move on to comparison

    '''
    
    print(segwayFile)
    print(chmmFile)

    segwayAccession = re.search('ENCFF.*_', segwayFile)[0][:-1]
    chmmAccession = re.search('ENCFF.*\.', chmmFile)[0][:-1]

    print(segwayAccession)
    print(chmmAccession)

    segway_cluster_count = len(segway_cluster_list)
    chmm_class_count = len(chmm_class_list)

    overlap_mat = np.zeros((segway_cluster_count, chmm_class_count))

    segwayLineInd = 1
    segLine = linecache.getline(segwayFile, segwayLineInd)
    sfields = segLine.split()
    schr = sfields[0]
    sstart = int(sfields[1])
    send = int(sfields[2])

    chmmLineInd = 1
    chmmLine = linecache.getline(chmmFile, chmmLineInd)
    cfields = chmmLine.split()
    cchr = cfields[0]
    cstart = int(cfields[1])
    cend = int(cfields[2])

    # while not at the end 
    while cchr != 'chr22' and schr != 'chr22': # while not at the end

        #print(chmmLineInd)
        #print(segwayLineInd)
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
                        chmmLine = linecache.getline(chmmFile, chmmLineInd)
                        cfields = chmmLine.split()
                        cchr = cfields[0]
                        cstart = int(cfields[1])
                        cend = int(cfields[2])

                    else: # cend == send:
                        chmmLineInd += 1
                        chmmLine = linecache.getline(chmmFile, chmmLineInd)
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
                    #print('here')
                    segwayLineInd += 1
                    segLine = linecache.getline(segwayFile, segwayLineInd)
                    sfields = segLine.split()
                    schr = sfields[0]
                    sstart = int(sfields[1])
                    send = int(sfields[2])

                if overlap == 0:
                    chmmLineInd += 1
                    chmmLine = linecache.getline(chmmFile, chmmLineInd)
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
                #print('there')
                segwayLineInd += 1
                segLine = linecache.getline(segwayFile, segwayLineInd)
                sfields = segLine.split()
                schr = sfields[0]
                sstart = int(sfields[1])
                send = int(sfields[2])

            while cchr == current_chr:
                chmmLineInd += 1
                chmmLine = linecache.getline(chmmFile, chmmLineInd)
                cfields = chmmLine.split()
                cchr = cfields[0]
                cstart = int(cfields[1])
                cend = int(cfields[2])

    linecache.clearcache()
        
        # save the matrix - keep the count

    overlap_df = pd.DataFrame(overlap_mat, index = segway_cluster_list, columns = chmm_class_list)

    outputFileName = 'overlap_segway_%s_chmm_%s.pkl' %(segwayAccession, chmmAccession)
    outputFile = dataFolder + dataSubFolder + sampleFolder + '/' + outputFileName
    with open(outputFile, 'wb') as f:
        pickle.dump(overlap_df, f)


#########################################
# 4. get the plots - the heatmaps for overlap matrix
#########################################

# for each file load the matrix and the heatmap -

# the plot can be made with the background percentage of overlap as normalization. Can we tell who is favoring who? I need a panel of plots for each thing, also the track value plots.

# but perhaps we need another plot for the values. Just to get the outcome. 

# TODO : first I do segway label correction and then I do the overlap thing. For the Segway label correction I also do it with the ccre 
######
# for each sample folder, fetch the chmm, fetch the segway, fix the chmm



