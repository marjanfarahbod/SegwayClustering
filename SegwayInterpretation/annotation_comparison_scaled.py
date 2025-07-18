# scaling up the code annotation_comparison.py for current folders
#
# ITEMS IN THE CODE:
# ########################################
# 0. Initials
# 1. Unzip and sort the chmm files for each sample
# 2. Get the summary info for each chrom file (we dont need it now)
# 3. do the comparison, save the matrix and plot
# 4. get the plots - the two heatmaps for now
# 5. get the plots - with the updated classifier labels 
#
#
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
from scipy import stats
from QC_transcriptionComparison_util import Gene, Exon, Annotation, AnnotationClass

import glob

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

# segwayLabels = ['Quiescent', 'ConstitutiveHet', 'FacultativeHet', 'Transcribed', 'Promoter', 'Enhancer', 'RegPermissive', 'Bivalent', 'LowConfidence'] # obsolete
segwayLabels = ['Enhancer', 'EnhancerLow', 'Promoter', 'PromoterFlanking', 'Transcribed', 'CTCF', 'K9K36', 'Bivalent', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent']    

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

    # unzip the chmm file
    # os.system('gunzip %s' %(chmmFile + '.gz'))

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
                sindex = int(sfields[3].split('_')[0])
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
    
inputFile = dataFolder + dataSubFolder + 'biosample_tissue_info.pkl'
with open(inputFile, 'rb') as f:
    tissue_info = pickle.load( f)
    
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/testBatch105/'
for sampleFolder in sampleFolder_list:
    
    # load the overlap dataframe
    sampleFolderAdd = dataFolder + dataSubFolder + sampleFolder + '/' 
    fileList = list(glob.iglob(sampleFolderAdd + 'overlap_segway_*_chmm_*.pkl', recursive = True))

    if len(fileList) == 0:
        print('no overlap file')
        continue
    
    inputFile = fileList[0]

    with open(inputFile, 'rb') as f:
        overlap_df = pickle.load(f)

    # do the normalization:
    # sum of each column or row will be the total bp for chmm/segway. total sum of columns should equal total soum of rows.
    overlap_mat = overlap_df.to_numpy()

    # total bp count that is covered by both annotations
    total_bp = np.sum(overlap_mat)

    # fraction of base pairs in each label to the total bp count
    chmm_labelFraction = np.sum(overlap_mat, axis = 0) / total_bp
    segway_labelFraction = np.sum(overlap_mat, axis = 1) / total_bp

    # get the observed versus the expected fraction
    expectedFraction_mat = np.zeros(overlap_mat.shape)
    for i in range(len(segway_labelFraction)):
        for j in range(len(chmm_labelFraction)):
            expectedFraction_mat[i][j] = segway_labelFraction[i] * chmm_labelFraction[j]

    observedFraction_mat = np.zeros(overlap_mat.shape)
    for i in range(len(segway_labelFraction)):
        for j in range(len(chmm_labelFraction)):
            observedFraction_mat[i][j] = overlap_mat[i][j] / total_bp

    obs_exp = np.divide(observedFraction_mat, expectedFraction_mat)

    # get the normalized value
    overlap_mat_colNorm = overlap_mat / np.sum(overlap_mat, axis=0)[np.newaxis, :] # chmm
    overlap_mat_rowNorm = overlap_mat / np.sum(overlap_mat, axis=1)[:, np.newaxis] # segway

    # get the axis labels
    segway_cluster_list = list(overlap_df.index.values)
    chmm_class_list = list(overlap_df.columns.values)

    # add the ratio to the axis labels
    segwayAxis_list = []
    for i, cluster in enumerate(segway_cluster_list):
        segwayAxis_list.append(cluster + '_' + str(round(segway_labelFraction[i], 4)))

    chmmAxis_list = []
    for i, chmmclass in enumerate(chmm_class_list):
        chmmAxis_list.append(chmmclass + '_' + str(round(chmm_labelFraction[i], 4)))

    # plot the three heatmaps

    fig, axs = plt.subplots(2, 2, figsize=(12,8))

    obs_exp_log = np.log10(obs_exp, out=np.zeros_like(obs_exp), where=(obs_exp!=0))
    obs_exp_log = np.where(obs_exp_log < 0, 0, obs_exp_log)
    
    # heatmap 1: the observed over expected ratio
    h1 = pd.DataFrame(obs_exp_log, index = segwayAxis_list, columns = chmmAxis_list)
    #cmap = sns.diverging_palette(240, 10, s=100, l=30, as_cmap=True)
    cmap = sns.diverging_palette(240, 10, s=80, l=30, as_cmap=True)
    g1 = sns.heatmap(h1, center = 0,cmap=cmap, ax=axs[0,0])
    #g1.set_ylabel('this')
    #g1.set_xlabel('that')
    g1.set_title('ratio - log (observed to expected)')
    #g1.set_xticklabels(g1.get_xticklabels(), rotation=45)
    #plt.ylabel('segway')
    #plt.ylabel('chmm')
    #plt.subplots_adjust(left=.075, right=.95, top=.9, bottom=.25)
    sns.set(font_scale=1)
    plt.tight_layout()
    #plt.title('ratio - observed vs expected')
    #plt.show()

    # h2: the segway ratio covered
    h2 = pd.DataFrame(overlap_mat_rowNorm, index = segwayAxis_list, columns = chmmAxis_list)
    cmap = sns.diverging_palette(240, 10, s=80, l=30, as_cmap=True)
    g2 = sns.heatmap(h2, center = 0,cmap=cmap, ax=axs[0,1])
    #g1.set_ylabel('this')
    #g1.set_xlabel('that')
    g2.set_title('ratio of bp overlap - Segway')
    #plt.ylabel('segway')
    #plt.ylabel('chmm')
    #plt.subplots_adjust(left=.075, right=.95, top=.9, bottom=.25)
    plt.tight_layout()

    #plt.show()

    # h3: track average signal (for the tracks in the thing only)
    # get the mapping from the cluster label , get the mapping from track assay, build the matrix and do the meeting
    mapping_file = sampleFolderAdd + 'trackname_assay.txt'
    signal_file = sampleFolderAdd + 'signal_distribution.tab.txt'

    # read the mapping_file
    track_assay_map = {}
    inputTrack_list = [] # list of the track type for the sample 
    with open(mapping_file) as inputFile:
        for line in inputFile:
            fields = line.strip().split()
            track_assay_map[fields[0]] = fields[1]
            inputTrack_list.append(fields[1])

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

    
    cluster_class_map = {}
    for cluster_label in segway_cluster_list:
        int_label = cluster_label.split('_')[0]
        cluster_class_map[int_label] = cluster_label


    # filling the signal_dist_mat for the plot 
    signal_dist_mat = np.zeros((len(segway_cluster_list), segway_track_count))
    with open(signal_file, 'r') as inputFile:
        header = inputFile.readline()
        for line in inputFile:
            fields = line.strip().split()
            track = track_assay_map[fields[1]]
            track_ind = inputTrack_list.index(track) # track list order
            segway_cluster = cluster_class_map[fields[0]]
            cluster_ind = segway_cluster_list.index(segway_cluster) # cluster list order
            signal_dist_mat[cluster_ind][track_ind] = round(float(fields[2]), 4)
            
    z_signal_dist_mat = stats.zscore(signal_dist_mat, axis = 0)

    # make the dataframe
    h3 = pd.DataFrame(z_signal_dist_mat, index = segway_cluster_list, columns = inputTrack_list)
    cmap = sns.diverging_palette(240, 10, s=80, l=30, as_cmap=True)
    g3 = sns.heatmap(h3, center = 0,cmap=cmap, ax=axs[1,0])
    #g1.set_ylabel('this')
    #g1.set_xlabel('that')
    g3.set_title('mean track value - zscore')
    #plt.ylabel('segway')
    #plt.ylabel('chmm')
    #plt.subplots_adjust(left=.075, right=.95, top=.9, bottom=.25)
    plt.tight_layout()
    #plt.show()

    # h4: probs for each segway label
    #fig, axs = plt.subplots(2, 2, figsize=(12,8))

    # reorder the dataframe 
    
    probs_file = sampleFolderAdd + 'probs.txt'
    if sampleFolder == 'ENCSR614XBK':
        print('no probs')
        continue
    h4 = pd.read_table(probs_file)
    h4.drop('label', inplace=True, axis=1)
    for int_label in list(cluster_class_map.keys()):
        h4.rename(index={int(int_label) : cluster_class_map[int_label]}, inplace=True)

    h4 = h4.reindex(index = segway_cluster_list)
    h4 = h4.reindex(columns = segwayLabels[0:8])
    cmap = sns.diverging_palette(240, 10, s=80, l=30, as_cmap=True)
    g4 = sns.heatmap(h4, center = 0,cmap=cmap, ax=axs[1,1])
    #g1.set_ylabel('this')
    #g1.set_xlabel('that')
    g4.set_title('probs - classifier')
    #plt.ylabel('segway')
    #plt.ylabel('chmm')
    #plt.subplots_adjust(left=.075, right=.95, top=.9, bottom=.25)
    plt.tight_layout()

    plt.subplots_adjust( top=.9)

    fig.suptitle(sampleFolder + ' - ' +tissue_info[sampleFolder][0] + ' - ' + tissue_info[sampleFolder][1])
    #plt.show()
    plotFolder_add = plotFolder + sampleFolder + '/'
    figFile = plotFolder_add + 'the_panel_02.pdf'
    print(figFile)
    plt.savefig(figFile)
    plt.close('all')
    
    
    # H: a whole other plot for total track values? 
    
    book = pd.DataFrame(overlap_mat_rowNorm, index = segway_cluster_list, columns = chmm_class_list)
    book = pd.DataFrame(overlap_mat_colNorm, index = segway_cluster_list, columns = chmm_class_list)
    sns.heatmap(book)
    plt.ylabel('segway')
    plt.ylabel('chmm')
    plt.title('segway ratios')
    plt.show()

    # plot the probs and the heatmap of comparison

    # add the count and other info
    
#########################################
# 5. get the plots - with the updated classifier labels 
#########################################

segwayLabels = ['Enhancer_low', 'Enhancer', 'Promoter_flanking', 'Promoter', 'Transcribed', 'CTCF', 'K9K36', 'Bivalent', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent', 'Unclassified']

track_order = [ 'H3K4me3', 'H3K27ac', 'H3K4me1', 'H3K36me3', 'H3K27me3', 'H3K9me3', 'CTCF', 'DNase-seq', 'ATAC-seq', 'POLR2A', 'EP300']

inputFile = dataFolder + dataSubFolder + 'biosample_tissue_info.pkl'
with open(inputFile, 'rb') as f:
    tissue_info = pickle.load( f)

inputFileName = 'all_annInfo_list.pkl'
inputFile = dataFolder + dataSubFolder + inputFileName
with open(inputFile, 'rb') as f:
    annInfo_list = pickle.load(f)

plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/testBatch105/'
for ann in annInfo_list:
    
    index = ann['index']
    print(index)
    sampleFolder =  ann['accession']
    
    sampleFolderAdd = dataFolder + dataSubFolder + sampleFolder + '/'

    # load the updated mnemonics
    label_term_mapping = {}
    mnemonics_file = sampleFolderAdd + 'mnemonics_v02.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term
    
    # load the overlap dataframe
    sampleFolderAdd = dataFolder + dataSubFolder + sampleFolder + '/' 
    fileList = list(glob.iglob(sampleFolderAdd + 'overlap_segway_*_chmm_*.pkl', recursive = True))

    if len(fileList) == 0:
        print('no overlap file')
        continue
    
    inputFile = fileList[0]

    with open(inputFile, 'rb') as f:
        overlap_df = pickle.load(f)

    # do the normalization:
    # sum of each column or row will be the total bp for chmm/segway. total sum of columns should equal total soum of rows.
    overlap_mat = overlap_df.to_numpy()

    # total bp count that is covered by both annotations
    total_bp = np.sum(overlap_mat)

    # fraction of base pairs in each label to the total bp count
    chmm_labelFraction = np.sum(overlap_mat, axis = 0) / total_bp
    segway_labelFraction = np.sum(overlap_mat, axis = 1) / total_bp

    # get the observed versus the expected fraction
    expectedFraction_mat = np.zeros(overlap_mat.shape)
    for i in range(len(segway_labelFraction)):
        for j in range(len(chmm_labelFraction)):
            expectedFraction_mat[i][j] = segway_labelFraction[i] * chmm_labelFraction[j]

    observedFraction_mat = np.zeros(overlap_mat.shape)
    for i in range(len(segway_labelFraction)):
        for j in range(len(chmm_labelFraction)):
            observedFraction_mat[i][j] = overlap_mat[i][j] / total_bp

    obs_exp = np.divide(observedFraction_mat, expectedFraction_mat)

    # get the normalized value
    overlap_mat_colNorm = overlap_mat / np.sum(overlap_mat, axis=0)[np.newaxis, :] # chmm
    overlap_mat_rowNorm = overlap_mat / np.sum(overlap_mat, axis=1)[:, np.newaxis] # segway

    # get the axis labels
    segway_cluster_list_old = list(overlap_df.index.values)
    chmm_class_list = list(overlap_df.columns.values)

    # change segway_cluster_list_updated based on the mnemonics
    segway_cluster_list = []
    for cluster in segway_cluster_list_old:
        label = cluster.split('_')[0]
        term = label_term_mapping[label]
        segway_cluster_list.append(label + '_' + term)

    # reorder the cluster list
    segway_cluster_list_reordered = []
    for label in segwayLabels:
        for cluster in segway_cluster_list:
            if cluster.split('_')[1] == label:
                segway_cluster_list_reordered.append(cluster)

    # add the ratio to the axis labels
    segwayAxis_list = []
    for i, cluster in enumerate(segway_cluster_list):
        segwayAxis_list.append(cluster + '_' + str(round(segway_labelFraction[i], 4)))


    segwayAxis_list_reordered = []
    for label in segwayLabels:
        for axis in segwayAxis_list:
            if axis.split('_')[1] == label:
                segwayAxis_list_reordered.append(axis)

    chmmAxis_list = []
    for i, chmmclass in enumerate(chmm_class_list):
        chmmAxis_list.append(chmmclass + '_' + str(round(chmm_labelFraction[i], 4)))

    # plot the three heatmaps

    fig, axs = plt.subplots(2, 2, figsize=(12,8))

    obs_exp_log = np.log10(obs_exp, out=np.zeros_like(obs_exp), where=(obs_exp!=0))
    obs_exp_log = np.where(obs_exp_log < 0, 0, obs_exp_log)
    
    # heatmap 1: the observed over expected ratio
    h1 = pd.DataFrame(obs_exp_log, index = segwayAxis_list, columns = chmmAxis_list)
    h1 = h1.reindex(segwayAxis_list_reordered)
    #cmap = sns.diverging_palette(240, 10, s=100, l=30, as_cmap=True)
    cmap = sns.diverging_palette(240, 10, s=80, l=30, as_cmap=True)
    g1 = sns.heatmap(h1, center = 0,cmap=cmap, ax=axs[0,0])
    #g1.set_ylabel('this')
    #g1.set_xlabel('that')
    g1.set_title('ratio - log (observed to expected)')
    #g1.set_xticklabels(g1.get_xticklabels(), rotation=45)
    #plt.ylabel('segway')
    #plt.ylabel('chmm')
    #plt.subplots_adjust(left=.075, right=.95, top=.9, bottom=.25)
    sns.set(font_scale=.8)
    plt.tight_layout()
    #plt.title('ratio - observed vs expected')
    #plt.show()
    # h2: the segway ratio covered
    h2 = pd.DataFrame(overlap_mat_rowNorm, index = segwayAxis_list, columns = chmmAxis_list)
    h2 = h2.reindex(segwayAxis_list_reordered)
    cmap = sns.diverging_palette(240, 10, s=80, l=30, as_cmap=True)
    g2 = sns.heatmap(h2, center = 0,cmap=cmap, ax=axs[0,1])
    #g1.set_ylabel('this')
    #g1.set_xlabel('that')
    g2.set_title('ratio of bp overlap - Segway')
    #plt.ylabel('segway')
    #plt.ylabel('chmm')
    #plt.subplots_adjust(left=.075, right=.95, top=.9, bottom=.25)
    plt.tight_layout()

    #plt.show()

    # h3: track average signal (for the tracks in the thing only)
    # get the mapping from the cluster label , get the mapping from track assay, build the matrix and do the meeting
    mapping_file = sampleFolderAdd + 'trackname_assay.txt'
    signal_file = sampleFolderAdd + 'signal_distribution.tab.txt'

    # read the mapping_file
    track_assay_map = {}
    inputTrack_list = [] # list of the track type for the sample 
    with open(mapping_file) as inputFile:
        for line in inputFile:
            fields = line.strip().split()
            track_assay_map[fields[0]] = fields[1]
            inputTrack_list.append(fields[1])

    # reorder the input track list
    inputTrack_list_ordered = []
    for track in track_order:
        if track in inputTrack_list:
            inputTrack_list_ordered.append(track)

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

    
    cluster_class_map = {}
    for cluster_label in segway_cluster_list:
        int_label = cluster_label.split('_')[0]
        cluster_class_map[int_label] = cluster_label


    # filling the signal_dist_mat for the plot 
    signal_dist_mat = np.zeros((len(segway_cluster_list), segway_track_count))
    with open(signal_file, 'r') as inputFile:
        header = inputFile.readline()
        for line in inputFile:
            fields = line.strip().split()
            track = track_assay_map[fields[1]]
            track_ind = inputTrack_list_ordered.index(track) # track list order
            segway_cluster = cluster_class_map[fields[0]]
            cluster_ind = segway_cluster_list.index(segway_cluster) # cluster list order
            signal_dist_mat[cluster_ind][track_ind] = round(float(fields[2]), 4)
            
    z_signal_dist_mat = stats.zscore(signal_dist_mat, axis = 0)

    # make the dataframe
    h3 = pd.DataFrame(z_signal_dist_mat, index = segway_cluster_list, columns = inputTrack_list_ordered)
    h3 = h3.reindex(segway_cluster_list_reordered)
    cmap = sns.diverging_palette(240, 10, s=80, l=30, as_cmap=True)
    g3 = sns.heatmap(h3, center = 0,cmap=cmap, ax=axs[1,0])
    #g1.set_ylabel('this')
    #g1.set_xlabel('that')
    g3.set_title('mean track value - zscore')
    #plt.ylabel('segway')
    #plt.ylabel('chmm')
    #plt.subplots_adjust(left=.075, right=.95, top=.9, bottom=.25)
    plt.tight_layout()
    #plt.show()

    # h4: probs for each segway label
    #fig, axs = plt.subplots(2, 2, figsize=(12,8))

    # reorder the dataframe 
    
    probs_file = sampleFolderAdd + 'probs_v02.csv'
    #if sampleFolder == 'ENCSR614XBK':
    #    print('no probs')
    #    continue
    #h4 = pd.read_table(probs_file)
    h4 = pd.read_csv(probs_file)
    # h4.drop('label', inplace=True, axis=1)
    h4.drop(columns=h4.columns[0], inplace=True, axis=1)
    for int_label in list(cluster_class_map.keys()):
        h4.rename(index={int(int_label) : cluster_class_map[int_label]}, inplace=True)
        
    #fig, axs = plt.subplots(2, 2, figsize=(12,8))
    #h4 = h4.reindex(index = segway_cluster_list)
    h4 = h4.reindex(segway_cluster_list_reordered)
    #h4 = h4.reindex(columns = segwayLabels[0:8])
    sns.set(font_scale=1)
    cmap = sns.diverging_palette(240, 10, s=80, l=30, as_cmap=True)
    g4 = sns.heatmap(h4, center = 0,cmap=cmap, ax=axs[1,1], yticklabels=True)
    #g1.set_ylabel('this')
    #g1.set_xlabel('that')
    g4.set_title('probs - classifier')
    #plt.ylabel('segway')
    #plt.ylabel('chmm')
    #plt.subplots_adjust(left=.075, right=.95, top=.9, bottom=.25)
    plt.tight_layout()

    plt.subplots_adjust( top=.9)

    fig.suptitle(sampleFolder + ' - ' +tissue_info[sampleFolder][0] + ' - ' + tissue_info[sampleFolder][1])
    #plt.show()
    plotFolder_add = plotFolder + sampleFolder + '/'
    figFile = plotFolder_add + 'the_panel_03.pdf'
    print(figFile)
    plt.savefig(figFile)
    plt.close('all')

# the plot can be made with the background percentage of overlap as normalization. Can we tell who is favoring who? I need a panel of plots for each thing, also the track value plots.

# but perhaps we need another plot for the values. Just to get the outcome. 

# TODO : first I do segway label correction and then I do the overlap thing. For the Segway label correction I also do it with the ccre 
######
# for each sample folder, fetch the chmm, fetch the segway, fix the chmm


# Note: for segway, I need to sort the clusters for each file since the cluster id (the number) is assigned randomly and clusters are identified by this number but they also have a class label. But for chromhmm that is not needed

# 6. getting overlap for others 
########################################
# load the list of chromhmm classes
book =  'EnhA1 EnhA2 EnhBiv EnhG1 EnhG2 EnhWk Het Quies ReprPC ReprPCWk TssA TssBiv TssFlnk TssFlnkD TssFlnkU Tx TxWk ZNF/Rpts'
chmm_class_list = book.split()

inputFile = dataFolder + dataSubFolder + 'metaInfo.pkl'
with open(inputFile, 'rb') as f:
    annMeta = pickle.load(f)

annAccessionList = list(annMeta.keys())

count = 0
for annAccession in annAccessionList[1:]: # for each sample

    print(count)
    count += 1

    '''
    1. get the segway stuff: cluster list, annotation file name

     '''
    
    # get the name of the segway bed file from annMeta
    originalBedFile = annMeta[annAccession]['bedFile']
    originalBed_accession = originalBedFile.split('.')[0]

    segway_ann_sum_file = dataFolder + dataSubFolder + annAccession + '/' + originalBed_accession + '_annotationSummary.pkl'

    # get the segway annotation summary file
    with open(segway_ann_sum_file, 'rb') as pickledFile:
        segway_anns = pickle.load(pickledFile)

    segway_cluster_list = list(segway_anns['clusters'].keys())

    # get the segway annotation file
    segwayFile = dataFolder + dataSubFolder + annAccession + '/' + originalBed_accession + '_filteredSorted.bed'

    os.system('gunzip %s.gz' %(segwayFile))

    '''
    2. the chmm file

    '''
    chmmCount = 0
    for nonSelfAccession in annAccessionList:
        if nonSelfAccession == annAccession:
            continue

        if nonSelfAccession == 'ENCSR313QGL' or nonSelfAccession == 'ENCSR592IOP' or nonSelfAccession == 'ENCSR721USS' or nonSelfAccession == 'ENCSR273LUT' or nonSelfAccession == 'ENCSR699DMW':
            continue
    
        chmmFile = chmmFile_dict[nonSelfAccession]

        if chmmFile == 'none':
            continue

        print(chmmCount)
        chmmCount+=1

        #if chmmFile.endswith('.gz'):
        #    chmmFile = chmmFile[:-3]

        '''
        3. we have both chmmFile and segwayFile, so we can move on to comparison

        '''
    
        print(segwayFile)
        print(chmmFile)

        os.system('gunzip %s.gz' %(chmmFile))

        # unzip the chmm file
        # os.system('gunzip %s' %(chmmFile + '.gz'))

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
                    sindex = int(sfields[3].split('_')[0])
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

        outputFileName = 'overlap_segway_%s_chmm_%s_nonself.pkl' %(segwayAccession, chmmAccession)
        outputFile = dataFolder + dataSubFolder + annAccession + '/' + outputFileName
        with open(outputFile, 'wb') as f:
            pickle.dump(overlap_df, f)

    os.system('gzip %s' %(segwayFile))



