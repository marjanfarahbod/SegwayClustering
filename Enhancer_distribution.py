# This code is to get insight about the distribution of regions marked as enhancer throughout the genome.
# the plan is to just pick one file, one chromosome, and get some sense of this distribution. We will do chromosome 19.
# Items TODO:
# 1. Pick a sample (3, to catch the behavior)
# 2. Extract rows from chromosome 19 into another file
# 3. Load the mnemonics for the sample
# 4. Parse the file based on the algorithm. We have 1500 genes on chromosome 19.


import linecache
import pickle
import re
import numpy as np
import pandas as pd
import subprocess as sp
import seaborn as sns
import matplotlib.pyplot as plt
from QC_transcriptionComparison_util import Gene, Exon, Annotation, AnnotationClass


# open the folder and look up the stuff
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/testBatch105/'
runFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/runs/'

dataSubFolder = 'testBatch105/fromAPI/'

inputFile = dataFolder + dataSubFolder + 'metaInfo.pkl'
with open(inputFile, 'rb') as f:
    annMeta = pickle.load(f)

inputFileName = 'all_annInfo_list.pkl'
inputFile = dataFolder + dataSubFolder + inputFileName
with open(inputFile, 'rb') as f:
    ann_info_list = pickle.load(f)

# get the mapping from accession to index
accession_index_map = {}
for ann in ann_info_list:
   accession_index_map[ ann['accession'] ] = ann['index']

#classifier_input_file = runFolder + "model_296exp_reg0.067_auc0.77on.32test_classifier_labels.pickle"
#with open(classifier_input_file, "rb") as f:
    #classifier_labels = pickle.load(f)

#book = classifier_labels['cluster_list']
#interpretation_terms_set = set(classifier_labels['interpretation_term'])

#segwayLabels = ['Quiescent', 'ConstitutiveHet', 'FacultativeHet', 'Transcribed', 'Promoter', 'Enhancer', 'RegPermissive', 'Bivalent', 'LowConfidence']
segwayLabels = ['Enhancer_low', 'Enhancer', 'Promoter_flanking', 'Promoter', 'Transcribed', 'CTCF', 'K9K36', 'Bivalent', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent', 'Unclassified']


inputFile = dataFolder + dataSubFolder + 'biosample_tissue_info.pkl'
with open(inputFile, 'rb') as f:
    tissue_info = pickle.load( f)

'''
# GTF data structure
fileName = dataFolder + '/geneLists.pkl'
with open(fileName, 'rb') as pickleFile:
    geneListsAndIDs = pickle.load(pickleFile)

geneList = geneListsAndIDs[0]
geneIDList = geneListsAndIDs[1]
del geneListsAndIDs
'''

annAccessionList = list(annMeta.keys())
annAccession = annAccessionList[104]
print(annAccession)

annAccession = 'ENCSR566HBT'
index = accession_index_map[annAccession]
print(index)
    
sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
print(sampleFolderAdd)

# get the label from label from mnemonics: a dictionary from labels to terms
label_term_mapping = {}
mnemonics_file = sampleFolderAdd + 'mnemonics_v02.txt'
with open(mnemonics_file, 'r') as mnemonics:
    for line in mnemonics:
        #print(line)
        label = line.strip().split()[0]
        term = line.strip().split()[1]
        label_term_mapping[label] = term


# get the .bed file chr19 extract. 
originalBedFile = annMeta[annAccession]['bedFile']
originalBed_accession = originalBedFile.split('.')[0]
segwayFile = dataFolder + dataSubFolder + annAccession + '/' + originalBed_accession + '_filteredSorted.bed'
# I got the address of the .bed file from here and extracted the chromosome 19 annotation from it.


inputFile = dataFolder + 'testBatch105/' + 'runID_accession_map_105run.pkl'
with open(inputFile, 'rb') as f:
    runID_map = pickle.load(f)

accession_IDmap = {y: x for x, y in runID_map.items()}

print(accession_IDmap[annAccession])


fileAdd = dataFolder + dataSubFolder + annAccession + '/chr19.bed'
reg_list = []
flag = False
enCounter = 0
region_length = 3000 # it was 1000 before
with open(fileAdd, 'r') as f:
    for line in f:

        #if enCounter == 4:
        #    break

        items = line.split('\t')
        label = items[3].split('_')[0]
        term = label_term_mapping[label]
        region_s = int(items[1]) # just this line/item
        distance = 0
        region_dict = {}
        region_dict['start'] = region_s
        region_dict['length'] = int(items[2]) - region_s
        region_dict['label'] = label

        
        if term == 'Enhancer' and not(flag): # the first enhancer
            print(enCounter)
            enCounter += 1
            flag = True
            enhancer_region_list = []
            module_s = int(items[2]) # the whole region starts
            

        if flag:
            if (region_s - module_s) < region_length - 200: # read the region,
                #print(line)
                enhancer_region_list.append(region_dict)
            else:
                reg_list.append(enhancer_region_list)
                flag = False


# get the difference start reg lists:
reg_region_dis = np.zeros(len(reg_list))
first_region = reg_list[0]
s1 = first_region[0]['start']
for i in range(1,len(reg_list)):
    region = reg_list[i]
    s2 = region[0]['start']
    reg_region_dis[i-1] = s2-s1
    s1 = s2

    
plt.hist(np.log10(reg_region_dis+1), bins = 20)
plt.show()
'''
only 80 (3%) had less than 1500 distance from each other, which I consider safe distance as going through a new region
'''

# filling the reg_region_mat
colCount = int(region_length/100)
reg_region_mat = np.zeros((len(reg_list), colCount))
for i in range(len(reg_list)):
    reg_module = reg_list[i]
    resCount = 0
    for region in reg_module:
        region_len = region['length']
        region_label = region['label']
        region_resCount = min(int(region_len/100), colCount-resCount) # this was 10-resCount, I introduced colCount to change region length
        for j in range(region_resCount):
            reg_region_mat[i, j+resCount] = region_label

        resCount = resCount + region_resCount

#### MAKING THE TREE
#order_mat = np.zeros((120,8))
#remain = 0
#final_order = np.zeros(len(reg_region_mat))
maxNodeCount = 250
ratio_thr = int(.01*2664)+1
#layer, label, parent, parentind, count
branches = np.zeros((5, maxNodeCount)) - 1
branches[:, 0] = [0, 7, -1, -1, 2664]
walked = np.zeros(maxNodeCount) - 1
nci = 1 # node walk index, walks on - root is already filled
layer = 1 # initial - root is zero
parent = 0 # initial - for layer 1
parentind = 0 # initial - for layer 1
#walkedCount = 0 I think I don't need these, criteria should be if there are any 0 in walked: nodes generated but not fully walked
#nodeCount = 0
node_indsLists = {}
node_excess_indLists = {}
child_indsLists = {}
visited_leaf = []
walked_node = []
select_label_ind = 0
max_layer = colCount - 3 # it was 7 before (wwhen I had 1000 bp region, now I have 17)
while sum(walked == 0) > 0 or sum(walked==-1)==maxNodeCount: # while there are node generated, but not walked

    #if walked[11] == 1:
    #    break

    # get the counts and creating nodes: 
    if sum(walked==0)==0: # if we are at the root: have not visited or walked any node yet 
        temp_indsList = {}
        counts = np.zeros((16))
        excess_list = np.array([])
        for i in range(16):
            temp_inds = np.where(reg_region_mat[:,2] == i)[0] # where returns a tuple
            count = sum(reg_region_mat[:,2]==i)
            if count >= ratio_thr:
                counts[i] = count
                temp_indsList[i] = temp_inds
            else:
                excess_list = np.concatenate([temp_inds, excess_list])
            
                
    else: # we are not in the root, and we are on a path from some node

        parentind = select_label_ind # the parent of the label nodes generated in this iteration
        parent = branches[1, select_label_ind] # the parent vlue of the current label sets

        # TODO: do we need to fix the layer here
        layer = branches[0, select_label_ind] + 1 # the layer of the nodes generated in this iteration
        
        filter_ind = int(layer) + 2 # which index of the label matrix should we look at in this iteration
        # list of indices in the mat we are going to work with:
        ind_list = node_indsLists[select_label_ind] # list of the indices that we used for the filter

        # get the sub_mat from the whole mat
        subMat = reg_region_mat[ind_list, ]
        
        # do the counts for the sub mat. Remember to do the indexing from the indexing.
        temp_indsList = {}
        counts = np.zeros((16))
        excess_list = np.array([])
        for i in range(16):
            temp_inds = ind_list[np.where(subMat[:, filter_ind] == i)[0]]
            count = len(temp_inds)
            if count >= ratio_thr:
                counts[i] = count
                temp_indsList[i] = temp_inds
            else:
                excess_list = np.concatenate([temp_inds, excess_list])

    node_excess_indLists[select_label_ind] = excess_list
    sorted_labels = np.argsort(counts)[::-1][0:sum(counts>0)]

    if len(sorted_labels) == 0 or layer == max_layer: # go up the tree

        # mark the current node as walked,
        # walked[select_label_ind] = 1 # this node has no child 
        # parent = branches[2, select_label_ind]

        # leaf_visited
        visited_leaf.append(select_label_ind)

        # return to a layer above, pick the next one to go down
        inds = np.array([]) # inds of unvisited siblings

        current_ind = select_label_ind
        while len(inds)==0 and layer >-1: 
            walked[current_ind] = 1 # none of the childs are remaining, so it is walked
            walked_node.append(current_ind)
            print('walked %d' %(current_ind))
            current_ind = int(branches[3, current_ind]) # going one layer up to the parent's node. 
            layer = branches[0, current_ind] # going one layer up, to the parent's layer
            #label = branches[1,:]
            #label_inds = (branches[0,:] == layer)
            child_inds = (branches[3,:] == current_ind) # indices of siblings
            sib = np.logical_and(walked==0, child_inds == True) # indices of parent's children that were not walked
            inds = np.where(sib==True)[0] # the first element of the tuple

        if layer > -1:
            select_label_ind = np.min(inds) # the label that hasn't been walked,

    else: # build the tree and go down, add the child and fill the data 
        child_list = []
        for label in sorted_labels:
            #layer, label, parent, count
            branches[0, nci] = layer
            branches[1, nci] = label
            branches[2, nci] = parent
            branches[3, nci] = parentind
            branches[4, nci] = counts[label]
            node_indsLists[nci] = temp_indsList[label]
            walked[nci] = 0
            child_list.append(nci)
            
            nci += 1
            print(nci)
            # node_excess_indLists[layer] = excess_list # this is wrong

        # from the branch, in this layer, pick one with the most count that hasn't been walked, to continue to walk on
        label_inds = (branches[0,:] == layer)
        sib = np.logical_and(walked==0, label_inds==True)
        inds = np.where(sib==True)

        select_label_ind = np.min(inds) # the label that hasn't been walked,

        child_indsLists[parentind] = child_list


final_inds = np.array([]).astype(int)
for node in walked_node:
    print('walked_node_nci %d' %(node))
    if node in visited_leaf:
        print('leaf')
        all_inds = set(node_indsLists[node])
        excess_inds = set(node_excess_indLists[node])
        final_inds = np.concatenate([final_inds, np.array(list(all_inds - excess_inds))])
        final_inds = np.concatenate([final_inds, node_excess_indLists[node]])
        print(len(final_inds))
    else:
        print('node')
        final_inds = np.concatenate([final_inds, node_excess_indLists[node]])
        print(len(final_inds))

    
kado = np.unique(final_inds)
book = final_inds.astype(int)
sib = reg_region_mat[book, :]

label_term_mapping

label_order = [10, 14, 2, 7, 9, 13, 11, 12, 6, 4, 15, 0, 3, 8, 1, 5] # sample 02
label_color_mapping = {}
label_color_mapping[10] = [255, 0, 0] # promoter
label_color_mapping[14] = [255, 69, 0] # promoter
label_color_mapping[2] = [220, 210, 0] # enhancer
label_color_mapping[7] = [195, 180, 0] # enhancer
label_color_mapping[9] = [255, 255, 0] # enhancer
label_color_mapping[13] = [180, 160, 0] # enhancer
label_color_mapping[11] = [102, 205, 170] # CTCF
label_color_mapping[12] = [0, 128, 0] # Transcribed
label_color_mapping[6] = [0, 170, 0] # Transcribed
label_color_mapping[4] = [200, 10, 200] # FacHet
label_color_mapping[15] = [150, 5, 150] # FacHet
label_color_mapping[0] = [138, 145, 208] # consHet
label_color_mapping[3] = [180, 195, 250] # consHet
label_color_mapping[8] = [100, 120, 180] # consHet
label_color_mapping[1] = [240, 240, 240] # Quis
label_color_mapping[5] = [230, 230, 230] # Quis

label_order = [2,6,12,13,14, 4,3,7,9,11, 15,0, 10,1,5,8] # sample 104
label_color_mapping = {}
label_color_mapping[2] = [255, 0, 0] # promoter
label_color_mapping[6] =  [180, 160, 0] # enhancer
label_color_mapping[12] = [220, 210, 0] # enhancer
label_color_mapping[13] = [195, 180, 0] # enhancer
label_color_mapping[14] = [255, 255, 0] # enhancer
label_color_mapping[4] = [0, 170, 0]# transcribed
label_color_mapping[3] =  [0, 150, 0]# Transcribed
label_color_mapping[7] = [0, 128, 0] # Transcribed
label_color_mapping[9] =  [220, 14, 220]# FacHet
label_color_mapping[11] = [200, 10, 200] # FacHet
label_color_mapping[15] = [150, 5, 150] # FacHet
label_color_mapping[0] = [138, 145, 208] # consHet
label_color_mapping[10] = [180, 195, 250] # consHet
label_color_mapping[1] = [200, 200, 200] # Quis
label_color_mapping[5] = [240, 240, 240] # Quis
label_color_mapping[8] = [230, 230, 230] # Quis
# I will make the mapping myself here as I want each of them to have a different color

# making the colormap:
colormap = label_color_mapping[label_order[0]]
colorlist = []
colorlist.append(np.array(label_color_mapping[label_order[0]])/255)
for label in label_order[1:]:
    colormap = np.vstack([colormap, label_color_mapping[label]])
    colorlist.append(np.array(label_color_mapping[label])/255)

import matplotlib.colors as colors
cmap = colors.ListedColormap(colorlist)
boundaries = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
norm = colors.BoundaryNorm(boundaries, cmap.N, clip=True)

# we have to do a change in sib.
sibplus = np.copy(sib)

for i,label in reversed(list(enumerate(label_order))):
    book = np.where(sib == label)
    sibplus[book] = i


sns.heatmap(sibplus, cmap=cmap, norm=norm)
plt.show()

kado = np.zeros(sibplus.shape)
for i in range(30):
    halva = np.sort(sibplus[:, i])
    kado[:,i] = halva

sns.heatmap(kado, cmap=cmap, norm=norm)
plt.show()



# for the 30 columns, plot the coverage of different labels.
covmat = np.zeros((16, 30))
# for each column, for each label
for i in range(30):
    for j in range(16):
        covmat[j, i] = sum(sibplus[:,i]==j)

# do a line plot sth.
print(np.sum(covmat, axis=0))

covmat = covmat/2044
for i in range(16):
    
    plt.plot(covmat[i,2:], color = cmap.colors[i])
plt.show()

sns.heatmap(covmat[:, 3:])
plt.show()

covmat = covmat/200

fig, ax = plt.subplots()

for i in range(16):
    ax.bar(range(30), covmat[i,:], color = cmap.colors[i])
    print(i)

plt.show()

plt.close('all')



### TODO: read the tree and make the matrix
# TODO: add the child inds 
branches 
node_indsLists 
layer_excess_indLists
child_indsList
reg_region_mat.shape



for i in range(14):
    count = sum(reg_region_mat[:,2]==i)
    if count/2664 > .05:
        print('%d: %d' %(i,count))

inds = reg_region_mat[:,2]==7
sib = reg_region_mat[inds, :]

for i in range(14):
    count = sum(sib[:,3]==i)
    if count/2664 > .01:
        print('*%d: %d' %(i,count))
    else:
        print('%d: %d' %(i,count))


# mygoald is to extract information regarding the Enhancer regions. My prior is that our labels for Enahncer happen in clusters and are interrupted by other relevant or non relevant (to be examined) labels. For this sample, I have 

with open(fileAdd, 'r') as f:
    


for annAccession in annAccessionList:

    index = accession_index_map[annAccession]
    
    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
    print(sampleFolderAdd)

    # get the label from label from mnemonics: a dictionary from labels to terms
    label_term_mapping = {}
    mnemonics_file = sampleFolderAdd + 'mnemonics_v02.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term




