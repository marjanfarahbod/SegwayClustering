#########################################
import linecache
import pickle
import re
import numpy as np
import pandas as pd
import os
import shutil
import subprocess


dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
dataSubFolder = 'testBatch105/fromAPI/'

inputFile = dataFolder + dataSubFolder + 'metaInfo.pkl'
with open(inputFile, 'rb') as f:
    annMeta = pickle.load(f)

sample_count = len(annMeta)

annAccessionList = list(annMeta.keys())
annAccession = annAccessionList[104]
print(annAccession)

inputFileName = 'all_annInfo_list.pkl'
inputFile = dataFolder + dataSubFolder + inputFileName
with open(inputFile, 'rb') as f:
    ann_info_list = pickle.load(f)

# get the mapping from accession to index
accession_index_map = {}
for ann in ann_info_list:
   accession_index_map[ ann['accession'] ] = ann['index']

# just get chr19 file
# create the en file
for annAccession in annAccessionList:

    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
    print(sampleFolderAdd)

    index = accession_index_map[annAccession]
    print(index)

    # get the .bed file chr19 extract. 
    originalBedFile = annMeta[annAccession]['bedFile']
    originalBed_accession = originalBedFile.split('.')[0]
    segwayFile = dataFolder + dataSubFolder + annAccession + '/' + originalBed_accession + '_filteredSorted.bed.gz'

    os.system('gunzip %s' %(segwayFile))

    command = "grep -E 'chr19.*' %s" %(segwayFile[0:-3])
    print(command)
    out = sampleFolderAdd + 'chr19.bed'
    f = open(out, 'w')
    subprocess.run(command, shell=True, stdout=f)

    os.system('gzip %s' %(segwayFile[0:-3]))

    print('chr19 copied')


for annAccession in annAccessionList:

    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
    print(sampleFolderAdd)

    index = accession_index_map[annAccession]
    print(index)
    # load the mnomics, and extract the enhancer labels for mnemonics
    label_term_mapping = {}
    mnemonics_file = sampleFolderAdd + 'mnemonics_v02.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.split()[0]
            term = line.split()[1]
            label_term_mapping[label] = term

    chr19_file = sampleFolderAdd + 'chr19.bed'
    browser_file = sampleFolderAdd + 'chr19_updated.bed'
    with open(chr19_file, 'r') as old, open(browser_file, 'w') as updated:
        for line in old:
            fields = line.split()
            label = fields[3].split('_')[0]
            term = label_term_mapping[label]
            f3 = label + '_' + term
            updated.write('%s\t%s\t%s\t%s\n' %(fields[0], (fields[1]), (fields[2]), f3))

        
        

