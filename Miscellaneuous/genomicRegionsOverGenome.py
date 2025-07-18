# this is a code to record all genomic regions around and on gene bodies. The reason I need this code (instead of just getting regions around the genes) is that I want to count overlapping genes once only. 

# the next step after this one, I want to record what percentage of the enhancers identified in annotations happen in these regions. Region 1: genomic regions: gene body + 2, 5 and 10kb around the gene. 

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
import pyBigWig

import glob
 
# General data folder
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/'
figureFolder = 'figure02/'


inputFileName = 'all235Annot_meta_corrected.pkl'
#inputFileName = 'all235Annot_meta.pkl'
inputFile = dataFolder +  inputFileName
with open(inputFile, 'rb') as f:
    allMeta = pickle.load(f)

# GTF data structure
fileName = dataFolder + '/geneLists.pkl'
with open(fileName, 'rb') as pickleFile:
    geneListsAndIDs = pickle.load(pickleFile)

book =  'EnhA1 EnhA2 EnhBiv EnhG1 EnhG2 EnhWk Het Quies ReprPC ReprPCWk TssA TssBiv TssFlnk TssFlnkD TssFlnkU Tx TxWk ZNF/Rpts'
chromLabels = book.split()

segwayStates = ['Promoter', 'PromoterFlanking', 'Enhancer', 'EnhancerLow', 'Bivalent', 'CTCF', 'Transcribed', 'K9K36', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent']
print(len(segwayStates))
segwayStateCount = len(segwayStates)

inputFileName = 'all235Annot_meta_corrected.pkl'
#inputFileName = 'all235Annot_meta.pkl'
inputFile = dataFolder +  inputFileName
with open(inputFile, 'rb') as f:
    allMeta = pickle.load(f)

accessionList = list(allMeta.keys())
accessionList.remove(accessionList[205])

# GTF data structure
fileName = dataFolder + '/geneLists.pkl'
with open(fileName, 'rb') as pickleFile:
    geneListsAndIDs = pickle.load(pickleFile)

geneList = geneListsAndIDs[0]
geneIDList = geneListsAndIDs[1]
del geneListsAndIDs

extensions = [0, 2000, 5000, 10000]
allRegionLists = {}
totalgenebody = 0
for ext in extensions:
    
    previous_gene_start = 0
    previous_gene_end = 0
    previous_gene_chr = 'chr'
    regionList = []
    regionInd = 0
    for cgi in range(len(geneList)):
    
        geneID = geneIDList[cgi]
        gene_chr = geneList[geneIDList[cgi]].chrom

        gene_start = np.max([0, geneList[geneID].start - ext])
        gene_end = geneList[geneID].end + ext
        totalgenebody += gene_end - gene_start

        if gene_start < previous_gene_end and gene_chr == previous_gene_chr:
            previous_gene_start = regionList[regionInd-1][1]
            regionList[regionInd -1] = (gene_chr, previous_gene_start, gene_end)
            #print('_____', cgi)
            #print(regionInd)
            #print('here')
            #break
            continue

        regionList.append((gene_chr, gene_start, gene_end))
        previous_gene_chr = gene_chr
        previous_gene_end = gene_end
        cgi+=1
        regionInd +=1

    allRegionLists[ext] = regionList
    

file = dataFolder + 'genomicRegionCoverage.pkl'
with open(file, 'wb') as f:
    pickle.dump(allRegionLists, f)

# obtain total coverage of these regions
totalbp = {}
for ext in extensions:
    tbp = 0
    myRegions = allRegionLists[ext]
    for region in myRegions:
        tbp += region[2] - region[1]

    totalbp[ext] = tbp

geneCounts = {}
for ext in extensions:
    myRegions = allRegionLists[ext]
    for r, region in enumerate(myRegions):
        if region[0] == 'chrX':
            print(region, r)
            geneCounts[ext] = r
            break

# almost half of the 3bp regions is covered when we extend the genomic region by 10kb.
for accession in accessionList[183:]:
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']
    annFile = annotation['bedFile']
    mnemFile = annotationFolder + 'mnemonics_v04.txt'
    print(accession)

    # >>>>>>>>>> Segway, ev = 2000
    for ext in extensions[1:]:
        regionList = allRegionLists[ext][0:geneCounts[ext]]
        annFile = annotation['bedFile']
        mnemFile = annotationFolder + 'mnemonics_v04.txt'
        Segway_label_coverage(annotationFolder,
                                           annFile, mnemFile, regionList, ext)


def Segway_label_coverage(annotationFolder, annFile, mnemFile, regionList, ext):

    totalGeneCount = len(regionList)
    chrIndex = {'chr1':1, 'chr2':2, 'chr3':3, 'chr4':4, 'chr5':5, 'chr6':6, 'chr7':7, 'chr8':8, 'chr9':9,
                'chr10':10, 'chr11':11, 'chr12':12, 'chr13':13, 'chr14':14, 'chr15':15, 'chr16':16,
                'chr17':17, 'chr18':18, 'chr19':19, 'chr20':20, 'chr21':21, 'chr22':22}

    # prepare ann
    if annFile.endswith('.gz'):
        os.system('gunzip %s' %(annFile))
    else:
        os.system('gunzip %s.gz' %(annFile))

    print(annFile)

    annLineCount = int(sp.getoutput('wc -l %s' %(annFile)).split()[0])

    # get the mnemonics
    label_term_mapping = {}
    with open(mnemFile, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term

    clusterCount = len(label_term_mapping)
    sortedClusterList = list(range(0,clusterCount))

    print(sortedClusterList)

    cgi = 0 # walks on the geneIDList
    ann_start = 0 # the start of the current annotation
    ann_end = 0 # the end of the current annotation
    ann_chr = 'chr'
    ann_line_count = 0 # this is just to check the progress through the annotation file
    previous_class = ''

    #labelExpMat = np.zeros((len(geneIDList), clusterCount))
    #labelPromoterMat = np.zeros((len(geneIDList), clusterCount))

    labelExpMat = np.zeros((totalGeneCount, clusterCount))
        
    annLineInd = 1 # this keeps the annotation line index, the first line is empty, thus the index starts at 1
    genomic_region_start_annLineInd = 0
    
    previous_gene_chr = 'chr1'
    previous_extension_end = 0 # doesn't matter since chr condition is never true until the first gene is processed
    previous_ann_chr = 'chr1'

    while cgi <  totalGeneCount:# (count of genes before chrX and chrY) #23897:#21468: #26017: # >>>>>>>>>> MAIN

        'we are in the gene territory, we are walking on the genes'
        
        gene_start = regionList[cgi][1]
        gene_end = regionList[cgi][2]
        gene_chr = regionList[cgi][0]

        gene_length = gene_end - gene_start
        gene_coverage = 0

        # reading the next annotation
        line = linecache.getline(annFile, annLineInd)
        annLineInd +=1
        ann_line_count += 1
        fields = line.strip().split()
            
        previous_ann_chr = ann_chr
        ann_chr = fields[0]
        ann_start = int(fields[1])
        ann_end = int(fields[2])

        #>>>>>>>>>>>>>>>>>>>> the fix code with chromosome index. If gene moved to the next, walk on in the ann. If ann moved to the next, walk down in the ann

        if (chrIndex[gene_chr] > chrIndex[ann_chr]): # in case of chromosome change because gene moved to the next chromosome
            print('gene chr not equal to previous gene chr')
            print('gene_chr %s' %(gene_chr))
            print('p_gene_chr %s' %(previous_gene_chr))
            print('ann_chr %s' %(ann_chr))
                
            while (ann_chr != gene_chr): # move on in the annotation until we reach to the next chromosome
                #print('reading annotations until it is equal') # >>>> test
                line = linecache.getline(annFile, annLineInd)
                fields = line.strip().split()
                annLineInd +=1
                ann_line_count += 1

                previous_ann_chr = ann_chr                
                ann_chr = fields[0]
                ann_start = int(fields[1])
                ann_end = int(fields[2])

            print(ann_chr)
            print(previous_ann_chr)

                
        if (chrIndex[gene_chr] < chrIndex[ann_chr]): # if annotation moved to the next chromosome, but gene has not yet moved to the next chromosome
            annLineInd = annLineInd - 2
            line = linecache.getline(annFile, annLineInd)
            annLineInd +=1
            ann_line_count += 1
            fields = line.strip().split()

            ann_chr = fields[0]
            ann_start = int(fields[1])
            ann_end = int(fields[2])
            print('marker 01')
        

        while (ann_start > gene_start) and (gene_chr == ann_chr): # in case of overlapping genes
                
            annLineInd = max(annLineInd - 2,1)
            line = linecache.getline(annFile, annLineInd)
            annLineInd +=1
            ann_line_count += 1
            fields = line.strip().split()

            ann_chr = fields[0]
            ann_start = int(fields[1])
            ann_end = int(fields[2])

        while ((ann_start < gene_end) and (gene_chr == ann_chr)) and gene_coverage < gene_length: #and geneMatWalkIndex < 160: # while we still have annotation before the gene

            # in the next sections we are processing the annotation (if there is overlap)
            if not((ann_end < gene_start) or ann_start > gene_end):

                ''' We are in the genonimc region (with extension)'''
                ''' Taking off for filling the matrices... '''
                        
                adjusted_ann_start = max(0, ann_start - gene_start)
                adjusted_ann_end = min(ann_end - gene_start, gene_end - gene_start)
                adjusted_ann_length = adjusted_ann_end - adjusted_ann_start

                ann_cluster = fields[3].split('_')[0]
                clusterInd = int(ann_cluster)

                if gene_coverage < gene_length and adjusted_ann_length > 0:
                    new_coverage = min(adjusted_ann_length, gene_length - gene_coverage)
                    labelExpMat[cgi, clusterInd] += new_coverage
                    adjusted_ann_length = adjusted_ann_length - new_coverage
                    gene_coverage += new_coverage
                    #coverPromoter
 

            # if there is no overlap
            if gene_coverage < gene_length:
                line = linecache.getline(annFile, annLineInd)
                annLineInd +=1

                ann_line_count += 1
                fields = line.strip().split()

                ann_chr = fields[0]
                ann_start = int(fields[1])
                ann_end = int(fields[2])
                #print(ann_start)
                #print(ann_end)

        if np.sum(labelExpMat[cgi, :]) == 0 and not(gene_length <= 0):
            print('------------- zero')
            print(cgi)
            break
        
        cgi += 1 # next gene
        previous_gene_end = gene_end
        previous_gene_chr = gene_chr

    linecache.clearcache()
    os.system('gzip %s' %(annFile))

    myOutPut = labelExpMat
    outputFile = annotationFolder  + 'Segway_genomicRegionLabelCoverage_%d.pkl' %(ext)
    with open(outputFile, 'wb') as f:
        pickle.dump(myOutPut, f)

    print(outputFile)

# >>>>> Analysis: now let's get the coverage of our enhancer labels - and the coverage of fantom5 enhancers?

# Get the count of bps in chr1-22 for our annotations.

# Get the coverage of each label through the genome. 

# get the coverage of the enhancer label within our regions.

# >>>>> Analysis: get the coverage of Fantom5 enhancers - how much do they overlap?



