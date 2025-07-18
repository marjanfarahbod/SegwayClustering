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

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/'
figureFolder = 'figure01/'

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
print(len(segwayStates))
segwayStateCount = len(segwayStates)

trackList = ['H3K36me3', 'H3K27me3', 'H3K4me3', 'H3K9me3', 'H3K27ac', 'H3K4me1', 'CTCF', 'DNase-seq', 'ATAC-seq', 'POLR2A', 'EP300']

chrList = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
           'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']


snpFile = dataFolder + 'trait_snp_linkage_for_marjan.csv'

snpLineCount = int(sp.getoutput('wc -l %s' %(snpFile)).split()[0])

line = linecache.getline(snpFile, 12345)
snpCount = 0
snpList = []
for i in range(snpLineCount):
    line = linecache.getline(snpFile, i)
    trait = line.split(',')[0]

    if trait == 'blood metabolite levels':
        snp = line.split(',')[1]
        snpList.append(snp)
        snpCount += 1

# find the main accession
for accession in accessionList:
    annotation = allMeta[accession]
    print(annotation['tissueInfo'])
    if (annotation['tissueInfo'][0] == 'Liver tissue male adult (31 years)'):
        print(accession)
    
accession = 'ENCSR111ABE'

# define the matrix for collecting
chrList = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
           'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']


annotation = allMeta[accession]
annotationFolder = annotation['folder']
annFile = annotation['bedFile']

if annFile.endswith('gz'):
    os.system('gunzip %s' %(annFile))
    annFile = annFile[0:-3]
else:
    os.system('gunzip %s.gz' %(annFile))
    
annLineCount = int(sp.getoutput('wc -l %s' %(annFile)).split()[0])
    
mnemFile = annotationFolder + 'mnemonics_v04.txt'
    
label_term_mapping = {}
with open(mnemFile, 'r') as mnemonics:
    for line in mnemonics:
        #print(line)
        label = line.strip().split()[0]
        term = line.strip().split()[1]
        label_term_mapping[label] = term

W = 200
snpMat = np.zeros((snpCount, W)) # 20k around snp
for s,snp in enumerate(snpList):

    print('>>>>>>>>>>>>>>>')
    regionChr = snp.split('_')[1]
    regionChrInd = chrList.index(regionChr)
    scor = int(snp.strip().split('_')[2]) - (W*50)
    
    slineInd = 1
    elineInd = annLineCount
    lineInd = int(elineInd/2)
    line = linecache.getline(annFile, lineInd)
    chr = line.split()[0]
    while chrList.index(chr) != regionChrInd:
    
        if chrList.index(chr) < regionChrInd:
            slineInd = lineInd
            lineInd = int(slineInd + ((elineInd-slineInd)/2))
            line = linecache.getline(annFile, lineInd)
            chr = line.split()[0]
            print(chr)
        else:
            if chrList.index(chr) > regionChrInd:
                elineInd = lineInd
                lineInd = int(slineInd + ((elineInd-slineInd)/2))
                line = linecache.getline(annFile, lineInd)
                chr = line.split()[0]

    sind = int(line.split()[1])
    if sind > scor: # staring coordinate
        eind = int(line.split()[2])
        while (sind > scor) or (eind < scor):
            lineInd = lineInd - 1
            line = linecache.getline(annFile, lineInd)
            sind = int(line.split()[1])
            eind = int(line.split()[2])
    else:
        if sind < scor:
            eind = int(line.split()[2])
            while (sind > scor) or (eind < scor):
                lineInd = lineInd + 1
                line = linecache.getline(annFile, lineInd)
                sind = int(line.split()[1])
                eind = int(line.split()[2])
            

    walker = 0
    steps = int((eind - scor)/100)
    color = label_term_mapping[line.split()[3].split('_')[0]]
    colorInd = segwayStates.index(color)
    snpMat[s, walker:walker+steps-1] = colorInd
    walker += steps
    print(steps)
    while walker < W:
        print(line)
        print(steps)
        lineInd +=1
        print(walker)
        print(color)
        line = linecache.getline(annFile, lineInd)
        sind = int(line.split()[1])
        eind = int(line.split()[2])
        steps = min(int((eind-sind)/100), W-walker)
        color = label_term_mapping[line.split()[3].split('_')[0]]
        colorInd = segwayStates.index(color)
        print(colorInd)
        snpMat[s, walker:walker+steps] = colorInd
        walker += steps
        
linecache.clearcache()
os.system('gzip %s' %(annFile))


mycolors = [[255,0,0],
          [255,68,0],
          [255,195,77],
          [255,255,0],
          [189,183,107],
          [196,225,5],
          [0,128,0],
          [102,205,170],
          [128,0,128],
          [138,145,208],
          [255,255,255]]

colorList = []
for color in mycolors:
    colorList.append(np.array(color)/256)

import matplotlib.colors as colors
cmap = colors.ListedColormap(colorList)
boundaries = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
norm = colors.BoundaryNorm(boundaries, cmap.N, clip=True)


df = pd.DataFrame(snpMat)

sib = sns.clustermap(df, method = 'ward', col_cluster=False)
rowInds_cluster = sib.dendrogram_row.reordered_ind

sns.heatmap(snpMat[rowInds,:], cmap=cmap, norm=norm)
plt.show()

copySnpMat = np.copy(snpMat)

for i in range(snpCount):
    copySnpMat[i,:] = np.sort(copySnpMat[i,:])

df = pd.DataFrame(copySnpMat)

sib = sns.clustermap(df, method = 'ward', col_cluster=False)
rowInds_sortCluster = sib.dendrogram_row.reordered_ind

sns.heatmap(snpMat[rowInds_sortCluster,:], cmap=cmap, norm=norm)
sns.heatmap(snpMat[rowInds_sortCluster,:], cmap=cmap, norm=norm)
plt.title('snps for blood metabolite levels')
figFile = plotFolder + figureFolder + 'MAIN_%s_%s_sortClusterluster.pdf' %(tissue, accession)
plt.savefig(figFile)
plt.close('all')

plt.show()

########################################
figureFolder = 'liver_SNP_regionAnnotation/'

for accession in accessionList[1:]:

    print(accessionList.index(accession), accession)
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']
    annFile = annotation['bedFile']
    tissue = annotation['tissueInfo'][0]

    if annFile.endswith('gz'):
        os.system('gunzip %s' %(annFile))
        annFile = annFile[0:-3]
    else:
        os.system('gunzip %s.gz' %(annFile))
    
    annLineCount = int(sp.getoutput('wc -l %s' %(annFile)).split()[0])
    
    mnemFile = annotationFolder + 'mnemonics_v04.txt'
    
    label_term_mapping = {}
    with open(mnemFile, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term

    W = 200
    snpMat = np.zeros((snpCount, W)) # 20k around snp
    for s,snp in enumerate(snpList):

        #print('>>>>>>>>>>>>>>>')
        regionChr = snp.split('_')[1]
        regionChrInd = chrList.index(regionChr)
        scor = int(snp.strip().split('_')[2]) - (W*50)
    
        slineInd = 1
        elineInd = annLineCount
        lineInd = int(elineInd/2)
        line = linecache.getline(annFile, lineInd)
        chr = line.split()[0]
        while chrList.index(chr) != regionChrInd:
    
            if chrList.index(chr) < regionChrInd:
                slineInd = lineInd
                lineInd = int(slineInd + ((elineInd-slineInd)/2))
                line = linecache.getline(annFile, lineInd)
                chr = line.split()[0]
                #print(chr)
            else:
                if chrList.index(chr) > regionChrInd:
                    elineInd = lineInd
                    lineInd = int(slineInd + ((elineInd-slineInd)/2))
                    line = linecache.getline(annFile, lineInd)
                    chr = line.split()[0]

        sind = int(line.split()[1])
        if sind > scor: # staring coordinate
            eind = int(line.split()[2])
            while (sind > scor) or (eind < scor):
                lineInd = lineInd - 1
                line = linecache.getline(annFile, lineInd)
                sind = int(line.split()[1])
                eind = int(line.split()[2])
        else:
            if sind < scor:
                eind = int(line.split()[2])
                while (sind > scor) or (eind < scor):
                    lineInd = lineInd + 1
                    line = linecache.getline(annFile, lineInd)
                    sind = int(line.split()[1])
                    eind = int(line.split()[2])
                    

        walker = 0
        steps = int((eind - scor)/100)
        color = label_term_mapping[line.split()[3].split('_')[0]]
        colorInd = segwayStates.index(color)
        snpMat[s, walker:walker+steps-1] = colorInd
        walker += steps
        #print(steps)
        while walker < W:
            #print(line)
            #print(steps)
            lineInd +=1
            #print(walker)
            #print(color)
            line = linecache.getline(annFile, lineInd)
            sind = int(line.split()[1])
            eind = int(line.split()[2])
            steps = min(int((eind-sind)/100), W-walker)
            color = label_term_mapping[line.split()[3].split('_')[0]]
            colorInd = segwayStates.index(color)
            #print(colorInd)
            snpMat[s, walker:walker+steps] = colorInd
            walker += steps
        
    linecache.clearcache()
    os.system('gzip %s' %(annFile))

    sns.heatmap(snpMat[rowInds_sortCluster,:], cmap=cmap, norm=norm)
    plt.title('snps for blood metabolite levels')
    figFile = plotFolder + figureFolder + '%s_%s_sortCluster.pdf' %(tissue, accession)
    plt.savefig(figFile)
    plt.close('all')

    sns.heatmap(snpMat[rowInds_cluster,:], cmap=cmap, norm=norm)
    plt.title('snps for blood metabolite levels')
    figFile = plotFolder + figureFolder + '%s_%s_cluster.pdf' %(tissue, accession)
    plt.savefig(figFile)
    plt.close('all')

    copySnpMat = np.copy(snpMat)
    for i in range(snpCount):
        copySnpMat[i,:] = np.sort(copySnpMat[i,:])

    df = pd.DataFrame(copySnpMat)

    sib = sns.clustermap(df, method = 'ward', col_cluster=False)
    rowInds_sortCluster_self = sib.dendrogram_row.reordered_ind
    plt.close('all')
    sns.heatmap(snpMat[rowInds_sortCluster_self,:], cmap=cmap, norm=norm)
    plt.title('snps for blood metabolite levels')
    figFile = plotFolder + figureFolder + '%s_%s_selfSortCluster.pdf' %(tissue, accession)
    plt.savefig(figFile)
    plt.close('all')
    

# in main sample:
# for THE trait, for all the snps, get the heatmap of the regional annotation
# cluster the heatmap
# record the clustering index
# get the same plot for the rest of the samples

