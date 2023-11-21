# check the code SNP_regionAnnotation.py

import gzip
import linecache
import pickle
import os
import numpy as np
import pandas as pd
import subprocess as sp
import seaborn as sns
import matplotlib.pyplot as plt


# General data folder
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/'

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
for accession in accessionList:
    bedFile = allMeta[accession]['bedFile']

    myFile = bedFile + '.gz'
    with gzip.open(myFile, 'rb') as myzip:
            #ml = myzip.readline().decode()
            ml = myzip.readlines()


########################################
# find the tissue-trait with high level of association
########################################

snpFile = dataFolder + 'trait_snp_linkage_for_marjan.csv'

snpLineCount = int(sp.getoutput('wc -l %s' %(snpFile)).split()[0])

line = linecache.getline(snpFile, 12345)
line = linecache.getline(snpFile, 848)
snpCount = 0
snpList = []
for i in range(snpLineCount):
    line = linecache.getline(snpFile, i)
    trait = line.split(',')[0]

    #if trait == 'blood metabolite levels':
    #if trait == 'cholesterol to total lipids ratio in small hdl':
    if trait == 'eczema':
        snp = line.split(',')[1]
        snpList.append(snp)
        snpCount += 1

# find the accession in 
accession = 'ENCSR111ABE' # liver?
accession = 'ENCSR900BXG' # liver tissue female adult 25
accession = 'ENCSR372MII' # CD4 positive, alpha, beta T cell
accession = 'ENCSR607YQA' # CD4 positive, alpha beta memory Tcell

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

    print('>>>>>>>>>>>>>>>', s)
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

sns.heatmap(snpMat[rowInds_cluster,:], cmap=cmap, norm=norm)
plt.show()

### get the count of enhancer or enhancer low in the snpMat

enhancerFilter = (snpMat==2).astype(int) + (snpMat==3).astype(int)
promoFilter = (snpMat==0).astype(int) + (snpMat==1).astype(int)
transFilter = (snpMat == 6).astype(int)

mainSnpMat = np.copy(snpMat)

print(np.sum(enhancerFilter)/ (snpMat.shape[0] * snpMat.shape[1]))
print(np.sum(promoFilter)/ (snpMat.shape[0] * snpMat.shape[1]))
print(np.sum(transFilter)/ (snpMat.shape[0] * snpMat.shape[1]))
# let's get these for other tissues and count the enhancer, promoter and transcribed coverage.

enhOverlap = np.zeros((len(accessionList)))
#promOverlap = np.zeros((len(accessionList), 1))
transOverlap = np.zeros((len(accessionList)))
for a, accession in enumerate(accessionList):
#for accession in accessionList:
#    a  = accessionList.index(accession)

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

    print('here')
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

        print('>>>>>>>>>>>>>>>', s)
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
        chrLineInd = lineInd
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
                    lineInd = lineInd +1
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

    print('here')
    linecache.clearcache()
    os.system('gzip %s' %(annFile))

    thisEnh = ((snpMat==2).astype(int) + (snpMat==3).astype(int))
    enhOverlap[a] = np.sum((np.multiply(thisEnh, enhancerFilter) > 0).astype(int))
    #promOverlap = np.zeros((len(accessionList), 1))
    thisTrans = (snpMat == 6).astype(int)
    transOverlap[a] = np.sum((np.multiply(snpMat, transFilter) > 0).astype(int))

res ={}
res['enhOverlap'] = enhOverlap
res['transOverlap'] = transOverlap
expTermFile = annotationFolder + 'overlapfor_liverFemal25_.pkl'
expTermFile = dataFolder + 'overlapfor_CD4Pos_%s.pkl' %(accession)
print(expTermFile)
with open(expTermFile, 'wb') as f:
    pickle.dump(res, f)

expTermFile = annotationFolder + 'overlapfor_liverFemal25_.pkl'
print(expTermFile)
with open(expTermFile, 'rb') as f:
    res = pickle.load(f)


kado = np.sort(res['enhOverlap'])
kado = np.sort(res['transOverlap'])
plt.plot(kado)
plt.show()

kado = np.argsort(res['enhOverlap'])
for k in range(210, 234):
    print(accessionList[int(kado[k])])
    print(allMeta[accessionList[kado[k]]]['tissueInfo'])
    print(res['enhOverlap'][kado[k]])

# for one tissue of interest, for it's significant trait,, get the agreement of the other tissues with it regarding the 3 labels: enh, promoter, transcribed. This can be all recorded in a matrix of 3 by tissue count.

# sort the most relevant samples and see what shows. The first choices must be the same tissue.

# Do it for a few diseases. 
            
