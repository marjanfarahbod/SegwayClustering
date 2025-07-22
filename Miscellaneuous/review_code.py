# The code, plots and items for addressing the review comments.
# 1. Get the transcriptomic data for Habib
# 2. Get the plot that Max asked for
# 3. Get the table that you have noted (the training data)
# 4. Get a better Fig01F version
# 5. Print a list of accession vs tissueList

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
from QC_transcriptionComparison_util import Gene, Exon, Annotation, AnnotationClass

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/'

# load the meta file 
inputFileName = 'all235Annot_meta_corrected.pkl'
#inputFileName = 'all235Annot_meta.pkl'
inputFile = dataFolder +  inputFileName
with open(inputFile, 'rb') as f:
    allMeta = pickle.load(f)

# GTF data structure
fileName = dataFolder + '/geneLists.pkl'
with open(fileName, 'rb') as pickleFile:
    geneListsAndIDs = pickle.load(pickleFile)
    
geneList = geneListsAndIDs[0]
geneIDList = geneListsAndIDs[1][0:24771]
del geneListsAndIDs

accessionList = list(allMeta.keys())
accessionList.remove(accessionList[205])
# Segway states:
segwayStates = ['Promoter', 'PromoterFlanking', 'Enhancer', 'EnhancerLow', 'Bivalent', 'CTCF', 'Transcribed', 'K9K36', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent']
print(len(segwayStates))
segwayStateCount = len(segwayStates)

trackList = ['H3K36me3', 'H3K27me3', 'H3K4me3', 'H3K9me3', 'H3K27ac', 'H3K4me1', 'CTCF', 'DNase-seq', 'ATAC-seq', 'POLR2A', 'EP300']

##############################
# 1. Get the transcriptomic data for Habib
##############################

# 1.0: get the gene object dictionary (geneList) to a geneID dictionary
# 1.1: Make a table of expression from each expression data

count = 0
accessionExpressionArray = {}
for accession in accessionList:
    annotation = allMeta[accession]

    annotationFolder = annotation['folder']

    if not(annotation['RNAseqFile'] == 'none') and not(annotation['chromFile'] == 'none'):
        count+=1

        #if count == 10:
        #    break

        if len(annotation['RNAseqFile'][0]) > 1:
            RNAFile = annotation['RNAseqFile'][0]
            print(RNAFile)
        else:
            RNAFile = annotation['RNAseqFile']
            print(RNAFile)

        expAccession = RNAFile[-15:-4]
        expFile = annotationFolder +  'geneExp_dict_' + expAccession + '.pkl'
        
        with open(expFile, 'rb') as pickledFile:
            expression = pickle.load(pickledFile) # load the expression file

        accessionExpressionArray[accession] = expression


file = dataFolder + 'accessionToGeneExpression_88.pkl'
with open(file, 'wb') as f:
    pickle.dump(accessionExpressionArray, f)

##############################
# 2. Get the plot that Max asked for
##############################

# index of samples and labels from the training set
file = dataFolder + 'label_mappings_trainSelect - label_mappings_trainSelect_Jan2024.tsv'

# read lines

import linecache

sample_label_list = []
for li in range(296, 385):

    line = linecache.getline(file, li)
    s = int(line.strip().split('\t')[0].split('_')[1])
    l = int(line.strip().split('\t')[1])
    sample_label_list.append((s, l)) # the original sample index.

# add 3 promoter flanking to this
sample_label_list.append((10,6))
sample_label_list.append((14,3))
sample_label_list.append((15,0))

# Here make the allBarValues from paperFigure02.py. DONE
# Now walk through the sample_label_list and make a list of all the training labels for each of the 11 Segway labels. Then, for each of the 11 labels, fetch that list and get it


# the product of this loop will be a set of sample index and the label index, but the label index is from Segway 11 labels. 
# to get the label: accessionList ID + 38
allBarValues_sampleLabelInd = []
for sl in sample_label_list:
    sample = sl[0]
    print(sample)
    accession = accessionList[sample+38] # sample index in the all sample list
    if accession in rnaAccessionList:
        # find the sample index in the allBarValues
        sampleBarIndex = rnaAccessionList.index(accession)

        # find the label index of the sample from 11 Segway labels (this way, I can fetch the bars from the allBarValues)
        labelInd = sl[1]

        annotation = allMeta[accession]
        annotationFolder = annotation['folder']

        # 1.1 get the mnemonics - we need this fraction coverage
        # >>>>>>>>>>>>>>>>>>>>>>>>>>
        label_term_mapping = {}
        mnemonics_file = annotationFolder + 'mnemonics_v04.txt'
        with open(mnemonics_file, 'r') as mnemonics:
            for line in mnemonics:
                print(line)
                label = line.strip().split()[0]
                term = line.strip().split()[1]
                label_term_mapping[label] = term


        stateInd = segwayStates.index(label_term_mapping[str(labelInd)])
        
        allBarValues_sampleLabelInd.append((sampleBarIndex, stateInd))
        

sumBars = np.zeros((3, 11, 160)) # define the size
stateLabelCount = np.zeros((11)) # for the mean divide
for bar in allBarValues_sampleLabelInd:
    # read the bar, add the values to the meanBars
    sampleInd = bar[0]
    barsInd = bar[1]

    sumBars[0, barsInd, :] += allBarValues[sampleInd, 0, barsInd, :]
    sumBars[1, barsInd, :] += allBarValues[sampleInd, 1, barsInd, :]
    sumBars[2, barsInd, :] += allBarValues[sampleInd, 2, barsInd, :]

    stateLabelCount[barsInd]+= 1

meanBars = np.zeros((3, 11, 160))
for bar in range(11):
    if stateLabelCount[bar] > 0:
        meanBars[0, bar, :] = sumBars[0, bar, :]/ stateLabelCount[bar]
        meanBars[1, bar, :] = sumBars[1, bar, :]/ stateLabelCount[bar]
        meanBars[2, bar, :] = sumBars[2, bar, :]/ stateLabelCount[bar]

# do the plot

fig, axs = plt.subplots(segwayStateCount, 3, figsize=(12,8))

xticks = [30, 130]
xticksLabels = ['TSS', 'TTS']

indexList = np.array(list(range(160)))

#thisMat = expMats['clusterMats'][0]

#logMat = np.log10(thisMat, out=np.zeros_like(thisMat), where=(thisMat!=0))
for m in range(3):
    thisMeanBars = meanBars[m, :, :]
    for s in range(segwayStateCount):
        positiveInds = indexList[thisMeanBars[s,:] >= 0]
        negativeInds = indexList[thisMeanBars[s,:] < 0]
        posBar = np.copy(thisMeanBars[s, :])
        posBar[negativeInds]=0
        axs[s,m].bar(range(160), posBar, color=[249/256, 80/256, 97/256,1], width=1)
        negBar = np.copy(thisMeanBars[s, :])
        negBar[positiveInds]=0
        axs[s,m].bar(range(160), negBar, color=[50/256,164/256,249/256,1], width=1)
        axs[s,m].set_ylim((-2,2))
        ylabel = segwayStates[s]
        #axs[i,0].text(60, .55, ylabel, fontsize=8)
        axs[s,m].set_xticks(xticks)
        axs[s,m].set_yticks([-1, 1])
        
        if m == 0:
            axs[s,m].set_yticklabels([-1, 1], fontsize=8)
            axs[s,m].set_ylabel(ylabel, rotation=0, fontsize=10, labelpad=3, ha='right', va='center')

    axs[segwayStateCount-1,m].set_xticklabels(xticksLabels)
    axs[segwayStateCount-1,m].set_xlabel('Position relative to gene')


axs[0,0].set_title('Enrichment of labels at genes with \nzero expression (log10(observed/expected))', fontsize = 9)
axs[0,1].set_title('Enrichment of labels at genes with \nbottom 30% expression (log10(observed/expected))', fontsize = 9)
axs[0,2].set_title('Enrichment of labels at genes with \ntop 70% expression (log10(observed/expected))', fontsize = 9)

#figFile = plotFolder + figureFolder + 'transcriptomicPlot_median_QartileLines.pdf'
figureFolder = 'figure02/'
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/'
figFile = plotFolder + figureFolder + 'transcriptomicPlot_mean_trainingDataOnly.pdf'
plt.savefig(figFile)
plt.close('all')

    

# looks good, let's finish and wrap it up. 

##############################
# 3. The table
##############################



#########################################
# 4. Get a better Fig01F version
#########################################

file = dataFolder + 'accessionToGeneExpression_88.pkl'
with open(file, 'rb') as f:
    expData = pickle.load(f)

accessions = list(expData.keys())

for accession in accessions: # all have the same length
    print(len(expData[accession]))

# getting the expression levels that I care about (ones in the geneIDList)

expMat = np.zeros((len(geneIDList), len(accessions)))
for a, accession in enumerate(accessions):
    myExpData = expData[accession]
    expVals = [myExpData[gene] for gene in geneIDList]
    expMat[:, a] = np.array(expVals)

# count of genes is too high, let's remove the genes which are expressed in <3 samples.
# how to get expression filter:

myLowQs = np.zeros((88,1))
for i in range(expMat.shape[1]):
    myExps = expMat[:, i]
    myExps = myExps[myExps > 0]
    #print(np.quantile(myExps, .3))
    myLowQs[i] = np.quantile(myExps,.3)

print(np.mean(myLowQs)) # roughly .3

expBool = expMat > .3
expSum = np.sum(expBool, 1)
fil1 = expSum > 4 # filter1: genes that are expressed in 2 or more samples
print(np.sum(fil1))

# fil2, gene length
genesLength = np.asarray([geneList[g].end - geneList[g].start for g in geneIDList])
print(len(genesLength))

plt.plot(np.sort(genesLength))
plt.show()

print(np.sum(np.asarray(genesLength) < 1000))
fil2 = genesLength > 1000
print(np.sum(fil2))

print(np.median(genesLength))
print(np.quantile(genesLength, .75))

thr = np.median(genesLength) + 2*(np.quantile(genesLength, .75) - np.median(genesLength))
fil3 = genesLength < thr
print(np.sum(fil3))

tf = np.logical_and(fil1, fil2)
finalFil = np.logical_and(tf, fil3)
print(np.sum(finalFil ==1)) # count of genes which passed the filter

filExpMat = expMat[finalFil, :] # *** got the array
geneIDListArray = np.asarray(geneIDList) 
filGeneIDs = geneIDListArray[finalFil] # *** got the IDs

# just do a hierarchical clustering
import seaborn as sns

mydf = pd.DataFrame(filExpMat[4000:8000,])
mydf = pd.DataFrame(filExpMat[0:4000,])
mydf = pd.DataFrame(filExpMat[8000:12000,])
sib = sns.clustermap(mydf, figsize=(6, 150), method = 'average', dendrogram_ratio = [.25,.01])
sns.set(font_scale = .4)

figFile = '%sclustering_genes_filtered14806_0_4000_lowExpFil.pdf' %(plotFolder)
figFile = '%sclustering_genes_filtered14806_4000_8000_lowExpFil.pdf' %(plotFolder)
print(figFile)
plt.savefig(figFile)
plt.close('all')

myInds_tissues = sib.dendrogram_col.reordered_ind
myInds_genes = sib.dendrogram_row.reordered_ind

for ind in myInds:
    print(allMeta[accessions[ind]]['tissueInfo'])
# note: they are absolutely clustered - now let's get the genes.

print(myInds_genes.index(1434)) # plot that region with gene and tissue names, this was gene 195
print(myInds_genes.index(2434)) # plot that region with gene and tissue names, this was gene 2434, index was 2214
print(myInds_genes.index(1784)) # index was 42
print(myInds_genes.index(736)) # index was 180

sg = 178
eg = 193
# from here, TNNC1 seems like a good candidate, let's get the ID and index

sg = 206
eg = 250

sg = 2214
eg = 2223

sg = 47
eg = 49

myTissueList = [allMeta[accessions[x]]['tissueInfo'][0] for x in myInds_tissues]
myTissueAccessions = [accessions[x] for x in myInds_tissues]
#myGeneIDs = [filGeneIDs[x+4000] for x in myInds_genes[sg:eg]] # for the second matrix
myGeneIDs = [filGeneIDs[x] for x in myInds_genes[sg:eg]] # for the second matrix
myGeneNames = [geneList[x].name for x in myGeneIDs]

subDF = mydf.loc[myInds_genes[sg:eg], myInds_tissues]
subDF = subDF.set_axis(myGeneNames, axis=0)
subDF = subDF.set_axis(myTissueList, axis=1)

sns.heatmap(subDF)
plt.show()

TS_geneListID = []
for g in myGeneIDs:
    print(geneList[g].ENS_ID)
    TS_geneListID.append(geneList[g].ENS_ID)

# name = TNNC1, ENS_ID = ENSG00000114854.7, gtype = protein_coding, chrom = chr3, start = 52451102, end = 52454070, strand = -, exon_count = 0
gid = 'ENSG00000114854.7'
gidc = geneList[gid].end - geneList[gid].start

#name = ACTN2, ENS_ID = ENSG00000077522.12, gtype = protein_coding, chrom = chr1, start = 236686454, end = 236764631, strand = +, exon_count = 0
gid = 'ENSG00000077522.12'
gidc = geneList[gid].end - geneList[gid].start

# name = CKMT2, ENS_ID = ENSG00000131730.15, gtype = protein_coding, chrom = chr5, start = 81233285, end = 81266397, strand = +, exon_count = 0
gid = 'ENSG00000131730.15'
gidc = geneList[gid].end - geneList[gid].start

# name = ACTA1, ENS_ID = ENSG00000143632.14, gtype = protein_coding, chrom = chr1, start = 229431245, end = 229434098, strand = -, exon_count = 0
gid = 'ENSG00000143632.14'
gidc = geneList[gid].end - geneList[gid].start

# name = TNNC1, ENS_ID = ENSG00000114854.7, gtype = protein_coding, chrom = chr3, start = 52451102, end = 52454070, strand = -, exon_count = 0
gid = 'ENSG00000114854.7'
gidc = geneList[gid].end - geneList[gid].start


import linecache
# define the matrix for collecting
chrList = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
           'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']


theList = ['ADPRHL1', 'MYL2', 'NRAP', 'TNNT2', 'CASQ2', 'ACTN2', 'ACTA1']
TS_geneListID = []
for gene in geneList:
    if geneList[gene].name in theList:
        TS_geneListID.append(gene)


for gene in TS_geneListID:
    print(geneList[gene])
    if geneList[gene].strand == '+':
        scor = geneList[gene].start - 50000
    else:
        scor = geneList[gene].end - 50000

    print(scor)
    geneSymbol = geneList[gene].name
    regionChr = (geneList[gene].chrom) # 'chr15'
    regionChrInd = chrList.index(regionChr)
    regionPlot(geneSymbol, scor, regionChr, regionChrInd)

def regionPlot(geneSymbol, scor, regionChr, regionChrInd):

    # scor: starting coordinate of the region
    import matplotlib.colors as colors
    
    W = 1000 #300
    annotMat = np.zeros((234, W))

    #regionChr = 'chr2' # 'chr15'
    #regionChr = (geneList[gid].chrom) # 'chr15'
    #regionChrInd = chrList.index(regionChr)
    #for accession in myTissueAccessions[0:40]:
    # for each of the accessions
    for accession in accessionList:
        
        print(allMeta[accession]['tissueInfo'][0])

        if accession == 'ENCSR424JDX':
            print(accessionList.index(accession))
            continue
    
        a = accessionList.index(accession)

        print(a, accession)
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

        print(sind)
        print(scor)

        walker = 0
        steps = int((eind - scor)/100)
        color = label_term_mapping[line.split()[3].split('_')[0]]
        colorInd = segwayStates.index(color)
        annotMat[a, walker:walker+steps-1] = colorInd
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
            annotMat[a, walker:walker+steps] = colorInd
            walker += steps
        
        linecache.clearcache()
        os.system('gzip %s' %(annFile))


    print(annotMat.shape)
    outputFileName = 'annotationMap234Cells_%s_allSamples_longRegion_100000.pkl' %(geneSymbol)
    outputFile = dataFolder + outputFileName
    with open(outputFile, 'wb') as f:
        pickle.dump(annotMat, f)


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


    cmap = colors.ListedColormap(colorList)
    boundaries = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    norm = colors.BoundaryNorm(boundaries, cmap.N, clip=True)

    df = pd.DataFrame(annotMat)

    sib = sns.clustermap(df, method = 'ward', col_cluster=False)
    rowInds = sib.dendrogram_row.reordered_ind
    plt.close('all')

    sns.heatmap(annotMat[rowInds,:], cmap=cmap, norm=norm)
    figFile = plotFolder + 'annotMat_%s_allSample_longRegion_100000.pdf' %(geneSymbol)
    plt.savefig(figFile)
    plt.close('all')


geneSymbol = 'NRAP'
inputFileName = 'annotationMap234Cells_%s_allSamples_longRegion_100000.pkl' %(geneSymbol)
inputFile = dataFolder + inputFileName
with open(inputFile, 'rb') as f:
    annotMat = pickle.load(f)

    
df = pd.DataFrame(annotMat)
sib = sns.clustermap(df, method = 'ward', col_cluster=False)
rowInds = sib.dendrogram_row.reordered_ind
plt.close('all')

for i,sind in enumerate(rowInds):

    accession = accessionList[sind]
        
    print(i, sind,allMeta[accession]['tissueInfo'][0])

    if accession == 'ENCSR424JDX':
        print(accessionList.index(accession))
        continue
    
    a = accessionList.index(accession)

    print(a, accession)
    annotation = allMeta[accession]


sns.heatmap(annotMat[rowInds,:], cmap=cmap, norm=norm)
figFile = plotFolder + 'annotMat_%s_allSample_longRegion_100000_test.pdf' %(geneSymbol)
plt.savefig(figFile)
plt.close('all')

# ########################################
# 5. Print a list of accession vs tissueList
# ########################################

        label_term_mapping = {}
        mnemonics_file = annotationFolder + 'mnemonics_v04.txt'
        with open(mnemonics_file, 'r') as mnemonics:
            for line in mnemonics:
                print(line)
                label = line.strip().split()[0]
                term = line.strip().split()[1]
                label_term_mapping[label] = term

tissueFile = dataFolder + 'accession_tissue_info.txt'
with open(tissueFile, 'w') as outputFile:
    for i in range(len(accessionList))


