# 0. Initials
# 1. Enhancer file exploring
# 2. Matching tissues
# 3. Get the feature matrix
# 4. Prediction
# 5. Print the matching info for the supplement table
# 6. Get the coverage prediction like transcriptomic
# 7. Get the plot Max mentioned


# so here we are.
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
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/'
figureFolder = 'figure02/'


inputFileName = 'all235Annot_meta_corrected.pkl'
#inputFileName = 'all235Annot_meta.pkl'
inputFile = dataFolder +  inputFileName
with open(inputFile, 'rb') as f:
    allMeta = pickle.load(f)

fanData = dataFolder + 'Fantom5Enh_data/human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix2.txt'
fandf = pd.read_csv(fanData, sep='\t')

fanData_sorted = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/Fantom5Enh_data/human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix2_sorted.txt'

fanSampleInfo = dataFolder + 'Fantom5Enh_data/Human.sample_name2library_id.txt'
metadf = pd.read_csv(fanSampleInfo, sep='\t')


########################################
# 1. Enhancer file exploring
########################################
# count of regions: 65423, 1828

# get list of all rows in the ID column (the header is actually refereing to the list of sample IDs as column header and is not a header for the coordinate column, but here it works like this)
book = fandf['Id'].values.tolist()
len(book)

lengths = np.zeros(len(book))

import re
for i,coor in enumerate(book):
    fields = re.split(':|-', coor)
    lengths[i] = int(fields[2]) - int(fields[1])

sl = np.sort(lengths)

plt.plot(sl)
plt.show()

a = sum(lengths > 1000)/len(lengths)
b = sum(lengths < 100)/len(lengths)

print(sum(lengths > 500)/len(lengths))
print(sum(lengths < 100)/len(lengths))

print(1 - (a+b))

# 84% of the enhancers identified are between 100-500bps
# 90% are between 100-100bps

# how many tissues do enhancer appear in? sum of the rows
# how many enhancers does a tissue have? sum of the columns

colList = fandf.columns.values.tolist()
print(len(colList))

fanArr = fandf[colList[2:]].to_numpy()
print(fanArr.shape)

# how many enhancers does a tissue have? sum of the columns
########################################
fanArr_present = fanArr.copy()
fanArr_present[fanArr_present > 0] = 1
enhCount = np.sum(fanArr_present, axis=0)

plt.plot(np.sort(enhCount))
plt.show()

enhCount = np.sort(enhCount)

a = sum(enhCount < 2000)/len(enhCount)
b = sum(enhCount > 15000)/len(enhCount)

print(1 - (a+b))
# 81% have between 2k-15k enhancers of theses

# how many tissues do enhancer appear in? sum of the rows
########################################
tissueCount = np.sum(fanArr_present, axis=1)
plt.plot(np.sort(tissueCount))
plt.show()

tissueCount = np.sort(tissueCount)

a = sum(tissueCount < 1)/len(tissueCount)
b = sum(tissueCount > 200)/len(tissueCount)

print(1 - (a+b))
# 91% of enhancers have appeared in < 350 tissues
# 74% appear in <200 tissues

########################################
# 2. Matching tissues
########################################

# get the list of tissues from the samples meta info file. Get the ID of my tissues, match them.
fanSampleInfo = dataFolder + 'Fantom5Enh_data/Human.sample_name2library_id.txt'
tissueNameID_map = {}
with open(fanSampleInfo, 'r') as f:
    for line in f:
        tissueNameID_map[line.split('\t')[0]] = line.split('\t')[1]

fanTissueKeyList = list(tissueNameID_map.keys())
fanTissueKeyList_lower = []
for tissue in fanTissueKeyList:
    fanTissueKeyList_lower.append(tissue.lower())
# do the thing

myTissueList = []
for accession in accessionList:
    annotation = allMeta[accession]
    myTissueList.append(annotation['tissueInfo'][0])

myTissueList_lower = []
for tissue in myTissueList:
    myTissueList_lower.append(tissue.lower())

fullInds = []
maxFanListI = []
maxFanListT = []
for t, tissue in enumerate(myTissueList_lower):

    if not(t in myFinalSet):
        continue
    print('*****************')
    print(tissue)
    print(t)
    tissueIndList = []
    tissueParts = tissue.split(' ')
    maxPart = 0
    maxFanInd = 0
    for fanTissue in fanTissueKeyList_lower:
        partCount = 0
        fanInd = fanTissueKeyList_lower.index(fanTissue)
        for part in tissueParts:
            if part != 'cell' and len(part)>1:
                if part in fanTissue:
                    partCount +=1
                    print('--', fanInd ,fanTissue)
                    book = (fanInd, fanTissue)
                    tissueIndList.append(book)
                    if part == 'heart':
                        partCount +=15
                    if part == 'liver':
                        partCount +=15
                    if part == 'lung':
                        partCount +=15
                
        if partCount > maxPart:
            maxPart = partCount
            maxFanInd = fanInd
            maxFanT = fanTissue
    maxFanListI.append(maxFanInd)
    maxFanListT.append(maxFanT)
    print('>>>>>>>', maxFanT)

                
    fullInds.append(tissueIndList)

# get the tissue file.

# difficult indexes:
lookInds = [2, 11, 12, 14, 18, 21, 23, 24, 27, 33, 35, 37, 38, 40, 42, 43, 44, 45, 47, 48, 49, 52, 54, 55, 56, 58, 60, 61, 62, 68, 71, 74,  75, 76, 79, 81, 83, 85, 88, 89, 90, 91, 95, 97, 98, 100, 101, 105, 109, 110, 111, 112, 116, 119, 127, 129, 130, 131, 132, 134, 137, 144, 145, 146, 150,151, 152,154, 156,157,160,161,165,166,167,169,170,175,176, 177, 178, 179, 180,181, 182,183,184,185,189 ,190,193 ,196,197,198,199,202 ,208,209,210,211,212,213,215,221 ,227,229, 231,232]

kado = range(234)
mySet = set(kado)
myFinalSet = mySet - set(lookInds)

# to fix:
book = {16: 1712,19: 356, 26: 307, 30:1649, 31:1813,32:1654,36:'r', 39:1712, 41:'r',
        46:440, 51:1681, 63:1681, 66:1660, 72:1712, 77:'r', 78:1712, 82:1712, 93:1712,
        92:'r',113:1712,120:1651,124:1660, 125:'r',133:1712,136:1712,138:323,139:1712,
        143: 'r', 153:1712,155:564, 171: 'r',187: 1712,188: 1712,194:1681, 201:360,205: 'r',
        206: 'r', 219: 1712,220: 1681}

decided = list(book.keys())

#fullInds = []
#maxFanListI = []
#maxFanListT = []
allMappingInfo = {}
for t, tissue in enumerate(myTissueList_lower):

    if not(t in myFinalSet):
        continue
    
    accession = accessionList[t]
    if t in decided and book[t]!='r':
        info = (t, tissue, book[t], fanTissueKeyList[book[t]])
        print('>>>>>>>>>')
        print(info)
        print(accession)
        print(allMeta[accession]['tissueInfo'])
    else: 
        tissueParts = tissue.split(' ')
        maxPart = 0
        maxFanInd = 0
        for fanTissue in fanTissueKeyList_lower:
            partCount = 0
            fanInd = fanTissueKeyList_lower.index(fanTissue)
            for part in tissueParts:
                if part != 'cell' and len(part)>1:
                    if part in fanTissue:
                        partCount +=1
                        #print('--', fanInd ,fanTissue)
                        if part == 'heart':
                            partCount +=15
                        if part == 'liver':
                            partCount +=15
                        if part == 'lung':
                            partCount +=15
                
            if partCount > maxPart:
                maxPart = partCount
                maxFanInd = fanInd
                maxFanT = fanTissue
        maxFanListI.append(maxFanInd)
        maxFanListT.append(maxFanT)
        print('>>>>>>>', maxFanT)
        info = (t, tissue, maxFanInd, maxFanT)
        print(info)
        print(accession)
        print(allMeta[accession]['tissueInfo'])

    allMappingInfo[accession] = info
    #fullInds.append(tissueIndList)

outputFileName = 'fantom5Mapping.pkl'
outputFile = dataFolder + outputFileName
print(outputFile)
with open(outputFile, 'wb') as f:
    pickle.dump(allMappingInfo, f)

inputFileName = 'fantom5Mapping.pkl'
inputFile = dataFolder + inputFileName
print(inputFile)
with open(inputFile, 'rb') as f:
    allMappingInfo = pickle.load(f)

########################################
# 3. Get the feature matrix
########################################

inputFileName = 'fantom5Mapping.pkl'
inputFile = dataFolder + inputFileName
print(inputFile)
with open(inputFile, 'rb') as f:
    allMappingInfo = pickle.load(f)

enhAccessionList = list(allMappingInfo.keys())
    
# what is the expression distribution?
book = fandf[colList[10]]
book = np.array(book)
print(sum(book>.1))
plt.hist(book)
plt.show()

book = np.arange(0, 60000, 100)
for i in book:
    print(fandf.loc[i].at["Id"])
    
import subprocess as sp
import linecache
chrList = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
           'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
           'chr18', 'chr19', 'chr20', 'chr21', 'chr22']

annLineCount = int(sp.getoutput('wc -l %s' %(fanData)).split()[0])

line = linecache.getline(fanData, 3)

pervChr = 'chr1'
for i in range(4, annLineCount):
    ml = linecache.getline(fanData, i)[0:10]
    curChr = ml.split(':')[0]
    if curChr != pervChr:
        print(i, pervChr, curChr)
        pervChr = curChr
'''>>> 6003 chr1 chr10
9300 chr10 chr11
12263 chr11 chr12
15517 chr12 chr13
17219 chr13 chr14
19315 chr14 chr15
21497 chr15 chr16
23597 chr16 chr17
26519 chr17 chr18
28075 chr18 chr19
30102 chr19 chr2
35605 chr2 chr20
37530 chr20 chr21
38589 chr21 chr22
39740 chr22 chr3
43846 chr3 chr4
46814 chr4 chr5
50430 chr5 chr6
54714 chr6 chr7
57988 chr7 chr8
61290 chr8 chr9
63993 chr9 chrX
65367 chrX chrY'''

chrLines = [[2, 6002],
            [30102, 35604],
            [39740, 43845],
            [43846, 46813],
            [46814, 50429],
            [50430, 54713],
            [54714, 57987],
            [57988, 61289],
            [61290, 63992],
            [6003, 9299],
            [9300, 12262],
            [12263, 15516],
            [15517, 17218],
            [17219, 19314],
            [19315, 21496],
            [21497, 23596],
            [23597, 26518],
            [26519, 28074],
            [28075, 30101],
            [35605, 37529],
            [37530, 38588],
            [38589, 39739]]
            
for arr in chrLines:
    print(linecache.getline(fanData, arr[0])[0:25])
    print(linecache.getline(fanData, arr[1])[0:25])

# there is an offset of 2 for the dataframe        
print(fandf.loc[6000].at["Id"])
print(fandf.loc[6003].at["Id"])

fanData_sorted = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/Fantom5Enh_data/human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix2_sorted.txt'

# 3.1 sorting the fanData
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# read the lines from the fanData according to the indices, write them in the fanData_sorted.
with open(fanData_sorted, 'w') as output:
    ml = linecache.getline(fanData, 0)        
    output.write(ml)
    ml = linecache.getline(fanData, 1)        
    output.write(ml)

    for arr in chrLines:
        s = arr[0]
        e = arr[1]
        for i in range(s, e+1):
            ml = linecache.getline(fanData, i)
            output.write(ml)

# test
pervChr = 'chr1'
for i in range(2, annLineCount):
    ml = linecache.getline(fanData_sorted, i)[0:10]
    curChr = ml.split(':')[0]
    if curChr != pervChr:
        print(i, pervChr, curChr)
        pervChr = curChr

# 3.2 read the enhancers into a tupple, dictionary based on the line index
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
annLineCount = int(sp.getoutput('wc -l %s' %(fanData_sorted)).split()[0])
enh_dict = {}
with open(fanData_sorted, 'r') as f:
    for i in range(2,annLineCount-2):
        ml = linecache.getline(fanData_sorted, i)[0:30].split('\t')[0]
        sml = ml.split(':')
        chr = sml[0]
        s = int(sml[1].split('-')[0])
        e = int(sml[1].split('-')[1])
        enh_dict[i] = (chr, s, e)

# 3.3 make another dict that will give me the line number based on the tuple
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
'''
in the out of gcsc, the start position of the old coordinates have the offset of exactly 1, either higher or lower.
The end position remains the same, so I am only going to use that as the key
'''
enh_dict_lineInd = {}
for i in range(2, annLineCount-2):
    myKey = enh_dict[i]
    myNewKey = (myKey[0], myKey[2])
    enh_dict_lineInd[myNewKey] = i

# 3.4 convert to the hg20 genome assembly
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

'''
I will get the enhancer coverages based on the hg20 sorted enhancer coordinates.
To fetch the expression values later, I will find the line index by mapping of hg20 to hg19 to line index 
'''

# 3.4.1 just writing the coordinates in the file to upload to gcsc to get new coordinates
# ..............................
annLineCount = int(sp.getoutput('wc -l %s' %(fanData_sorted)).split()[0])
outputFile = dataFolder + '/Fantom5Enh_data/bedCoordinates.txt'
with open(outputFile, 'w') as of:
    for i in range(0,annLineCount-2):
        ml = linecache.getline(fanData_sorted, i+2)[0:30].split('\t')[0]
        sml = ml.split(':')
        chr = sml[0]
        s = int(sml[1].split('-')[0])
        e = int(sml[1].split('-')[1])
        of.write('%s\t%d\t%d\n' %(chr, s, e))

# 3.4.2 Getting the line index for the new coordinates 
# ..............................
# based on the enh chromosome and start, give the coordinates.
updateCoordinates = dataFolder + '/Fantom5Enh_data/hglft_genome_39b08_f8a700.bed'
lineCount = int(sp.getoutput('wc -l %s' %(updateCoordinates)).split()[0])
assemblyMap = {} # map from the old coordinates to the new coordinates
assemblyMap2 = {} # map from new coordinates to old coordinates
fantomLineIndMap = {} # map from the new coordinates to the line index of the original file
for i in range(0, lineCount-4):
    ml = linecache.getline(updateCoordinates, i+1)
    sml = ml.split('\t')
    oldChr = sml[3].split(':')[0]
    newChr = sml[0]
    oldEnd = int(sml[3].split('-')[1])
    newEnd = int(sml[2])
    #oldStart = int(sml[3].split('-')[0].split(':')[1])
    newStart = int(sml[1])
    oldCoors = (oldChr, oldEnd)
    newCoors = (newChr, newStart, newEnd)
    if newChr in chrList:
        assemblyMap[oldCoors] = newCoors
        assemblyMap2[newCoors] = oldCoors
        fantomLineIndMap[newCoors] = enh_dict_lineInd[oldCoors]


# checking to see if all match.
# this will show me if the new coordiates get me to the right line index in fanData_sorted
# it is important since my expressions are in fanData_sorted, and should be obtained from there
# the order of the 'newCoordinates' doesn't matter anymore, since I know where to find their expression
i = 30000
ml = linecache.getline(updateCoordinates, i)
print(ml)
sml = ml.split('\t')
newChr = sml[0]
newEnd = int(sml[2])
newStart = int(sml[1])
newCoors = (newChr, newStart, newEnd)
print(newCoors)
print(fantomLineIndMap[newCoors])
l = fantomLineIndMap[newCoors]
originalLine = linecache.getline(fanData_sorted, l)
print(originalLine[0:100])

# 3.4.3 sorting the new coordinates before getting into the annotation file
# ..............................
sortedCoorFile = dataFolder + '/Fantom5Enh_data/sorted_hglft_genome_39b08_f8a700.bed'
os.system('sort -V -k1,1 -k2,2 %s > %s' %(updateCoordinates, sortedCoorFile))

# checking the lines and if sort worked - that is fine - looks fine
preS = 0
with open(sortedCoorFile, 'r') as f:
    for line in f:
        s = int(line.split('\t')[1])
        if s < preS:
            print(line)
        preS = s

with open(sortedCoorFile, 'r') as f:
    for line in f:
        chr = line.split('\t')[0]
        if not(chr in chrList):
            print(line)


# 3.4.5 reading the new sorted coordinates 
# ..............................
newEnhCoors = {}
l = 0
with open(sortedCoorFile, 'r') as f:
    for line in (f):
        myChr = line.split('\t')[0]
        if myChr in chrList:
            myStart = int(line.split('\t')[1])
            myEnd = int(line.split('\t')[2])
            newEnhCoors[l] = (myChr, myStart, myEnd)
            l +=1


# getting the line index for each of the coors:
enhInds = np.zeros(len(newEnhCoors))
for i in range(len(newEnhCoors)):
    newCoors = newEnhCoors[i]
    enhInds[i] = fantomLineIndMap[newCoors]-2

###### that function for Segway

chrIndex = {'chr1':1, 'chr2':2, 'chr3':3, 'chr4':4, 'chr5':5, 'chr6':6, 'chr7':7, 'chr8':8, 'chr9':9,
            'chr10':10, 'chr11':11, 'chr12':12, 'chr13':13, 'chr14':14, 'chr15':15, 'chr16':16,
            'chr17':17, 'chr18':18, 'chr19':19, 'chr20':20, 'chr21':21, 'chr22':22}

chrList = list(chrIndex.keys())

# get the coverage of the region for each of the labels:
for accession in enhAccessionList:
    print(accession)
    annotation = allMeta[accession]
    annFolder = annotation['folder']

    # open the .bed file
    annFile = annotation['bedFile']

    # prepare ann
    if annFile.endswith('.gz'):
        os.system('gunzip %s' %(annFile))
        annFile = annFile[0:-3]
    else:
        os.system('gunzip %s.gz' %(annFile))

    # open the mnemnoics file to get the count    
    mnemFile = annFolder + 'mnemonics_v04.txt'
    label_term_mapping = {}
    with open(mnemFile, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term

    clusterCount = len(label_term_mapping)
    sortedClusterList = list(range(0,clusterCount))

    cgi = 0 # walks on the geneIDList
    ann_start = 0 # the start of the current annotation
    ann_end = 0 # the end of the current annotation
    ann_chr = 'chr'
    ann_line_count = 0 # this is just to check the progress through the annotation file
        
    previous_class = ''

    #labelExpMat = np.zeros((len(enh_dict), clusterCount))
    labelExpMat = np.zeros((len(enh_dict), segwayStateCount))
    
    annLineInd = 1 # this keeps the annotation line index, the first line is empty, thus the index starts at 1

    previous_gene_chr = 'chr1'
    previous_extension_end = 0 # doesn't matter since chr condition is never true until the first gene is processed
    previous_ann_chr = 'chr1'

    while cgi < len(newEnhCoors): # >>>>>>>>>> MAIN 24124

        'we are in the gene territory, we are walking on the genes'
        # enh_dict
        gene_chr = newEnhCoors[cgi][0]
            
        gene_start = newEnhCoors[cgi][1]
        gene_end = newEnhCoors[cgi][2]

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

        #if (gene_chr != previous_gene_chr): # in case of chromosome change because gene moved to the next chromosome
        if(chrIndex[gene_chr] > chrIndex[ann_chr]):
            #print('gene chr not equal to previous gene chr')
            #print('gene_chr %s' %(gene_chr))
            #print('p_gene_chr %s' %(previous_gene_chr))
            #print('ann_chr %s' %(ann_chr))
                
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

            #print(ann_chr)
            #print(previous_ann_chr)

                
        #if (ann_chr != gene_chr): # if annotation moved to the next chromosome, but gene has not yet moved to the next chromosome
        if (chrIndex[gene_chr] < chrIndex[ann_chr]):
            annLineInd = annLineInd - 2
            line = linecache.getline(annFile, annLineInd)
            annLineInd +=1
            ann_line_count += 1
            fields = line.strip().split()

            ann_chr = fields[0]
            ann_start = int(fields[1])
            ann_end = int(fields[2])
            #print('marker 01')


        ''' 
        if this gene starts somewhere in the preivous genomic region, 
        I will go back in the annotation file to the beginning of annotation for the previous gene
        changed this th the next while loop - basically we will go back until the annotations start before the gene territory

        '''

        while (ann_start > gene_start) and (gene_chr == ann_chr): # in case of overlapping genes
                
            annLineInd = max(annLineInd - 5,1)
            line = linecache.getline(annFile, annLineInd)
            annLineInd +=1
            ann_line_count += 1
            fields = line.strip().split()

            ann_chr = fields[0]
            ann_start = int(fields[1])
            ann_end = int(fields[2])

        while ((ann_start < gene_end) and (gene_chr == ann_chr)) and gene_coverage < gene_length: #and geneMatWalkIndex < 160: # while we still have annotation before the gene

            # in the next sections we are processing the annotation
            if not((ann_end < gene_start) or ann_start > gene_end):

                ''' We are in the genonimc region (with extension)'''

                ''' Taking off for filling the matrices... '''
                        
                adjusted_ann_start = max(0, ann_start - gene_start)
                adjusted_ann_end = min(ann_end - gene_start, gene_end - gene_start)
                adjusted_ann_length = adjusted_ann_end - adjusted_ann_start

                ann_cluster = fields[3].split('_')[0]
                clusterInd = int(ann_cluster)
                stateInd = segwayStates.index(label_term_mapping[str(clusterInd)])

                if gene_coverage < gene_length and adjusted_ann_length > 0:
                    new_coverage = min(adjusted_ann_length, gene_length - gene_coverage)
                    labelExpMat[cgi, stateInd] += new_coverage
                    adjusted_ann_length = adjusted_ann_length - new_coverage
                    gene_coverage += new_coverage
                    #coverPromoter
 
                
            if gene_coverage < gene_length:
                line = linecache.getline(annFile, annLineInd)
                annLineInd +=1

                ann_line_count += 1
                fields = line.strip().split()

                ann_chr = fields[0]
                ann_start = int(fields[1])
                ann_end = int(fields[2])

        cgi += 1 # next gene
        previous_gene_end = gene_end
        previous_gene_chr = gene_chr
            

    linecache.clearcache()
    os.system('gzip %s' %(annFile))

    outputFile = annFolder  + 'enhancer_labelCoverage_segway_states_hg20.pkl' 
    with open(outputFile, 'wb') as f:
        pickle.dump(labelExpMat, f)

    print(outputFile)
    print(enhAccessionList.index(accession))
    
##### that function for chromhmm

book =  'EnhA1 EnhA2 EnhBiv EnhG1 EnhG2 EnhWk Het Quies ReprPC ReprPCWk TssA TssBiv TssFlnk TssFlnkD TssFlnkU Tx TxWk ZNF/Rpts'
chromLabels = book.split()

for accession in enhAccessionList[26:]:

    annotation = allMeta[accession]

    if not(annotation['chromFile'] == 'none'):
    
        annFolder = annotation['folder']

        chmmFileName = annotation['chromFile'].split('/')[-1]

        if (('38batch' in annotation['folder']) or ('May11' in annotation['folder'])):
            annFile = annFolder + 'sorted_' + chmmFileName
        else:
            annFile = annFolder + chmmFileName

        # prepare ann
        if annFile.endswith('.gz'):
            os.system('gunzip %s' %(annFile))
            annFile = annFile[0:-3]
        else:
            os.system('gunzip %s.gz' %(annFile))

        clusterCount = len(chromLabels)
        annLineCount = int(sp.getoutput('wc -l %s' %(annFile)).split()[0])

        sortedClusterList = chromLabels

        #print(sortedClusterList)
            
        cgi = 0 # walks on the geneIDList
        ann_start = 0 # the start of the current annotation
        ann_end = 0 # the end of the current annotation
        ann_chr = 'chr'
        ann_line_count = 0 # this is just to check the progress through the annotation file
    
        previous_class = ''

        labelExpMat = np.zeros((len(enh_dict), clusterCount))
        
        annLineInd = 1 # this keeps the annotation line index, the first line is empty, thus the index starts at 1

        previous_gene_chr = 'chr1'
        previous_extension_end = 0 # doesn't matter since chr condition is never true until the first gene is processed
        previous_ann_chr = 'chr1'

        while cgi < len(newEnhCoors): # >>>>>>>>>> MAIN

            'we are in the gene territory, we are walking on the genes'
            #enh_dict
            gene_chr = newEnhCoors[cgi][0]
                
            gene_start = newEnhCoors[cgi][1]
            gene_end = newEnhCoors[cgi][2]

            gene_length = gene_end - gene_start
            gene_coverage = 0

            line = linecache.getline(annFile, annLineInd)
            annLineInd +=1
            ann_line_count += 1
            fields = line.strip().split()
        
            # reading the next annotation
            previous_ann_chr = ann_chr
            ann_chr = fields[0]
            ann_start = int(fields[1])
            ann_end = int(fields[2])

            #if (gene_chr != previous_gene_chr): # in case of chromosome change because gene moved to the next chromosome
            if (chrIndex[gene_chr] > chrIndex[ann_chr]):
            #print('gene chr not equal to previous gene chr')
            #print('gene_chr %s' %(gene_chr))
            #print('p_gene_chr %s' %(previous_gene_chr))
            #print('ann_chr %s' %(ann_chr))
                
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

            if (chrIndex[gene_chr] < chrIndex[ann_chr]): # if annotation moved to the next chromosome, but gene has not yet moved to the next chromosome
                annLineInd = annLineInd - 2
                line = linecache.getline(annFile, annLineInd)
                annLineInd +=1
                ann_line_count += 1
                fields = line.strip().split()

                ann_chr = fields[0]
                ann_start = int(fields[1])
                ann_end = int(fields[2])
                #print('marker 01')


            ''' 
            if this gene starts somewhere in the preivous genomic region, 
            I will go back in the annotation file to the beginning of annotation for the previous gene
            changed this th the next while loop - basically we will go back until the annotations start before the gene territory
            
            '''
            
            while (ann_start > gene_start) and (gene_chr == ann_chr): # in case of overlapping genes
                
                annLineInd = max(annLineInd - 5,1)
                line = linecache.getline(annFile, annLineInd)
                annLineInd +=1
                ann_line_count += 1
                fields = line.strip().split()

                ann_chr = fields[0]
                ann_start = int(fields[1])
                ann_end = int(fields[2])
                

            while ((ann_start < gene_end) and (gene_chr == ann_chr)) and gene_coverage < gene_length: #and geneMatWalkIndex < 160: # while we still have annotation before the gene

                if not((ann_end < gene_start) or ann_start > gene_end):
                
                    ''' We are in the genonimc region (with extension)'''

                    ''' Taking off for filling the matrices... '''
                        
                    adjusted_ann_start = max(0, ann_start - gene_start)
                    adjusted_ann_end = min(ann_end - gene_start, gene_end - gene_start)
                    adjusted_ann_length = adjusted_ann_end - adjusted_ann_start

                    ann_cluster = fields[3]
                    clusterInd = chromLabels.index(ann_cluster)

                    if gene_coverage < gene_length and adjusted_ann_length > 0:
                        new_coverage = min(adjusted_ann_length, gene_length - gene_coverage)
                        labelExpMat[cgi, clusterInd] += new_coverage
                        adjusted_ann_length = adjusted_ann_length - new_coverage
                        gene_coverage += new_coverage
                        #coverPromoter
 
                if gene_coverage < gene_length:
                    line = linecache.getline(annFile, annLineInd)
                    annLineInd +=1

                    ann_line_count += 1
                    fields = line.strip().split()

                    ann_chr = fields[0]
                    ann_start = int(fields[1])
                    ann_end = int(fields[2])
                

            cgi += 1 # next gene
            #print(cgi)
            previous_gene_end = gene_end
            previous_gene_chr = gene_chr
            
        linecache.clearcache()
        os.system('gzip %s' %(annFile))

        outputFile = annFolder  + 'enhancer_labelCoverage_chromHMM_hg20.pkl' 
        with open(outputFile, 'wb') as f:
            pickle.dump(labelExpMat, f)

        print(outputFile)
        print(enhAccessionList.index(accession))


########################################
# 4. Prediction
########################################

# get the fanTissueKyeList from 2, then fetch the column from the data for expression
# tissueNameID_map - 
# fanTissueKeyList

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
from sklearn.metrics import classification_report
from sklearn.metrics import accuracy_score

inputFileName = 'fantom5Mapping.pkl'
inputFile = dataFolder + inputFileName
print(inputFile)
with open(inputFile, 'rb') as f:
    allMappingInfo = pickle.load(f)

enhAccessionList = list(allMappingInfo.keys())

fanData_sorted = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/Fantom5Enh_data/human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix2_sorted.txt'

fandf = pd.read_csv(fanData_sorted, sep='\t') # sorted one

colList = fandf.columns.values.tolist()
print(len(colList))

fanArr = fandf[colList[2]].to_numpy()
print(fanArr.shape)
plt.plot(np.sort(np.log(fanArr+1)))
plt.show()

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> testing:
print('new Coors', newEnhCoors[30000])
print('index in the fandf', enhInds[30000])
print('old coors', assemblyMap2[newEnhCoors[30000]])
print('row ID in the fandf', fandf['Id'][enhInds[30000]])

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

aucListSegway = np.zeros(len(enhAccessionList))
for a, accession in enumerate(enhAccessionList):
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']

    inputFile = annotationFolder  + 'enhancer_labelCoverage_segway_states_hg20.pkl' 
    with open(inputFile, 'rb') as f:
        features = pickle.load(f)

    features = features[0:63985, :]

    fantomKey = allMappingInfo[accession][2]
    fantomTissue = fanTissueKeyList[fantomKey]
    tissueColKey = tissueNameID_map[fantomTissue].strip()
    print(allMappingInfo[accession], fantomTissue)

    enhExpTemp = np.array(fandf[tissueColKey])[0:63990]

    #TODO reorder the enhExpTemp according to the matrix
    enhExp = np.asarray([enhExpTemp[int(enhInds[x])] for x in range(len(newEnhCoors))])

    # >>>>>>>>>>>>>>>>>>>>
    # if we want to have a threshold for skipping some enhancers
    myThr = np.quantile(enhExp[enhExp>0], .3)
    print(np.sum(enhExp>0))
    myEnh = (enhExp > myThr) + (enhExp ==0)

    # <<<<<<<<<<<<<<<<<<<<
    
    #myEnh = enhExp>=0 # if we want to include everything (all)
    
    book = enhExp[myEnh]
    exp = book > 0

    #features = features[:, range(4)]
    X_train, X_test, y_train, y_test = train_test_split(features[myEnh, :], exp, test_size=0.25, random_state=11)
    logreg = LogisticRegression(random_state=11, class_weight='balanced')
    #logreg = LogisticRegression(random_state=11)

    logreg.fit(X_train, y_train)

    y_pred = logreg.predict(X_test)
    y_prob = logreg.predict_proba(X_test)
    cnf_matrix = metrics.confusion_matrix(y_test, y_pred) # a matrix
            
    print(cnf_matrix)
        
    target_names = ['expressed', 'not expressed']
    #print(classification_report(y_test, y_pred, target_names=target_names)) # an str
    #report =(classification_report(y_test, y_pred, target_names=target_names)) 

    auc = accuracy_score(y_test, y_pred) # a float
    print(auc)
    aucListSegway[a] = auc

kado = np.sort(np.array(aucListSegway))
plt.plot(kado)
plt.show()
# apply the transcriptomic pipeline to it.

plt.boxplot(kado[kado>0])
plt.boxplot(aucListSegway)
plt.ylim((.55, 1))
#plt.show()

# figFile = plotFolder + 'Segway_logisticRegression_enhancer_hg20_all.pdf'
figFile = plotFolder + 'Segway_logisticRegression_enhancer_hg20_.3Q.pdf'
print(figFile)
plt.savefig(figFile, bbox_inches='tight')
plt.close('all')


outputFile = dataFolder  + 'enhancerPredictionAUCs_segway_hg20_all.pkl' 
with open(outputFile, 'wb') as f:
    pickle.dump(aucListSegway, f)

inputFile =  dataFolder  + 'enhancerPredictionAUCs_segway_hg20.pkl'
with open(inputFile, 'rb') as f:
    enhancerAUCsSeg = pickle.load(f)

aucListChrom = np.zeros(len(enhAccessionList))
for a, accession in enumerate(enhAccessionList):

    annotation = allMeta[accession]
    annotationFolder = annotation['folder']
    print(a)

    if (annotation['chromFile'] == 'none'):
        continue

    inputFile = annotationFolder  + 'enhancer_labelCoverage_chromHMM_hg20.pkl' 
    with open(inputFile, 'rb') as f:
        features = pickle.load(f)
        
    features = features[0:63985, :]
    
    fantomKey = allMappingInfo[accession][2]
    fantomTissue = fanTissueKeyList[fantomKey]
    tissueColKey = tissueNameID_map[fantomTissue].strip()
    print(allMappingInfo[accession], fantomTissue)

    #enhExp = np.array(fandf[tissueColKey])[0:63990]

    enhExpTemp = np.array(fandf[tissueColKey])[0:63990]

    #TODO reorder the enhExpTemp according to the matrix
    enhExp = np.asarray([enhExpTemp[int(enhInds[x])] for x in range(len(newEnhCoors))])

    # >>>>>>>>>>>>>>>>>>>>
    # if we want to have a threshold for skipping some enhancers
    myThr = np.quantile(enhExp[enhExp>0], .3)
    print(np.sum(enhExp>0))
    myEnh = (enhExp > myThr) + (enhExp ==0)

    # <<<<<<<<<<<<<<<<<<<<
    
    #myEnh = enhExp>=0 # if we want to include everything (all)

    book = enhExp[myEnh]
    exp = book > 0

    X_train, X_test, y_train, y_test = train_test_split(features[myEnh, :], exp, test_size=0.25, random_state=11)
    logreg = LogisticRegression(random_state=11, class_weight='balanced', max_iter = 12000)
    #logreg = LogisticRegression(random_state=11)

    logreg.fit(X_train, y_train)

    y_pred = logreg.predict(X_test)
    y_prob = logreg.predict_proba(X_test)
    cnf_matrix = metrics.confusion_matrix(y_test, y_pred) # a matrix
            
    print(cnf_matrix)
        
    target_names = ['expressed', 'not expressed']
    #print(classification_report(y_test, y_pred, target_names=target_names)) # an str
    #report =(classification_report(y_test, y_pred, target_names=target_names)) 

    auc = accuracy_score(y_test, y_pred) # a float
    print(auc)
    aucListChrom[a] = auc

kchrom = np.sort(np.array(aucListChrom))
plt.plot(kchrom[kchrom>0])
plt.show()

plt.boxplot(kado[kado>0])
plt.boxplot(aucListChrom)
plt.ylim((.55,1))
plt.show()
# apply the transcriptomic pipeline to it.

#figFile = plotFolder + 'Chrom_logisticRegression_enhancer_hg20_all.pdf'
figFile = plotFolder + 'Chrom_logisticRegression_enhancer_hg20_.3Q.pdf'
print(figFile)
plt.savefig(figFile, bbox_inches='tight')
plt.close('all')

outputFile = dataFolder  + 'enhancerPredictionAUCs_chrom_hg20.pkl' 
with open(outputFile, 'wb') as f:
    pickle.dump(aucListChrom, f)

inputFile =  dataFolder  + 'enhancerPredictionAUCs_chrom_hg20.pkl'
with open(inputFile, 'rb') as f:
    enhancerAUCsChrom = pickle.load(f)

plt.boxplot(enhancerAUCsChrom)
plt.show()

plt.boxplot(enhancerAUCsSeg)
plt.show()




mymat = fandf.to_numpy()
book = mymat[:, 1:]
myAllSum = np.sum(book>0, axis=0)
myfourSum = np.sum(book>.4, axis=0)

plt.boxplot(myAllSum)
plt.show()

plt.boxplot(myfourSum)
plt.show()

ratio = [myfourSum[x]/myAllSum[x] for x in range(len(myAllSum))]
plt.boxplot(ratio)
plt.show()

########################################
# 5. Print the matching info for the supplement table
########################################

# load the mapping info file
inputFileName = 'fantom5Mapping.pkl'
inputFile = dataFolder + inputFileName
print(inputFile)
with open(inputFile, 'rb') as f:
    allMappingInfo = pickle.load(f)

suppTableFile = dataFolder + 'FantomMappingFile.tsv'
with open(suppTableFile, 'w') as f:
    f.write('SegwayAccession\tSegwayTissue\tFantomIndex\tFantomTissue\n')
    
    for ann in allMappingInfo:
        accession = ann
        f.write('%s\t%s\t%d\t%s\n')
    

########################################
# plot
########################################

plt.boxplot(kado[kado>0])
#plt.show()

figFile = plotFolder + 'Segway_logisticRegression_enhancer.pdf'
figFile = plotFolder + 'Chrom_logisticRegression_enhancer.pdf'
print(figFile)
plt.savefig(figFile, bbox_inches='tight')
plt.close('all')


########################################
# 6. Get the coverage prediction like transcriptomic
########################################

segwayAUCs = np.zeros((len(enhAccessionList), 11))
chromAUCs = np.zeros((len(enhAccessionList), 18))
for a, accession in enumerate(enhAccessionList):
    
    annotation = allMeta[accession]
    if (annotation['chromFile'] == 'none'):
        continue

    annotationFolder = annotation['folder']

    inputFile = annotationFolder  + 'enhancer_labelCoverage_segway_states.pkl'
    with open(inputFile, 'rb') as f:
        features = pickle.load(f)

    fantomKey = allMappingInfo[accession][2]
    fantomTissue = fanTissueKeyList[fantomKey]
    tissueColKey = tissueNameID_map[fantomTissue].strip()
    print(allMappingInfo[accession], fantomTissue)

    enhExp = np.array(fandf[tissueColKey])[0:63990]

    myEnh = (enhExp > .4) + (enhExp ==0)
    # myEnh = (enhExp > 0) + (enhExp ==0) 
    inEnh = enhExp[myEnh] 
    exp = inEnh > 0 # those that are expressed

    myFeatures = features[myEnh, :]

    for i in range(11):

        sib = np.argsort(myFeatures[:,i])[::-1]

        AUCCurve = np.zeros((len(sib), 2))
        xwalk = 0
        ywalk = 0
        area = 0 
        for j in range(len(sib)):
            if exp[sib[j]] == 0:
                xwalk +=1
                area+= ywalk
            else:
                ywalk +=1
                    
            AUCCurve[j, 0] = xwalk
            AUCCurve[j, 1] = ywalk

        auc = area / (xwalk*ywalk)
        print(auc)
        segwayAUCs[a,i] = auc


    inputFile = annotationFolder  + 'enhancer_labelCoverage_chromHMM.pkl' 
    with open(inputFile, 'rb') as f:
        features = pickle.load(f)

    fantomKey = allMappingInfo[accession][2]
    fantomTissue = fanTissueKeyList[fantomKey]
    tissueColKey = tissueNameID_map[fantomTissue].strip()
    print(allMappingInfo[accession], fantomTissue)

    enhExp = np.array(fandf[tissueColKey])[0:63990]

    myEnh = (enhExp > .4) + (enhExp ==0)
    # myEnh = (enhExp > 0) + (enhExp ==0) 
    inEnh = enhExp[myEnh] 
    exp = inEnh > 0 # those that are expressed

    myFeatures = features[myEnh, :]

    for i in range(18):

        sib = np.argsort(myFeatures[:,i])[::-1]

        AUCCurve = np.zeros((len(sib), 2))
        xwalk = 0
        ywalk = 0
        area = 0 
        for j in range(len(sib)):
            if exp[sib[j]] == 0:
                xwalk +=1
                area+= ywalk
            else:
                ywalk +=1
                    
            AUCCurve[j, 0] = xwalk
            AUCCurve[j, 1] = ywalk

        auc = area / (xwalk*ywalk)
        chromAUCs[a, i] = auc
        print(auc)


sib = np.sum(chromAUCs, axis=1)
kado = np.copy(chromAUCs)
kado = np.delete(kado, np.where(sib==0) ,axis=0)
sns.heatmap(kado)
plt.show()

sib = np.sum(chromAUCs, axis=1)
kado = np.copy(chromAUCs)
kado = np.delete(kado, np.where(sib==0) ,axis=0)
from scipy.stats import zscore
data_z = zscore(kado, axis=1)
df = pd.DataFrame(data_z, columns = chromLabels)
sns.heatmap(df)
figFile = plotFolder + 'promoterAUC_labelCoverage_chrom_zscore.pdf'
print(figFile)
plt.tight_layout()
plt.savefig(figFile, bbox_inches='tight')
plt.close('all')
plt.show()

sib = np.sum(segwayAUCs, axis=1)
kado = np.copy(segwayAUCs)
kado = np.delete(kado, np.where(sib==0) ,axis=0)
sns.heatmap(kado)
df = pd.DataFrame(kado, columns = segwayStates)
data_z = zscore(kado, axis=1)
df = pd.DataFrame(data_z, columns = segwayStates)
sns.heatmap(df)
plt.title('Segway FANTOM5 enhancer AUC predictions by state coverage')
figFile = plotFolder + 'promoterAUC_labelCoverage_segway_zscore.pdf'
print(figFile)
plt.tight_layout()
plt.savefig(figFile, bbox_inches='tight')
plt.close('all')
plt.show()

# they do very similar the two labels for the FANTOM5 enhancers
# TODO plot: perhaps get the disribution of stuff

########################################
# 7. Get the plot Max mentioned
########################################

# for a sample from the list:
# do this for accession index 1, 30, 50, 70, 100

# 7.1 Get the H3K4me1 data
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

accessionIndex = 70
accession = accessionList[accessionIndex]
#for accession in accessionList:
# get the accession for H3K4me1
annotation = allMeta[accession]
annotationFolder = annotation['folder']
# get the track file
assay_file = annotationFolder + 'trackname_assay.txt'
with open(assay_file, 'r') as assays:
    for line in assays:
        if line.split()[1] == 'H3K4me1':
            print(line)
            h3k4 = line.split()[0]
            print(accession)

# 7.2 get the bigwig file for chr20
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

file = '/Users/marjanfarahbod/Downloads/%s.bigWig' %(h3k4)
print(file)
import pyBigWig
endInd = 64000000
bw = pyBigWig.open(file)
vals = bw.values('chr20', 0, endInd)

vals=np.asarray(vals)
print(sum(vals>0))
plt.hist(vals)
plt.show()

annFile = annotation['bedFile']

# 7.3 prepare ann
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

if annFile.endswith('.gz'):
    os.system('gunzip %s' %(annFile))
else:
    os.system('gunzip %s.gz' %(annFile))

command = "grep -E 'chr20.*' %s" %(annFile)
print(command)
    
out = annotation['folder'] + 'chr20.bed'
f = open(out, 'w')
import subprocess
subprocess.run(command, shell=True, stdout=f)

os.system('gzip %s' %(annFile))

chr20File = annotation['folder'] + 'chr20.bed'

meanSig = np.zeros((int(endInd/100)))
for i in range(len(meanSig-1)):
    meanSig[i] = np.mean(vals[i*100:(i+1)*100])

plt.hist(meanSig)
plt.show()

# 7.4 Get the colors and labels right
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

mnemFile = annotationFolder + 'mnemonics_v04.txt'
label_term_mapping = {}
with open(mnemFile, 'r') as mnemonics:
    for line in mnemonics:
        #print(line)
        label = line.strip().split()[0]
        term = line.strip().split()[1]
        label_term_mapping[label] = term

colorStateMap = {'Bivalent' : np.asarray([189, 183, 107])/256,
                 'Transcribed' : np.asarray([0, 128, 0])/256,
                 'CTCF' : np.asarray([196, 225, 5])/256,
                 'K9K36' : np.asarray([102, 205, 170])/256,
                 'Enhancer' : np.asarray([255, 195, 77])/256,
                 'PromoterFlanking' : np.asarray([255, 68, 0])/256,
                 'ConstitutiveHet' : np.asarray([138, 145, 208])/256,
                 'FacultativeHet' : np.asarray([128, 0, 128])/256,
                 'Promoter' : np.asarray([255, 0, 0])/256,
                 'Quiescent' : np.asarray([190, 190, 190])/256,
                 'EnhancerLow' : np.asarray([255, 255, 0])/256}


segStates = np.zeros(len(meanSig))
i = 0
annLineInd = 1
while i < len(segStates):
    line = linecache.getline(chr20File, annLineInd)
    annStart = int(line.split('\t')[1])
    annEnd = int(line.split('\t')[2])
    annState = int(line.split('\t')[3].split('_')[0])
    step = int((annEnd - annStart) / 100)
    segStates[i:i+step] = annState
    i = i + step
    annLineInd +=1

segStatesToTerms = np.zeros(len(meanSig))
for i in range(len(segStates)):
    segStatesToTerms[i] = segwayStates.index(label_term_mapping[str(int(segStates[i]))])

# color the scatter hahahaha
colorString = []
for state in segStatesToTerms:
    colorString.append(colorStateMap[label_term_mapping[str(int(state))]])

plt.scatter(range(2000), meanSig[10000:12000], c=colorString[10000:12000])
plt.show()

# 7.5 Get the heatmap and the bar plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# meanSig varies from 0 to 4 and greater. It is fold change, so we care about 0, 1, 2, 3, 4
stateCoverage = np.zeros((5, 11))
meanSigCopy = meanSig.copy()
meanSigCopy[meanSigCopy >= 5] = 4.99
for i in range(5):
    #inds = np.where(meanSigCopy > i & meanSigCopy<=(i+1)[0])
    inds = np.where(np.logical_and(meanSigCopy >= i, meanSigCopy<(i+1)))[0]
    myStates = segStatesToTerms[inds]
    for j in range(11):
        stateCoverage[i, j] = np.sum(myStates==j)

stateCoverage = stateCoverage/np.sum(stateCoverage, axis=1)[:, np.newaxis]

sns.heatmap(stateCoverage)
plt.show()

for i in range(11):
    print(i)
    if i >0:
        plt.bar(range(5), stateCoverage[:,i], bottom=np.sum(stateCoverage[:, 0:i], axis=1), color=colorStateMap[segwayStates[i]])
    else:
        plt.bar(range(5), stateCoverage[:,i], color=colorStateMap[segwayStates[i]])
    print(stateCoverage[:,i])

titleString = ''
for i in range(5):
    val = np.sum((np.logical_and(meanSigCopy >= i, meanSigCopy<(i+1)))/len(meanSigCopy))
    print(i, val)
    #ratioList.append(np.sum((np.logical_and(meanSigCopy >= i, meanSigCopy<(i+1)))/len(meanSigCopy)))
    titleString = titleString + '   %.3f' %(val)

figFile = plotFolder + 'segwayLabelCoverage_H3K4me1_%s.pdf' %(accession)
print(figFile)
plt.title(titleString)
plt.tight_layout()
plt.savefig(figFile, bbox_inches='tight')
plt.close('all')
plt.show()


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#  Get the chrom plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

chromFile = annotation['chromFile']
print(chromFile)

fileName = chromFile.split('/')[-1][0:-3]
print(fileName)

os.system('gunzip %s' %(chromFile))

command = "grep -E 'chr20.*' %s" %(chromFile[0:-3])
print(command)
    
out = annotation['folder'] + 'chr20_chrom.bed'
f = open(out, 'w')
import subprocess
subprocess.run(command, shell=True, stdout=f)

os.system('gzip %s' %(chromFile[0:-3]))

chr20FileChrom = annotation['folder'] + 'chr20_chrom.bed'

i = 0
annLineInd = 1
chromStates = np.zeros(int(endInd/100))

pend = 0
while i < len(chromStates):
    #print(line)
    line = linecache.getline(chr20FileChrom, annLineInd)
    annStart = int(line.split('\t')[1])
    annEnd = int(line.split('\t')[2])
    annState = chromLabels.index((line.split('\t')[3]))
    #print(annState)
    if not(annStart == pend):
        step = int((annStart-pend) / 100)
        chromStates[i:i+step] = 18
        i = i+step
        print(annStart, pend)
        print(annLineInd)
        print(step)
        print('----------------')
        
    step = int((annEnd - annStart) / 100)
    chromStates[i:i+step] = annState
    i = i + step
    annLineInd +=1
    pend = annEnd

colorStateMapChrom = {'11' : np.asarray([189, 183, 107])/256,
                 '15' : np.asarray([0, 128, 0])/256,
                      '16' : np.asarray([0, 200, 0])/256,
                 '17' : np.asarray([102, 205, 170])/256,
                 '0' : np.asarray([255, 195, 77])/256,
                '1' : np.asarray([255, 195, 77])/256,
                      '2' :np.asarray([189, 183, 107])/256,
                      '3' : np.asarray([255, 195, 77])/256,
                      '4' : np.asarray([255, 195, 77])/256,
                 '12' : np.asarray([255, 68, 0])/256,
                      '13' : np.asarray([255, 68, 0])/256,
                      '14' : np.asarray([255, 68, 0])/256,
                 '6' : np.asarray([138, 145, 208])/256,
                 '8' : np.asarray([128, 0, 128])/256,
                 '9' : np.asarray([128, 0, 128])/256,                      
                 '10' : np.asarray([255, 0, 0])/256,
                 '7' : np.asarray([180, 180, 180])/256,
                      '5' : np.asarray([255, 255, 0])/256,
                      '18': np.asarray([0,0,0])/256}

colorStringChrom = []
for state in chromStates:
    colorStringChrom.append(colorStateMapChrom[str(int(state))])

plt.scatter(range(2000), meanSig[10000:12000], c=colorStringChrom[10000:12000])
plt.show()

stateCoverageChrom = np.zeros((5, 19))
meanSigCopy = meanSig.copy()
meanSigCopy[meanSigCopy >= 5] = 4.99
for i in range(5):
    #inds = np.where(meanSigCopy > i & meanSigCopy<=(i+1)[0])
    inds = np.where(np.logical_and(meanSigCopy >= i, meanSigCopy<(i+1)))[0]
    myStates = chromStates[inds]
    for j in range(19):
        stateCoverageChrom[i, j] = np.sum(myStates==j)

stateCoverageChrom = stateCoverageChrom/np.sum(stateCoverageChrom, axis=1)[:, np.newaxis]

sns.heatmap(stateCoverageChrom)
plt.show()


orderList = [10, 12, 13, 14, 0, 1, 3, 4, 5, 2, 11, 15, 16, 17, 8, 9, 6, 7, 18] # MAKE THE TODO LIST AND MAKE IT PRETTY
for i in range(19):
    print(i)
    if i >0:
        bottom = np.sum(stateCoverageChrom[:, orderList[0:i]], axis=1)
        plt.bar(range(5), stateCoverageChrom[:,orderList[i]], bottom=bottom, color=colorStateMapChrom[str(int(orderList[i]))])
        #plt.bar(range(5), stateCoverageChrom[:,i], bottom=np.sum(stateCoverageChrom[:, 0:i], axis=1), color=colorStateMapChrom[str(int(i))])
    else:
        plt.bar(range(5), stateCoverageChrom[:, orderList[i]], color=colorStateMapChrom[str(int(orderList[i]))])
    print(stateCoverageChrom[:,orderList[i]])

figFile = plotFolder + 'ChromLabelCoverage_H3K4me1_%s.pdf' %(accession)
print(figFile)
plt.tight_layout()
plt.title(titleString)
plt.savefig(figFile, bbox_inches='tight')
plt.close('all')
plt.show()

########################################
# 8. The plot Max mentioned
########################################

# the plot that Max mentioned I am going to take it to enhancer_geneExpression.py

########################################
### DRAFT
########################################
# for samples 0 and 5 from FANTOMLIST
# for a, accession in enumerate(enhAccessionList)

# get the accession for the histone mod
# get the packages and stuff

# plot the labels for each bin for chrom19

accession = enhAccessionList[0]
for accession in enhAccessionList:
    # get the accession for H3K4me1
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']
    # get the track file
    assay_file = annotationFolder + 'trackname_assay.txt'
    with open(assay_file, 'r') as assays:
        for line in assays:
            if line.split()[0] == 'ENCFF096GQC':
                print(line)
                print(accession)

# get the accession for H3K4me1
accession = enhAccessionList[125]
annotation = allMeta[accession]
annotationFolder = annotation['folder']
# get the track file
assay_file = annotationFolder + 'trackname_assay.txt'
with open(assay_file, 'r') as assays:
    for line in assays:
        print(line)


'ENCFF096GQC' # ann 125

'ENCFF028VXO' # ann 5

file = '/Users/marjanfarahbod/Downloads/ENCFF096GQC.bigWig'
import pyBigWig
bw = pyBigWig.open(file)
vals = bw.values('chr19', 3000000, 5000000)
vals = bw.values('chr20', 0, 43000000)

vals=np.asarray(vals)
print(sum(vals>0))
plt.hist(vals)
plt.show()

annFile = annotation['bedFile']

# prepare ann
if annFile.endswith('.gz'):
    os.system('gunzip %s' %(annFile))
else:
    os.system('gunzip %s.gz' %(annFile))

command = "grep -E 'chr19.*' %s" %(annFile)
command = "grep -E 'chr20.*' %s" %(annFile)
print(command)
    
out = annotation['folder'] + 'chr19.bed'
out = annotation['folder'] + 'chr20.bed'
f = open(out, 'w')
import subprocess
subprocess.run(command, shell=True, stdout=f)

os.system('gzip %s' %(annFile))

print('chr19 enhancer regions selected')

chr19File = annotation['folder'] + 'chr19.bed'
chr20File = annotation['folder'] + 'chr20.bed'

annLineInd = 1000
annLineInd = 0
line = linecache.getline(chr19File, annLineInd)
while int(line.strip().split('\t')[1]) < 3000000:
    annLineInd+=1
    line = linecache.getline(chr19File, annLineInd)

print(line)
# now get the signal

# for this file it starts at 3,000,000 so we are good

meanSig = np.zeros((20000))
meanSig = np.zeros((500000))
meanSig = np.zeros((430000))
for i in range(len(meanSig-1)):
    meanSig[i] = np.mean(vals[i*100:(i+1)*100])

plt.hist(meanSig)
plt.show()

# get the mnemonics
mnemFile = annotationFolder + 'mnemonics_v04.txt'
label_term_mapping = {}
with open(mnemFile, 'r') as mnemonics:
    for line in mnemonics:
        #print(line)
        label = line.strip().split()[0]
        term = line.strip().split()[1]
        label_term_mapping[label] = term

colorStateMap = {'Bivalent' : np.asarray([189, 183, 107])/256,
                 'Transcribed' : np.asarray([0, 128, 0])/256,
                 'CTCF' : np.asarray([196, 225, 5])/256,
                 'K9K36' : np.asarray([102, 205, 170])/256,
                 'Enhancer' : np.asarray([255, 195, 77])/256,
                 'PromoterFlanking' : np.asarray([255, 68, 0])/256,
                 'ConstitutiveHet' : np.asarray([138, 145, 208])/256,
                 'FacultativeHet' : np.asarray([128, 0, 128])/256,
                 'Promoter' : np.asarray([255, 0, 0])/256,
                 'Quiescent' : np.asarray([255, 255, 255])/256,
                 'EnhancerLow' : np.asarray([255, 255, 0])/256}


segStates = np.zeros(20000)
segStates = np.zeros(len(meanSig))
i = 0
while i < len(segStates):
    line = linecache.getline(chr20File, annLineInd)
    annStart = int(line.split('\t')[1])
    annEnd = int(line.split('\t')[2])
    annState = int(line.split('\t')[3].split('_')[0])
    step = int((annEnd - annStart) / 100)
    segStates[i:i+step] = annState
    i = i + step
    annLineInd +=1

segStatesToTerms = np.zeros(20000)
segStatesToTerms = np.zeros(len(meanSig))
for i in range(len(segStates)):
    segStatesToTerms[i] = segwayStates.index(label_term_mapping[str(int(segStates[i]))])

# color the scatter hahahaha
colorString = []
for state in segStates:
    colorString.append(colorStateMap[label_term_mapping[str(int(state))]])

plt.scatter(range(1000), meanSig[0:1000], c=colorString[0:1000])
plt.show()


# meanSig varies from 0 to 4 and greater. It is fold change, so we care about 0, 1, 2, 3, 4
stateCoverage = np.zeros((5, 11))
meanSigCopy = meanSig.copy()
meanSigCopy[meanSigCopy > 5] = 5
for i in range(5):
    #inds = np.where(meanSigCopy > i & meanSigCopy<=(i+1)[0])
    inds = np.where(np.logical_and(meanSigCopy > i, meanSigCopy<=(i+1)))[0]
    myStates = segStatesToTerms[inds]
    for j in range(11):
        stateCoverage[i, j] = np.sum(myStates==j)

stateCoverage = stateCoverage/np.sum(stateCoverage, axis=1)[:, np.newaxis]

sns.heatmap(stateCoverage)
plt.show()

for i in range(11):
    print(i)
    if i >0:
        plt.bar(range(5), stateCoverage[:,i], bottom=np.sum(stateCoverage[:, 0:i], axis=1), color=colorStateMap[segwayStates[i]])
    else:
        plt.bar(range(5), stateCoverage[:,i], color=colorStateMap[segwayStates[i]])
    print(stateCoverage[:,i])

plt.show()

a
# do the same for chromfile
########################################
chromFile = annotation['chromFile']
print(chromFile)

fileName = chromFile.split('/')[-1][0:-3]
print(fileName)

os.system('gunzip %s' %(chromFile))


command = "grep -E 'chr20.*' %s" %(chromFile[0:-3])
print(command)
    
out = annotation['folder'] + 'chr20_chrom.bed'
f = open(out, 'w')
import subprocess
subprocess.run(command, shell=True, stdout=f)

chr20FileChrom = annotation['folder'] + 'chr20_chrom.bed'

annLineInd = 1000
line = linecache.getline(chr19File, annLineInd)
while int(line.strip().split('\t')[1]) < 3000000:
    annLineInd+=1
    line = linecache.getline(chr19File, annLineInd)

annLineInd = 1256
i = 0
chromStates = np.zeros(20000)
chromStates = np.zeros(430000)
line = linecache.getline(chr19FileChrom, annLineInd)
annStart = 3000000
annEnd = int(line.split('\t')[2])
annState = chromLabels.index((line.split('\t')[3]))
step = int((annEnd - annStart) / 100)
chromStates[i:i+step] = annState
i = i + step
annLineInd +=1

pend = 0
while i < len(chromStates):
    #print(line)
    line = linecache.getline(chr20FileChrom, annLineInd)
    annStart = int(line.split('\t')[1])
    annEnd = int(line.split('\t')[2])
    annState = chromLabels.index((line.split('\t')[3]))
    #print(annState)
    if not(annStart == pend):
        step = int((annStart-pend) / 100)
        chromStates[i:i+step] = 18
        i = i+step
        print(annStart, pend)
        print(annLineInd)
        print(step)
        print('----------------')
        
    step = int((annEnd - annStart) / 100)
    chromStates[i:i+step] = annState
    i = i + step
    annLineInd +=1
    pend = annEnd


colorStateMapChrom = {'11' : np.asarray([189, 183, 107])/256,
                 '15' : np.asarray([0, 128, 0])/256,
                      '16' : np.asarray([0, 200, 0])/256,
                 '17' : np.asarray([102, 205, 170])/256,
                 '0' : np.asarray([255, 195, 77])/256,
                '1' : np.asarray([255, 195, 77])/256,
                      '2' :np.asarray([189, 183, 107])/256,
                      '3' : np.asarray([255, 195, 77])/256,
                      '4' : np.asarray([255, 195, 77])/256,
                 '12' : np.asarray([255, 68, 0])/256,
                      '13' : np.asarray([255, 68, 0])/256,
                      '14' : np.asarray([255, 68, 0])/256,
                 '6' : np.asarray([138, 145, 208])/256,
                 '8' : np.asarray([128, 0, 128])/256,
                 '9' : np.asarray([128, 0, 128])/256,                      
                 '10' : np.asarray([255, 0, 0])/256,
                 '7' : np.asarray([180, 180, 180])/256,
                      '5' : np.asarray([255, 255, 0])/256,
                      '18': np.asarray([0,0,0])/256}

colorStringChrom = []
for state in chromStates:
    colorStringChrom.append(colorStateMapChrom[str(int(state))])


plt.scatter(range(1000), meanSig[0:1000], c=colorStringChrom[0:1000])
plt.show()

stateCoverageChrom = np.zeros((5, 19))
meanSigCopy = meanSig.copy()
meanSigCopy[meanSigCopy > 5] = 5
for i in range(5):
    #inds = np.where(meanSigCopy > i & meanSigCopy<=(i+1)[0])
    inds = np.where(np.logical_and(meanSigCopy > i, meanSigCopy<=(i+1)))[0]
    myStates = chromStates[inds]
    for j in range(19):
        stateCoverageChrom[i, j] = np.sum(myStates==j)

stateCoverageChrom = stateCoverageChrom/np.sum(stateCoverageChrom, axis=1)[:, np.newaxis]

sns.heatmap(stateCoverageChrom)
plt.show()

for i in range(19):
    print(i)
    if i >0:
        plt.bar(range(5), stateCoverageChrom[:,i], bottom=np.sum(stateCoverageChrom[:, 0:i], axis=1), color=colorStateMapChrom[str(int(i))])
    else:
        plt.bar(range(5), stateCoverageChrom[:,i], color=colorStateMapChrom[str(int(i))])
    print(stateCoverageChrom[:,i])

plt.show()


