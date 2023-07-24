# so here we are.

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
# 91% of enhancers have appear in <350 tissues.
# 74% appear in <200 tissues

# tissue recognition
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

# read the enhancers into a tupple
annLineCount = int(sp.getoutput('wc -l %s' %(fanData_sorted)).split()[0])
enh_dict = {}
with open(fanData_sorted, 'r') as f:
    for i in range(0,annLineCount-2):
        ml = linecache.getline(fanData_sorted, i+2)[0:30].split('\t')[0]
        sml = ml.split(':')
        chr = sml[0]
        s = int(sml[1].split('-')[0])
        e = int(sml[1].split('-')[1])
        enh_dict[i] = (chr, s, e)

###### that function for Segway

# get the coverage of the region for each of the labels:
for accession in enhAccessionList[90:]:
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

    labelExpMat = np.zeros((len(enh_dict), clusterCount))
    
    annLineInd = 1 # this keeps the annotation line index, the first line is empty, thus the index starts at 1

    previous_gene_chr = 'chr1'
    previous_extension_end = 0 # doesn't matter since chr condition is never true until the first gene is processed
    previous_ann_chr = 'chr1'


    while cgi < 63990: # >>>>>>>>>> MAIN

        'we are in the gene territory, we are walking on the genes'
        # enh_dict
        gene_chr = enh_dict[cgi][0]
            
        gene_start = enh_dict[cgi][1]
        gene_end = enh_dict[cgi][2]

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

        if (gene_chr != previous_gene_chr): # in case of chromosome change because gene moved to the next chromosome
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

                
        if (ann_chr != gene_chr): # if annotation moved to the next chromosome, but gene has not yet moved to the next chromosome
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
        previous_gene_end = gene_end
        previous_gene_chr = gene_chr
            

    linecache.clearcache()
    os.system('gzip %s' %(annFile))

    outputFile = annFolder  + 'enhancer_labelCoverage_segway.pkl' 
    with open(outputFile, 'wb') as f:
        pickle.dump(labelExpMat, f)

    print(outputFile)
    print(enhAccessionList.index(accession))
    
    
##### that function for chromhmm

book =  'EnhA1 EnhA2 EnhBiv EnhG1 EnhG2 EnhWk Het Quies ReprPC ReprPCWk TssA TssBiv TssFlnk TssFlnkD TssFlnkU Tx TxWk ZNF/Rpts'
chromLabels = book.split()

for accession in enhAccessionList[82:]:

    annotation = allMeta[accession]

    if not(annotation['chromFile'] == 'none'):
    
        annFolder = annotation['folder']

        chmmFileName = annotation['chromFile'].split('/')[-1]
        annFile = annFolder + 'sorted_' + chmmFileName
        #annFile = annFolder + chmmFileName # for 105 runs

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

        while cgi < 63990: # >>>>>>>>>> MAIN
        
            'we are in the gene territory, we are walking on the genes'
            #enh_dict
            gene_chr = enh_dict[cgi][0]
                
            gene_start = enh_dict[cgi][1]
            gene_end = enh_dict[cgi][2]

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

            if (gene_chr != previous_gene_chr): # in case of chromosome change because gene moved to the next chromosome
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

            if (ann_chr != gene_chr): # if annotation moved to the next chromosome, but gene has not yet moved to the next chromosome
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

        outputFile = annFolder  + 'enhancer_labelCoverage_chromHMM.pkl' 
        with open(outputFile, 'wb') as f:
            pickle.dump(labelExpMat, f)

        print(outputFile)
        print(enhAccessionList.index(accession))

# for each of my own tissues, get the list of tissues from the tissue file that has "any" of the words in the name in it. print it, pick it.

# once the list of tissues is made, for each tissue fetch the enhancer regions. We need the coverage of the labels for that region (much like the genes)

# apply the transcriptomic pipeline to it. 

# first task: get my sample tisssues, and search the directory of the tissues for my samples

# how many tissues do they have, what is the coverage, what is the overlap of the regions


