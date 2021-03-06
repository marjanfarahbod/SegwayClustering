# Scaling up the transcription data comparison

################################################################################
################################################################################
##### TODO: fix the gene thing
################################################################################
################################################################################

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
# TODO: not sure if it works
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

########################################
# 1. Specifics
########################################

# list of annotation folders
sampleList = os.listdir(dataFolder + dataSubFolder)
print(len(sampleList))

########################################
# 2. Pre processings 
########################################

# 2.0 annotation folder prep

# NOTE: we have a class for annotations, we are not using it at the moment

annMeta = {}
c = 1
for annFolder in sampleList:
    print(c)
    c += 1

    sampleFolderAdd = dataFolder + dataSubFolder + annFolder + '/'
    bedFile = 'none'
    expressionFile = 'none'

    files = os.listdir(dataFolder + dataSubFolder + annFolder)

    for file in files:
        if (file.endswith('.bed')) and file.startswith('ENC') and len(file) < 18:
            
            bedFile = file
        
        if file.endswith('bed.gz'):
            fileAdd = sampleFolderAdd + file
            os.system('gunzip %s' %(fileAdd))
            bedFile = file[0:-3]
        
        if file.startswith('preferred_default'):
            expressionFile = file
            
    annMeta[annFolder] = {"bedFile": bedFile, "expressionFile": expressionFile}

outputFile = dataFolder + dataSubFolder + 'metaInfo.pkl'
with open(outputFile, 'wb') as f:
    pickle.dump(annMeta, f)

# 2.1. preping expression data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

annAccessionList = list(annMeta.keys())
# annAccession = annAccessionList[3]

for annAccession in annAccessionList:

    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
    print(sampleFolderAdd)

    expFileName = annMeta[annAccession]['expressionFile']
    if (expFileName != 'none'):

#        mistakeFile = sampleFolderAdd + 'geneExp_dict_ENCFF776YSH.pkl'
#        print(mistakeFile)
#        os.remove(mistakeFile)
#        print('doneRemove')

        print('doing the expression')
        expFile = sampleFolderAdd + expFileName
        print(expFile)
        expression = {}
        with open(expFile, 'r') as file:
            line = file.readline() # there is one header line

            for line in file:
                fields = line.strip().split()
                geneID = fields[0]
                transcriptID = fields[1]

                if geneID in expression:
                    expression[geneID] += np.log10(float(fields[5]) + 1)
                else:
                    expression[geneID] = np.log10(float(fields[5]) + 1)

        # saving expression dict            
        expAccession = re.split('_|\.', expFileName)[2]
        outputFile = sampleFolderAdd + 'geneExp_dict_' + expAccession + '.pkl'
        print('printing this file %s' %(outputFile))
         
        with open(outputFile, 'wb') as f:
            pickle.dump(expression, f)


# loading expression dict
expAccession = re.split('_|\.', expFileName)[2]
inputFile = sampleFolderAdd + 'geneExp_dict_' + expAccession + '.pkl'
with open(inputFile, 'rb') as pickledFile:
    expression = pickle.load(pickledFile)


# 2.2. preping annotation data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# the function gets annotation bed file and documents the report on it for length and etc. for clusters, we also get the length distribution which is useful for classifier performance correction

# fetch the length distribution data from length_distribution

## >>>>>>>>>>>>>>>>>>>>>>>> code for one file
annAccession = annAccessionList[0]
sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
bedFileName = annMeta[annAccession]['bedFile']

bedFileAdd = sampleFolderAdd + bedFileName

annotation_generalInfo_clusters(bedFileAdd)

inputFile = sampleFolderAdd + 'ENCFF172LXH_annotationSummary.pkl'
with open(inputFile, 'rb') as pickledFile:
    annSum = pickle.load(pickledFile)

# the function
def annotation_generalInfo_clusters(bedFileAdd):

    splitFileName = bedFileAdd.split('/')
    bedFileName = splitFileName[-1].split('.')[0]
    index = bedFileAdd.index(bedFileName)
    outputFolder = bedFileAdd[0:index]

    previous_class = ''

    clusters = {}
    classes = {}
    c = 0
    with open(bedFileAdd, 'r') as annotations:

        # annotations have no header
        # header = annotations.readline()

#        f = open(bedFileAdd, 'r')
#        line = f.readline()
        
        for line in annotations:
            c += 1
            fields = line.strip().split()

            # doing the clusters first
            if fields[3] in clusters.keys():

                clusters[fields[3]].bp_count += int(fields[2]) - int(fields[1])
                clusters[fields[3]].region_count += 1
                
            else:
                
                called = fields[3]
                biolabel = fields[3].split('_')[1]
                cluster = fields[3].split('_')[0]
                color = fields[8]
                bp_count = int(fields[2]) - int(fields[1])
                region_count = 1
                region_dist = 1
                clusters[fields[3]] = Annotation(called, biolabel, cluster, color, bp_count, region_count, region_dist)


            # doing the class
            if previous_class == fields[3].split('_')[1]:
                classes[previous_class].bp_count +=  int(fields[2]) - int(fields[1])
            else:
                current_class = fields[3].split('_')[1]
                if current_class in classes.keys():
                    classes[current_class].bp_count +=  int(fields[2]) - int(fields[1])
                    classes[current_class].region_count += 1
                else:
                    clusterList = [] # not filling it now, it can be filled later using annotations 
                    biolabel = current_class
                    color = fields[8]
                    bp_count = int(fields[2]) - int(fields[1])
                    region_count = 1
                    region_dist = 1
                    classes[biolabel] = AnnotationClass(biolabel, clusterList, color, bp_count, region_count, region_dist)

            previous_class = fields[3].split('_')[1]

            
    # filling up the cluster distribution from distribution file
    
    segmentSizesFile = outputFolder + 'segment_sizes.tab.txt'
    with open(segmentSizesFile, 'r') as inputFile:
        lines = inputFile.readlines()
    
    clusterList = list(clusters.keys())
    for cluster in clusterList:
        cluster_number = int(cluster.split('_')[0])
        fields = lines[cluster_number + 2].strip().split('\t')
        info = {'num.segs': int(fields[1]), 'mean.len': float(fields[2]), 'median.len': float(fields[3]), 'stdev.len': float(fields[4]), 'num.bp': float(fields[5]), 'frac.bp': float(fields[6])}
        clusters[cluster].region_dist = info

    annotationSummary = {"classes": classes, "clusters": clusters}
    outputFile = outputFolder + bedFileName + '_annotationSummary.pkl'
    with open(outputFile, 'wb') as f:
        pickle.dump(annotationSummary, f)

    print('annotation summary saved in %s' %(outputFile))

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# code for all files:

for annAccession in annAccessionList:

    print(annAccession)
    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
    bedFileName = annMeta[annAccession]['bedFile']

    bedFileAdd = sampleFolderAdd + bedFileName

    files = os.listdir(sampleFolderAdd)
    for file in files:
        if file.endswith('filteredSorted.bed'):
            rmAdd = sampleFolderAdd + file
            os.remove(rmAdd)

    #annotation_generalInfo_clusters(bedFileAdd)

    # sort the bed file
    sortedBedFile = sampleFolderAdd + 'sortedBedFile.bed'
    os.system('sort -V -k1,1 -k2,2 %s > %s' %(bedFileAdd, sortedBedFile))

    # filter the bed file
    bedFileAccession = bedFileName.split('.')[0]
    filteredSortedBedFile = sampleFolderAdd + bedFileAccession +'_filteredSorted.bed'
    os.system("awk -F'\t' '$1 !~ \"_\"' %s > %s" %(sortedBedFile, filteredSortedBedFile))

    os.remove(sortedBedFile)


#############################################################    
# 3. Comparing the annotation labels to the genomic regions and transcriptomic data
#############################################################


#classList = list(classes.keys())
classList = segwayLabels
classCount = len(classList) # this is the default for our Segway annotation annotations - 9
classListIndMap = {}
for i in range(len(classList)):
   classListIndMap[classList[i]] = i

extension = 3000 # count of basepairs monitored before and after the gene coordinates

annAccessionCount = 0
for annAccession in annAccessionList:

    # for each annotation
    if (annMeta[annAccession]['expressionFile'] != 'none'):

        #  annAccession = annAccessionList[5]
        #  print(annMeta[annAccession]['expressionFile'] != 'none')

        print(annAccession)
        annAccessionCount += 1
        print(annAccessionCount)

        expFileName = annMeta[annAccession]['expressionFile']
        expAccession = re.split('_|\.', expFileName)[2]
        inputFile = dataFolder + dataSubFolder + annAccession + '/geneExp_dict_' + expAccession + '.pkl'
        print(inputFile)
        with open(inputFile, 'rb') as pickledFile:
            expression = pickle.load(pickledFile)


        #annFile = dataFolder + dataSubFolder + annAccession + '/' + annMeta[annAccession]['bedFile']
        annFile = dataFolder + dataSubFolder + annAccession + '/' + annMeta[annAccession]['bedFile'].split('.')[0] + '_filteredSorted.bed' # annotation file needed filtering and sorting
        
        bedFileName = annMeta[annAccession]['bedFile'].split('.')[0]
        #bedFileName = 'ENCFF153PCI'
        inputFile = dataFolder + dataSubFolder + annAccession + '/' + bedFileName + '_annotationSummary.pkl'
        with open(inputFile, 'rb') as pickledFile:
            summaryAnnotation = pickle.load(pickledFile)

        annLineCount = int(sp.getoutput('wc -l %s' %(annFile)).split()[0])

        #print(summaryAnnotation['classes'])
        #print(summaryAnnotation['clusters'])

        classes = summaryAnnotation['classes']
        clusters = summaryAnnotation['clusters']

        # TODO: this part should go to the annotation processing
        clusterList = list(clusters.keys())
        sortedClusterList = []
        #print(segwayLabels)
        for label in segwayLabels:
            #print(label)
            for item in clusterList:
                #print(item)
                if item.split('_')[1] == label:
                    sortedClusterList.append(item)


        # getting the cluster list index map
        clusterList = list(clusters.keys()) # this is per sample
        clusterCount = len(clusterList) # this could vary for each annotation file
        print(sortedClusterList)
        clusterListIndMap = {}
        for i in range(len(clusterList)):
            clusterListIndMap[sortedClusterList[i]] = i


        cgi = 0 # walks on the geneIDList
        ann_start = 0 # the start of the current annotation
        ann_end = 0 # the end of the current annotation
        ann_chr = 'chr'
        ann_line_count = 0 # this is just to check the progress through the annotation file
        previous_class = ''


        clusterMats = [np.zeros((clusterCount, 160)),np.zeros((clusterCount, 160)),np.zeros((clusterCount, 160))]
        classMats = [np.zeros((classCount, 160)),np.zeros((classCount, 160)),np.zeros((classCount, 160))]

        # this is to use for the negative strand genes - not using now
        #tempClusterMats = [np.zeros((clusterCount, 160)),np.zeros((clusterCount, 160)),np.zeros((clusterCount, 160))]
        #tempClassMats = [np.zeros((classCount, 160)),np.zeros((classCount, 160)),np.zeros((classCount, 160))]

        annLineInd = 1 # this keeps the annotation line index, the first line is empty, thus the index starts at 1
        genomic_region_start_annLineInd = 0
        #firstGenomicAnn = True

        previous_gene_chr = 'chr1'
        previous_extension_end = 0 # doesn't matter since chr condition is never true until the first gene is processed
        previous_ann_chr = 'chr'


        #while cgi < len(geneIDList) and annLineInd < annLineCount: # modify the condition for the test runs
        while cgi < 26017: # >>>>>>>>>> MAIN
        #while cgi < 5655:# >>>>>>>>>> TEST
            #print('cgi')
            #print(cgi)

            'we are in the gene territory, we are walking on the genes'
            #firstGenomicAnn = True
   
            geneID = geneIDList[cgi]
            gene_chr = geneList[geneIDList[cgi]].chrom
            #print(geneID)

            gene_strand = geneList[geneID].strand
            gene_start = geneList[geneID].start
            gene_end = geneList[geneID].end
            gene_length = gene_end - gene_start
            gene_length_unit = int(gene_length/100)
            if gene_length_unit == 0:
                #print('gene smaller than 100bp, skipping it')
                cgi = cgi + 1
                continue
            
            gene_length_last_unit = gene_length - (99* gene_length_unit)
            # TODO: something to fix: the gene_length_last_unit for negative strand versus positive strand

            extension_start = gene_start - extension
            extension_end = gene_end + extension

            #print('%ss, %s, %s' %(gene_chr, extension_start, extension_end))

            ''' 
            if this gene starts somewhere in the preivous genomic region, 
            I will go back in the annotation file to the beginning of annotation for the previous gene
            changed this th the next while loop - basically we will go back until the annotations start before the gene territory

            '''
            #if (gene_chr == previous_gene_chr) and (previous_extension_end > extension_start):
            #   annLineInd = genomicRegionStartAnnLineInd

            ''' picking the label/class matrix based on the gene expression level'''
            #TODO catch exception for when geneID is not in expression
            gene_exp = expression[geneID]
            if gene_exp == 0:
                expMatInd = 0
            elif gene_exp > 1:
                expMatInd = 2
            else:
                expMatInd = 1

            #print('expMatInd')
            #print(expMatInd) # >>>>> test
       
            geneMatWalkIndex = 0
            previous_fill = 0 # at the begining of each gene, annotations are full

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
                print('gene chr not equal to previous gene chr')
                print('gene_chr %s' %(gene_chr))
                print('p_gene_chr %s' %(previous_gene_chr))
                print('ann_chr %s' %(ann_chr))
                while (ann_chr != gene_chr): # move on in the annotation until we reach to the next chromosome
                    #print('reading annotations until it is equal')
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

            if (ann_chr != previous_ann_chr): # in case of chromosome change because annotation moved to the next chromosome
                if (ann_chr != gene_chr): # if annotation moved to the next chromosome, but gene has not yet moved to the next chromosome

        figFile = dataFolder + dataSubFolder + annAccession + '/exp2cluster_heatmap.pdf'
        plt.savefig(figFile)
        plt.close('all')

        book = pd.DataFrame(classMats[2])
        sns.heatmap(book)
        figFile = dataFolder + dataSubFolder + annAccession + '/exp2class_heatmap.pdf'
        plt.savefig(figFile)
        plt.close('all')

        book = pd.DataFrame(clusterMats[0])
        sns.heatmap(book)
        figFile = dataFolder + dataSubFolder + annAccession + '/exp0cluster_heatmap.pdf'
        plt.savefig(figFile)
        plt.close('all')

        book = pd.DataFrame(classMats[0])
        sns.heatmap(book)
        figFile = dataFolder + dataSubFolder + annAccession + '/exp0class_heatmap.pdf'
        plt.savefig(figFile)
        plt.close('all')


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# getting the enrichment
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



annAccession = 'ENCSR427OOB'

# expOverlapMats = {"clusterMats": clusterMats, "clussMats": classMats}
inputFile = dataFolder + dataSubFolder + annAccession + '/' + 'defaultExp_5kg_expSummary.pkl'
with open(inputFile, 'rb') as pickledFile:
    expOverlapMats = pickle.load(pickledFile)

book = pd.DataFrame(expOverlapMats['clussMats'][2])
sns.heatmap(book)
plt.show()
plt.close('all')

# they are the right matrix, get the matrix and the enrichment

myMat = expOverlapMats['clussMats'][2]

sumCol = book.sum(axis = 0)
divValue = max(sumCol)
bookDiv = book.div(divValue) # enrichment in the plot


bedFileName = annMeta[annAccession]['bedFile'].split('.')[0]
inputFile = dataFolder + dataSubFolder + annAccession + '/' + bedFileName + '_annotationSummary.pkl'
with open(inputFile, 'rb') as pickledFile:
    summaryAnnotation = pickle.load(pickledFile)

annClassList = list(summaryAnnotation['classes'].keys())

print(summaryAnnotation['classes'][annClassList[0]])

totalbp = 0
for ann in annClassList:
    totalbp += summaryAnnotation['classes'][ann].bp_count

annClassEnrichment = np.zeros(len(annClassList))
for i,ann in enumerate(annClassList):
    annClassEnrichment[i] = summaryAnnotation['classes'][ann].bp_count /totalbp

myMat = expOverlapMats['clussMats'][2] + 1
sumCol = myMat.sum(axis = 0)
maxVal = sumCol.max()
myMatDiv = myMat / maxVal

# divide the relevant classes with the relevant value
myMatEnr = myMadDiv
for i,ann in enumerate(annClassList):

    print(i)
    print(ann)
    row_number = classListIndMap[ann]
    annDivValue = annClassEnrichment[i]

    for j in range(myMatDiv.shape[1]):
        myMatEnr[row_number][j] = myMatDiv[row_number][j] / annDivValue


plt.plot(list(range(0,160)), np.log10(myMatDiv[3][:]))
plt.show()
  
x = np.arange(0, 160, 1)
y1 = np.log10(myMatEnr[3][:])
y2 = np.log10(myMatEnr[4][:])
y3 = np.log10(myMatEnr[5][:])
y4 = np.log10(myMatEnr[6][:])
y5 = np.log10(myMatEnr[0][:])

x = np.arange(0, 160, 1)
y1 =  (myMatDiv[3][:])
y2 =  (myMatDiv[4][:])
y3 =  (myMatDiv[5][:])
y4 =  (myMatDiv[6][:])
y5 =  (myMatDiv[0][:])         
y6 =  (myMatDiv[2][:])         

fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6,1, sharex=True)

ax1.fill_between(x, 0, y1)
ax1.set_ylim(0,1)

ax2.fill_between(x, 0, y2)
ax2.set_ylim(0,1)

ax3.fill_between(x, 0, y3)
ax3.set_ylim(0,1)

ax4.fill_between(x, 0, y4)
ax4.set_ylim(0,1)

ax5.fill_between(x, 0, y5)
ax5.set_ylim(0,1)

ax6.fill_between(x, 0, y6)
ax6.set_ylim(0,1)


figFile = dataFolder + dataSubFolder + annAccession + '/TPERQF_distPlot_exp2.pdf'
plt.savefig(figFile)
plt.show()


# >>>>>>>>>>>> copy for the expression zero

myMat = expOverlapMats['clussMats'][0] + 1
sumCol = myMat.sum(axis = 0)
maxVal = sumCol.max()
myMatDiv = myMat / maxVal

# divide the relevant classes with the relevant value
myMatEnr = myMadDiv
for i,ann in enumerate(annClassList):

    print(i)
    print(ann)
    row_number = classListIndMap[ann]
    annDivValue = annClassEnrichment[i]

    for j in range(myMatDiv.shape[1]):
        myMatEnr[row_number][j] = myMatDiv[row_number][j] / annDivValue


plt.plot(list(range(0,160)), np.log10(myMatDiv[3][:]))
plt.show()
  
x = np.arange(0, 160, 1)
y1 = np.log10(myMatEnr[3][:])
y2 = np.log10(myMatEnr[4][:])
y3 = np.log10(myMatEnr[5][:])
y4 = np.log10(myMatEnr[6][:])
y5 = np.log10(myMatEnr[0][:])

x = np.arange(0, 160, 1)
y1 =  (myMatDiv[3][:])
y2 =  (myMatDiv[4][:])
y3 =  (myMatDiv[5][:])
y4 =  (myMatDiv[6][:])
y5 =  (myMatDiv[0][:])
y6 =  (myMatDiv[2][:])


fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6,1, sharex=True)

ax1.fill_between(x, 0, y1)
ax1.set_ylim(0,1)

ax2.fill_between(x, 0, y2)
ax2.set_ylim(0,1)

ax3.fill_between(x, 0, y3)
ax3.set_ylim(0,1)

ax4.fill_between(x, 0, y4)
ax4.set_ylim(0,1)

ax5.fill_between(x, 0, y5)
ax5.set_ylim(0,1)

ax6.fill_between(x, 0, y6)
ax6.set_ylim(0,1)


figFile = dataFolder + dataSubFolder + annAccession + '/TPERQF_distPlot_exp0.pdf'
plt.savefig(figFile)
plt.show()
        
plotMat = pd.DataFrame(myMatDiv)
sns.heatmap(plotMat)

plt.show()

plt.plot()

