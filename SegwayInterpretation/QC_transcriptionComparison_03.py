# Scaling up the transcription data comparison - the previous version of this code QC_transcriptionComparison_02.py is messed up - I made this code using git history
################################################################################
################################################################################
# TODO: fetch the length distribution data from length_distribution section 2.2?
################################################################################
################################################################################
#
# 0. Initials 
# 1. Specifics
# 2. Pre processings
## 2.0 annotation folder prep
## 2.1. preping expression data
## 2.2. preping annotation data
# 3. Comparing the annotation labels to the genomic regions and transcriptomic data - OBSOLETE
# 4. Preprocessing transcriptomic data
# 5. process the annotation file
# 6. use the function to get the transcription enrichment for Segway

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

# 0.1:
geneIDsInfoFile = dataFolder  + 'geneIDsInfo.tsv'

myGenes = list(geneList.keys())
with open(geneIDsInfoFile, 'w') as f:
    f.write('geneSymbol\tENS_ID\ttype\tchr\tstart\tend\tstrand\n')
    for gene in myGenes:
        f.write('%s\t%s\t%s\t%s\t%d\t%d\t%s\n' %(geneList[gene].name,
                                                 geneList[gene].ENS_ID,
                                                 geneList[gene].gtype,
                                                 geneList[gene].chrom,
                                                 geneList[gene].start,
                                                 geneList[gene].end,
                                                 geneList[gene].strand))


########################################
# 1. Specifics
########################################

# Getting list of annotation folders
inputFile = dataFolder + dataSubFolder + 'metaInfo.pkl'
with open(inputFile, 'rb') as f:
    annMeta = pickle.load(f)

sample_count = len(annMeta)

annAccessionList = list(annMeta.keys())  # gettign list of annotationAccession
annAccession = annAccessionList[104]
print(annAccession)

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

inputFile = dataFolder + dataSubFolder + 'metaInfo.pkl'
with open(inputFile, 'rb') as f:
    annMeta = pickle.load(f)

 # 2.1. preping expression data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

annAccessionList = list(annMeta.keys())

## Please see QC_transcriptionComparison_util.py for function transcriptFile_preprocessing()
for annAccession in annAccessionList:
    
    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
    print(sampleFolderAdd)
    
    expFileName = annMeta[annAccession]['expressionFile']
    if (expFileName != 'none'):
        
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
# TODO: the code does the sort and filter for .bed files, it should be modified to call the function as well

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

classList = segwayLabels # I want them to be sorted the same way, so I will use the hard coded listed
classCount = len(classList) # this is the default for our Segway annotation annotations - 9
classListIndMap = {}
for i in range(len(classList)):
   classListIndMap[classList[i]] = i
   
extension = 3000 # count of basepairs monitored before and after the gene coordinates

annAccessionCount = 0
for annAccession in annAccessionList:
    
    # for each annotation
    if (annMeta[annAccession]['expressionFile'] != 'none'):
        
        #  annAccession = annAccessionList[5] #>>>> test
        #  print(annMeta[annAccession]['expressionFile'] != 'none') #>>>> test
        
        print(annAccession)
        annAccessionCount += 1
        print(annAccessionCount)
        
        expFileName = annMeta[annAccession]['expressionFile']
        expAccession = re.split('_|\.', expFileName)[2]
        inputFile = dataFolder + dataSubFolder + annAccession + '/geneExp_dict_' + expAccession + '.pkl'
        print(inputFile)
        with open(inputFile, 'rb') as pickledFile:
            expression = pickle.load(pickledFile)


        # using the sorted and filtered .bed file instead of the initial one
        annFile = dataFolder + dataSubFolder + annAccession + '/' + annMeta[annAccession]['bedFile'].split('.')[0] + '_filteredSorted.bed' 

        bedFileName = annMeta[annAccession]['bedFile'].split('.')[0]
        inputFile = dataFolder + dataSubFolder + annAccession + '/' + bedFileName + '_annotationSummary.pkl'
        with open(inputFile, 'rb') as pickledFile:
            summaryAnnotation = pickle.load(pickledFile)

        annLineCount = int(sp.getoutput('wc -l %s' %(annFile)).split()[0])

        #print(summaryAnnotation['classes']) # >>>> test
        #print(summaryAnnotation['clusters']) # >>>> test

        classes = summaryAnnotation['classes']
        clusters = summaryAnnotation['clusters']
        
        # Sorting the list of clusters TODO: this part should ideally go to the annotation processing and class creation. but I think it is better to do it here regardless just to be sure
        clusterList = list(clusters.keys())
        sortedClusterList = []
        #print(segwayLabels) # >>>> test
        for label in segwayLabels:
            #print(label)
            for item in clusterList:
                #print(item) # >>>> test
                if item.split('_')[1] == label:
                    sortedClusterList.append(item)

                    
        # getting the cluster list index map
        clusterList = list(clusters.keys()) # this is per sample
        clusterCount = len(clusterList) # note: this could vary for each annotation file
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
        

        annLineInd = 1 # this keeps the annotation line index, the first line is empty, thus the index starts at 1
        genomic_region_start_annLineInd = 0


        previous_gene_chr = 'chr1'
        previous_extension_end = 0 # doesn't matter since chr condition is never true until the first gene is processed
        previous_ann_chr = 'chr1'

        #while cgi < len(geneIDList) and annLineInd < annLineCount: # >>>> test: modify the condition for the test runs
        while cgi < 26017: # >>>>>>>>>> MAIN
        #while cgi < 18431:# >>>>
            
            print('cgi')
            print(cgi)

            'we are in the gene territory, we are walking on the genes'
   
            geneID = geneIDList[cgi]
            gene_chr = geneList[geneIDList[cgi]].chrom
            #print(geneID) # >>>> test
            
            gene_strand = geneList[geneID].strand
            gene_start = geneList[geneID].start
            gene_end = geneList[geneID].end
            gene_length = gene_end - gene_start
            gene_length_unit = int(gene_length/100)
            if gene_length_unit == 0:
                #print('gene smaller than 100bp, skipping it') # >>>> test
                cgi = cgi + 1
                continue

            gene_length_last_unit = gene_length - (99* gene_length_unit)
            # TODO: ideally: something to fix: the gene_length_last_unit for negative strand versus positive strand

            extension_start = gene_start - extension
            extension_end = gene_end + extension

            #print('%ss, %s, %s' %(gene_chr, extension_start, extension_end)) # >>>> test

            
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

                
            if (ann_chr != gene_chr): # if annotation moved to the next chromosome, but gene has not yet moved to the next chromosome
                annLineInd = annLineInd - 2
                line = linecache.getline(annFile, annLineInd)
                annLineInd +=1
                ann_line_count += 1
                fields = line.strip().split()

                ann_chr = fields[0]
                ann_start = int(fields[1])
                ann_end = int(fields[2])


            ''' 
            if this gene starts somewhere in the preivous genomic region, 
            I will go back in the annotation file to the beginning of annotation for the previous gene
            changed this th the next while loop - basically we will go back until the annotations start before the gene territory

            '''

            while (ann_start > extension_start) and (gene_chr == ann_chr): # in case of overlapping genes
                #print('ann start greater than extension start, getting back in annotation until it is not') # >>>> test
                #print(annLineInd) # >>>> test
                #print('%ss, %s, %s' %(ann_chr, ann_start, ann_end)) # >>>> test
                
                annLineInd = max(annLineInd - 5,1)
                line = linecache.getline(annFile, annLineInd)
                annLineInd +=1
                ann_line_count += 1
                fields = line.strip().split()

                ann_chr = fields[0]
                ann_start = int(fields[1])
                ann_end = int(fields[2])
                
                #print('overlapping genes here') # >>>> test
                #print(annLineInd) # >>>> test

            #while ((ann_start < extension_end) or not(gene_chr == ann_chr)) and geneMatWalkIndex < 160: 
            while ((ann_start < extension_end) and (gene_chr == ann_chr)) and geneMatWalkIndex < 160:
                
#            if (ann_chr != previous_ann_chr): # in case of chromosome change because annotation moved to the next chromosome
 #           if (ann_chr != gene_chr): # if annotation moved to the next chromosome, but gene has not yet moved to the next chromosome

                '''
                NOTE: the second condition is for when we are at the end of a gene's territory, and then at the end of a chromosome, so when we go to the next gene. At some point I changed the second condition from "or" to "and" 
                ann_start is not smaller than extension_end, and we need to read the annotations until we are on the same chromosome
                The condition to be in one gene's territory (and not the next one), and it is for reading the annotations
                while in THIS gene territory, read annotations until we reach to the next genomic region (plus extension)
      
                '''

                #print('in the gene region, reading annotations until we are out') # >>>> test
                #print('%ss, %s, %s' %(gene_chr, extension_start, extension_end)) # >>>> test

                # in the next sections we are processing the annotation
                if ann_chr == gene_chr: # if we are in the same choromosome
                    if ((ann_start < extension_end and ann_start > extension_start) or 
                        (ann_end < extension_end and ann_end > extension_start) or
                        (ann_start < extension_start and ann_end > extension_end)):


                        ''' We are in the genonimc region (with extension)'''
                        #print('annotation') # >>>> test
                        #print('%ss, %s, %s' %(ann_chr, ann_start, ann_end)) # >>>> test


                        ''' Taking off for filling the matrices... '''
                        
                        adjusted_ann_start = max(0, ann_start - extension_start)
                        adjusted_ann_end = min(ann_end - extension_start, extension_end - extension_start)
                        adjusted_ann_length = adjusted_ann_end - adjusted_ann_start

                        ann_cluster = fields[3]
                        clusterInd = clusterListIndMap[ann_cluster]

                        ann_class = ann_cluster.split('_')[1]
                        classInd = classListIndMap[ann_class]

                        # expMatInd

                        if gene_strand == '+':  # munching 100bp s from annotation, filling the geneMat
                            while(adjusted_ann_length > 0) and geneMatWalkIndex < 30:
                                if (adjusted_ann_length >= 100*(1- previous_fill)):
                                    classMats[expMatInd][classInd][geneMatWalkIndex] += 1 - previous_fill
                                    clusterMats[expMatInd][clusterInd][geneMatWalkIndex] += 1 - previous_fill
                                    adjusted_ann_length -= 100*(1- previous_fill)
                                    previous_fill = 0
                                    geneMatWalkIndex +=1
                                else:
                                    classMats[expMatInd][classInd][geneMatWalkIndex] += adjusted_ann_length/100
                                    clusterMats[expMatInd][clusterInd][geneMatWalkIndex] += adjusted_ann_length/100
                                    previous_fill += adjusted_ann_length/100
                                    adjusted_ann_length = 0
                                    if previous_fill > 1:
                                        previous_fill = 1

                            while(adjusted_ann_length > 0) and geneMatWalkIndex < 129:
                                if (adjusted_ann_length >= gene_length_unit*(1- previous_fill)):
                                    classMats[expMatInd][classInd][geneMatWalkIndex] += 1 - previous_fill
                                    clusterMats[expMatInd][clusterInd][geneMatWalkIndex] += 1 - previous_fill
                                    adjusted_ann_length -= gene_length_unit*(1 - previous_fill)
                                    previous_fill = 0
                                    geneMatWalkIndex +=1
                                else:
                                    classMats[expMatInd][classInd][geneMatWalkIndex] += adjusted_ann_length/gene_length_unit
                                    clusterMats[expMatInd][clusterInd][geneMatWalkIndex] += adjusted_ann_length/gene_length_unit
                                    previous_fill += adjusted_ann_length/gene_length_unit
                                    adjusted_ann_length = 0
                                    if previous_fill > 1:
                                        previous_fill = 1

                            if (adjusted_ann_length > 0) and geneMatWalkIndex == 129:
                                if (adjusted_ann_length >= gene_length_last_unit*(1- previous_fill)):
                                    classMats[expMatInd][classInd][geneMatWalkIndex] += 1 - previous_fill
                                    clusterMats[expMatInd][clusterInd][geneMatWalkIndex] += 1 - previous_fill
                                    adjusted_ann_length -= gene_length_last_unit*(1 - previous_fill)
                                    previous_fill = 0
                                    geneMatWalkIndex +=1
                                else:
                                    classMats[expMatInd][classInd][geneMatWalkIndex] += adjusted_ann_length/gene_length_last_unit
                                    clusterMats[expMatInd][clusterInd][geneMatWalkIndex] += adjusted_ann_length/gene_length_last_unit
                                    previous_fill += adjusted_ann_length/gene_length_last_unit
                                    adjusted_ann_length = 0
                                    if previous_fill > 1:
                                        previous_fill = 1

                            while(adjusted_ann_length > 0) and (geneMatWalkIndex < 160):
                                if (adjusted_ann_length >= 100*(1- previous_fill)):
                                    classMats[expMatInd][classInd][geneMatWalkIndex] += 1 - previous_fill
                                    clusterMats[expMatInd][clusterInd][geneMatWalkIndex] += 1 - previous_fill
                                    adjusted_ann_length -= 100*(1- previous_fill)
                                    previous_fill = 0
                                    geneMatWalkIndex +=1
                                else:
                                    classMats[expMatInd][classInd][geneMatWalkIndex] += adjusted_ann_length/100
                                    clusterMats[expMatInd][clusterInd][geneMatWalkIndex] += adjusted_ann_length/100
                                    previous_fill += adjusted_ann_length/100
                                    adjusted_ann_length = 0
                                    if previous_fill > 1:
                                        previous_fill = 1


                        if gene_strand == '-':
                           # print('here in the negative strand') # >>>> test
                            while(adjusted_ann_length > 0) and geneMatWalkIndex < 30:
                            #    print('first while') # >>>> test
                                if (adjusted_ann_length >= 100*(1- previous_fill)):
                                    classMats[expMatInd][classInd][159 - geneMatWalkIndex] += 1 - previous_fill
                                    clusterMats[expMatInd][clusterInd][159 - geneMatWalkIndex] += 1 - previous_fill
                                    adjusted_ann_length -= 100*(1- previous_fill)
                                    previous_fill = 0
                                    geneMatWalkIndex +=1
                                else:
                                    classMats[expMatInd][classInd][159 - geneMatWalkIndex] += adjusted_ann_length/100
                                    clusterMats[expMatInd][clusterInd][159 - geneMatWalkIndex] += adjusted_ann_length/100
                                    previous_fill += adjusted_ann_length/100
                                    adjusted_ann_length = 0
                                    if previous_fill > 1:
                                        previous_fill = 1

                            if (adjusted_ann_length > 0) and geneMatWalkIndex == 30:
                                if (adjusted_ann_length >= gene_length_last_unit*(1- previous_fill)):
                                    classMats[expMatInd][classInd][159 - geneMatWalkIndex] += 1 - previous_fill
                                    clusterMats[expMatInd][clusterInd][159 - geneMatWalkIndex] += 1 - previous_fill
                                    adjusted_ann_length -= gene_length_last_unit*(1 - previous_fill)
                                    previous_fill = 0
                                    geneMatWalkIndex +=1
                                else:
                                    classMats[expMatInd][classInd][159 - geneMatWalkIndex] += adjusted_ann_length/gene_length_last_unit
                                    clusterMats[expMatInd][clusterInd][159 - geneMatWalkIndex] += adjusted_ann_length/gene_length_last_unit
                                    previous_fill += adjusted_ann_length/gene_length_last_unit
                                    adjusted_ann_length = 0
                                    if previous_fill > 1:
                                        previous_fill = 1


                            while(adjusted_ann_length > 0) and geneMatWalkIndex < 129:
                                if (adjusted_ann_length >= gene_length_unit*(1- previous_fill)):
                                    classMats[expMatInd][classInd][159 - geneMatWalkIndex] += 1 - previous_fill
                                    clusterMats[expMatInd][clusterInd][159 - geneMatWalkIndex] += 1 - previous_fill
                                    adjusted_ann_length -= gene_length_unit*(1 - previous_fill)
                                    previous_fill = 0
                                    geneMatWalkIndex +=1
                                else:
                                    classMats[expMatInd][classInd][159 - geneMatWalkIndex] += adjusted_ann_length/gene_length_unit
                                    clusterMats[expMatInd][clusterInd][159 - geneMatWalkIndex] += adjusted_ann_length/gene_length_unit
                                    previous_fill += adjusted_ann_length/gene_length_unit
                                    adjusted_ann_length = 0
                                    if previous_fill > 1:
                                        previous_fill = 1


                            while(adjusted_ann_length > 0) and (geneMatWalkIndex < 160):
                                if (adjusted_ann_length >= 100*(1- previous_fill)):
                                    classMats[expMatInd][classInd][159 - geneMatWalkIndex] += 1 - previous_fill
                                    clusterMats[expMatInd][clusterInd][159 - geneMatWalkIndex] += 1 - previous_fill
                                    adjusted_ann_length -= 100*(1- previous_fill)
                                    previous_fill = 0
                                    geneMatWalkIndex +=1
                                else:
                                    classMats[expMatInd][classInd][159 - geneMatWalkIndex] += adjusted_ann_length/100
                                    clusterMats[expMatInd][clusterInd][159 - geneMatWalkIndex] += adjusted_ann_length/100
                                    previous_fill += adjusted_ann_length/100
                                    adjusted_ann_length = 0
                                    if previous_fill > 1:
                                        previous_fill = 1


                                # if gene strand is positive, or negative, flag which cluster list we use
                                # the only difference between the two strands is that I am going to reverse the index of

                if geneMatWalkIndex < 160: # we read annotations until we cover the genee
                    #print('geneMatWalkInd') # >>>> test
                    #print(geneMatWalkIndex) # >>>> test
                    #print('annLineInd') # >>>> test
                    #print(annLineInd) # >>>> test
                    line = linecache.getline(annFile, annLineInd)
                    annLineInd +=1
                    ann_line_count += 1
                    fields = line.strip().split()

                    ann_chr = fields[0]
                    ann_start = int(fields[1])
                    ann_end = int(fields[2])
                    #print('%ss, %s, %s' %(ann_chr, ann_start, ann_end)) # >>>> test

                if geneMatWalkIndex >= 159:
                    print('gene Done')
                    print(geneMatWalkIndex)


            cgi += 1 # next gene
            previous_extension_end = extension_end
            previous_gene_chr = gene_chr
            

        linecache.clearcache()

        expOverlapMats = {"clusterMats": clusterMats, "classMats": classMats}
        outputFile = dataFolder + dataSubFolder + annAccession + '/' + 'defaultExp_5kg_expSummary.pkl'
        with open(outputFile, 'wb') as f:
            pickle.dump(expOverlapMats, f)

        book = pd.DataFrame(clusterMats[2])
        sns.heatmap(book)
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


        print('annotation summary saved in %s' %(outputFile))

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

########################################
# 4. Preprocessing transcriptomic data
########################################

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'

# load the meta file 
inputFileName = 'all235Annot_meta_corrected.pkl'
inputFile = dataFolder +  outputFileName
with open(inputFile, 'rb') as f:
    allMeta = pickle.load(f)

accessionList = list(allMeta.keys())

count = 0
for accession in accessionList:
    annotation = allMeta[accession]

    if ((('38batch' in annotation['folder']) or ('May11' in annotation['folder'])) and not(annotation['RNAseqFile'] == 'none')):
        expFile = annotation['RNAseqFile'][0]
        count+=1

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
        sampleFolderAdd = allMeta[accession]['folder']
        expAccession = re.split('_|\.', expFile)[-2]
        outputFile = sampleFolderAdd + 'geneExp_dict_' + expAccession + '.pkl'
        print('printing this file %s' %(outputFile))
        with open(outputFile, 'wb') as f:
            pickle.dump(expression, f)


########################################
# 5. process the annotation file
########################################

inputFileName = 'all235Annot_meta_corrected.pkl'
inputFile = dataFolder +  outputFileName
with open(inputFile, 'rb') as f:
    allMeta = pickle.load(f)

accessionList = list(allMeta.keys())

for accession in accessionList:
    annotation = allMeta[accession]

    if ('38batch' in annotation['folder']):
        print('####### 38')
        sampleFolderAdd = annotation['folder']
        print(sampleFolderAdd)

        # get the bed file
        originalBedFile = sampleFolderAdd + 'segOutput/' + 'segway.bed.gz'
        print(originalBedFile)
         
        # unzip the bed file
        os.system('gunzip %s' %(originalBedFile))
        originalBedFile = originalBedFile[0:-3]
        print(originalBedFile)

        # do the modifications
        # sort the bed file
        sortedBedFile = sampleFolderAdd + 'sortedBedFile.bed'
        os.system('sort -V -k1,1 -k2,2 %s > %s' %(originalBedFile, sortedBedFile))

        # filter the bed file 
        filteredSortedBedFile = sampleFolderAdd + 'bedFile_filteredSorted.bed'
        os.system("awk -F'\t' '$1 !~ \"_\"' %s > %s" %(sortedBedFile, filteredSortedBedFile))

        tempFile = sampleFolderAdd + 'temp.bed'
        os.system('grep -v "chrM" %s > %s && mv %s %s' %(filteredSortedBedFile, tempFile, tempFile, filteredSortedBedFile))

        os.remove(sortedBedFile)
    
        # zip the .bed file
        os.system('gzip %s' %(originalBedFile))

        # zip the modified version
        os.system('gzip %s' %(filteredSortedBedFile))

    if ('May11' in annotation['folder']):
        print('######')
        sampleFolderAdd = annotation['folder'][0:-1]
        print(sampleFolderAdd)

        # get the bed file
        originalBedFile = sampleFolderAdd + 'call-segway/' + 'segway.bed.gz'
        print(originalBedFile)
         
        # unzip the bed file
        os.system('gunzip %s' %(originalBedFile))
        originalBedFile = originalBedFile[0:-3]
        print(originalBedFile)

        # do the modifications
        # sort the bed file
        sortedBedFile = sampleFolderAdd + 'sortedBedFile.bed'
        os.system('sort -V -k1,1 -k2,2 %s > %s' %(originalBedFile, sortedBedFile))

        # filter the bed file 
        filteredSortedBedFile = sampleFolderAdd + 'bedFile_filteredSorted.bed'
        os.system("awk -F'\t' '$1 !~ \"_\"' %s > %s" %(sortedBedFile, filteredSortedBedFile))

        tempFile = sampleFolderAdd + 'temp.bed'
        os.system('grep -v "chrM" %s > %s && mv %s %s' %(filteredSortedBedFile, tempFile, tempFile, filteredSortedBedFile))

        os.remove(sortedBedFile)
    
        # zip the .bed file
        os.system('gzip %s' %(originalBedFile))

        # zip the modified version
        os.system('gzip %s' %(filteredSortedBedFile))


########################################
# 6. use the function to get the transcription enrichment for Segway
########################################

from transcription_overlap import SegwayTranscriptionEnrichment
from transcription_overlap import SegwayGeneBodyEnrichment

# I need to do it for the 21 remaining samples with the transcriptomic data
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'

# GTF data structure
fileName = dataFolder + '/geneLists.pkl'
with open(fileName, 'rb') as pickleFile:
    geneListsAndIDs = pickle.load(pickleFile)
    
geneList = geneListsAndIDs[0]
geneIDList = geneListsAndIDs[1]
del geneListsAndIDs

# load the meta file 
inputFileName = 'all235Annot_meta_corrected.pkl'
#inputFileName = 'all235Annot_meta.pkl'
inputFile = dataFolder +  inputFileName
with open(inputFile, 'rb') as f:
    allMeta = pickle.load(f)

accessionList = list(allMeta.keys())

count = 0
for accession in accessionList:
    annotation = allMeta[accession]
    count+=1

    #if ((('38batch' in annotation['folder']) or ('May11' in annotation['folder'])) and not(annotation['RNAseqFile'] == 'none')):
    if(not(annotation['RNAseqFile'] == 'none')):
        print(count)
        if len(annotation['RNAseqFile']) > 5:
            RNAFile = annotation['RNAseqFile']
        else:
            RNAFile = annotation['RNAseqFile'][0]
        print(count)
        print(accession)
        print(RNAFile)

        annotationFolder = annotation['folder']
        print(annotationFolder)

        expAccession = RNAFile[-15:-4]
        RNAseqFile = annotationFolder +  '/geneExp_dict_' + expAccession + '.pkl'
        annFile = annotation['bedFile']
        mnemFile = annotationFolder + 'mnemonics_v04.txt'
        extension = 3000
        SegwayTranscriptionEnrichment(annotationFolder, annFile, RNAseqFile, extension, geneList, geneIDList, mnemFile)


# calling the genebody enrichment function

count = 0
for accession in accessionList[206:]:
    annotation = allMeta[accession]

    #if ((('38batch' in annotation['folder']) or ('May11' in annotation['folder'])) and not(annotation['RNAseqFile'] == 'none')):
    if((annotation['RNAseqFile'] == 'none')):
        count+=1
        print(count)
        print(accession)

        annotationFolder = annotation['folder']
        print(annotationFolder)

        annFile = annotation['bedFile']
        mnemFile = annotationFolder + 'mnemonics_v04.txt'
        #geneList
        #geneIDList =
        extension = 3000
        SegwayGeneBodyEnrichment(annotationFolder, annFile, extension, geneList, geneIDList, mnemFile)


        
         
