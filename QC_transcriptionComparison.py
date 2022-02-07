# Here I do checks and QC with relevance to the transcriptomic data
##
#For this task we need three input files:
#
#- This is the exact GTF: https://www.encodeproject.org/files/gencode.v29.primary_assembly.annotation_UCSC_names/ genomic coordinates used for Segtools
#- The RNA-seq data for the sample(s)
#- The annotations with class labels (the initial sample with the transcript data name:H1 ID:ENCSR938GXK)
## 0. initials
## 0.1 checking the fields for the GTF File/reading the file
## 0.2 preprocess and QC for the RNA-seq file
## 0.3 exploratory analysis for the annotation file
## 1. Parsing the annotation file, considering genomic coordinates
## 2. Getting the gene structures from GTF
## 3. Comparing the annotation labels to the genomic regions and transcriptomic data

###############################
# 0. initials
###############################

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle

# data folder
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'

# the GTF file
gftFile = dataFolder + 'testBatch/gencode.v29.primary_assembly.annotation_UCSC_names.gtf'

# sample list (by folder name)
sample = 'H1'

segwayLabels = ['Quiescent', 'ConstitutiveHet', 'FacultativeHet', 'Transcribed', 'Promoter', 'Enhancer', 'RegPermissive', 'Bivalent', 'LowConfidence']

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 0.1 checking the fields for the GTF File/reading the file
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>
# the GTF has multiple count of fields for different item identities in column 3, and 2,742,734 coordinate records.

#>>>>> some exploration for the files, columns and unique types etc.
with open(gftFile, 'r') as coors:
   line = coors.readline()
   while line[0] == '#': # skipping header lines
        print(line)
        line = coors.readline()
        
    # fetching uniqe gene_types from the file:
#    record = line.strip().split()
#    print(record)
#    record = coors.readline().strip().split()
#    print(record)
    
   geneTypeList = []
   geneTypeLineCount = [] # line count per 'gene_type'
   geneTypeGeneCount = [] # count of genes per 'gene_type' (we are not counting all lines)
   for line in coors: # reading
      index = line.index('gene_type')
      thisType = line[index:].split('"')[1]
      if not(thisType in geneTypeList):
         geneTypeList.append(thisType)
         geneTypeLineCount.append(0)
         geneTypeIndex = geneTypeList.index(thisType)
         geneTypeGeneCount.append(0)
      else:
         geneTypeIndex = geneTypeList.index(thisType)
         geneTypeLineCount[geneTypeIndex] += 1
         
         features = line.strip().split()
         if features[2] == 'gene':
            geneTypeGeneCount[geneTypeIndex] += 1

geneTypesFromGTF = [geneTypeList, geneTypeLineCount, geneTypeGeneCount]
for i in range(len(geneTypeList)):
    print('%s, count:%d' %(geneTypeList[i], geneTypeCount[i]))

for i in range(len(geneTypeList)):
    print('%s, count:%d' %(geneTypeList[i], geneTypeGeneCount[i]))

geneTypesFromGTF = [geneTypeList, geneTypeLineCount, geneTypeGeneCount]

fileName = dataFolder + '/geneTypesFromGTF.pkl'
with open(fileName, 'wb') as f:
    pickle.dump(geneTypesFromGTF, f)

# loading classes and labels
fileName = dataFolder + '/geneTypesFromGTF.pkl'
with open(fileName, 'rb') as pickledFile:
    geneTypesFromGTF = pickle.load(pickledFile)

geneTypeList = geneTypesFromGTF[0]    
    
#>>>>> getting the gene coordinates: based on the third column, I have 58,780 genes in the file
# I want the coordinates, gene name, gene ID, the strand, the type (both coding and non_coding)

chrom = []
start = []
end = []
strand = []
name = []
geneID = []
geneType = []
 
with open(gftFile, 'r') as coors:
    line = coors.readline()
    while line[0] == '#': # skipping header lines
        print(line)
        line = coors.readline()

    fields = line.strip().split()
    chrom.append(fields[0])
    start.append(fields[3])
    end.append(fields[4])
    strand.append(fields[6])
    name.append(fields[13])
    geneID.append(fields[9])
    geneType.append(fields[11])
    
    for line in coors:
        fields = line.strip().split()
        if fields[2] == 'gene':
            chrom.append(fields[0])
            start.append(fields[3])
            end.append(fields[4])
            strand.append(fields[6])
            name.append(fields[13])
            geneID.append(fields[9])
            geneType.append(fields[11])

    coors_frame = pd.DataFrame({'chrom': chrom,
                                'start': start,
                                'end': end,
                                'strand': strand,
                                'name': name,
                                'geneID': geneID,
                                'geneType': geneType})

#>>>>> save the gene coordinates
fileName = dataFolder + 'testBatch/coorsFile.pkl'
coors_frame.to_pickle(fileName)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 0.2 preprocess and QC for the RNA-seq file
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# I am doing this per sample, for the initial test the sample is H1
# from the RNAseq file, I want the geneID, transcript ID and then the TPM: 0,1,5

# reading the expression values into a dict
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

RNAseqFile = dataFolder + 'testBatch/' + sample + '/ENCFF432RPO.tsv'

expression = {}
with open(RNAseqFile, 'r') as expFile:
    line = expFile.readline() # there is one header line

    for line in expFile:
        fields = line.strip().split()
        geneID = fields[0]
        transcriptID = fields[1]

        if geneID in expression:
            expression[geneID] += np.log10(float(fields[5]) + 1)
        else:
            expression[geneID] = np.log10(float(fields[5]) + 1)
        

#TODO : save expression dictionary

# reading the expression values in a dataframe
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

RNAseqFile = dataFolder + 'testBatch/' + sample + '/ENCFF432RPO.tsv'

geneID = []
transcriptID = []
expression = []
with open(RNAseqFile, 'r') as expFile:
    line = expFile.readline() # there is one header line

    for line in expFile:
        fields = line.strip().split()
        geneID.append(fields[0])
        transcriptID.append(fields[1])
        expression.append(float(fields[5]))

    exp_frame = pd.DataFrame({'geneID': geneID,
                              'transcriptID': transcriptID,
                              'expression': expression})

fileName = dataFolder + 'testBatch/' + sample + '/expFile.pkl'
exp_frame.to_pickle(fileName)

#>>>>> just plotting the expression levels
print(sum(expression))
print(max(expression))
highExp = [x for x in expression if x > 2]
print(len(highExp)) # 12,166 genes expressed, as expected

arrayExp = np.array(expression)
logexp = numpy.log10(arrayExp + 1)

inds = np.argwhere(logexp > .5)

plt.hist(logexp[inds], bins = 100)
plt.xlabel('log10 TPM expression > 0.5')
plt.ylabel('gene counts')
plt.title('gene expression distribution')
figFile = dataFolder + 'testBatch/' + sample + '/exp_2TPMfilter_log10.pdf'
plt.savefig(figFile)
plt.show()

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 0.3 exploratory analysis for the annotation file
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#TODO: go through the file and count bps of different regions

########################################
# 1. Parsing the annotation file, considering genomic coordinates
########################################
# 1. Make the smaller .bed file for genomic regions
# 2. Get the class/cluster/annotation info from the .bed file

'''
Inputs for the section: 
1. pcgenes / genomic coordinates from the coordinate file
2. segway.bed file
'''

coors_frame
exp_frame
annFile = dataFolder + 'testBatch/' + sample + '/ccSegway.bed' 
annFileGR = dataFolder + 'testBatch/' + sample + '/ccSegway_genomicRegions.bed'

# get all the protein coding genes from the coors
pcgenes = coors_frame[coors_frame['geneType'].str.contains('protein_coding')]

# >>>>>>>>>>>>>>>>>>>>>>>>> DRAFT: I can either have a list of named tuples like this and then keep the index somewhere else
# unique combination of cluster + classes. I am keeping their info here
label = namedtuple('label', ['called', # the combination of the class + cluster
                             'category', # class label
                             'color', # color
                             'bp_count', # total count of basepairs
                             'region_count', # total count of regions
                             'region_dist']) # length distribution
ann_labels = [] # This is a list of labels.
# >>> or have them as dictionaries, with key the name and value the info tuple
info = namedtuple('info', ['called', # the combination of the class + cluster
                             'category', # class label
                             'color', # color
                             'bp_count', # total count of basepairs
                             'region_count', # total count of regions
                             'region_dist']) # length distribution # we have it from segtools, so I am keeping it empty for now
# name tuples won't work cause they are immuatable
# <<<<<<<<<<<<<<<<<<<<<<<<<< DRAFT


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Annotation and classes (biolabels)

# keeping the unique cluster_class labels and their info for annotations
# I defined my own classes since namedtuples are immutable - I preferred namedtuples to dict since they are defined with fields - 

class Annotation(object):
    def __init__(self, called, biolabel, cluster, color, bp_count, region_count, region_dist):
        self.called = called
        self.biolabel = biolabel
        self.cluster = cluster
        self.color = color
        self.bp_count = bp_count
        self.region_count = region_count
        self.region_dist = region_dist # not used currently

    def __str__(self):
        return 'called = %s, biolabel = %s, cluster = %s, color = %s, bp_count = %d, region_count = %d, region_dist = none' %(self.called, self.biolabel, self.cluster, self.color, self.bp_count, self.region_count)

class AnnotationClass(object):
    def __init__(self, biolabel, clusters, color, bp_count, region_count, region_dist):
        self.biolabel = biolabel
        self.clusters = clusters
        self.color = color
        self.bp_count = bp_count
        self.region_count = region_count # the merged count
        self.region_dist = region_dist # not used currently

    def __str__(self):
        return 'biolabel = %s, clusters = %s, color = %s, bp_count = %d, region_count = %d, region_dist = none' %(self.biolabel, self.clusters, self.color, self.bp_count, self.region_count)

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Annotation and classes (biolabels)


labels = {} # list of annotation labels
classes = {} # list of annotationClass
extension = 1500 # count of basepairs monitored before and after the gene coordinates

cgi = 0 # walks on the genomic region
ann_start = 0
ann_end = 0
ann_line_count = 0
previous_class = ''

with open(annFile, 'r') as annotations, open(annFileGR, 'w') as grFile:
    
    '''
    Walking on the genomic coordinates in the protein_coding dataframe <pcgenes> and the annotation file doing:
    1. Filling up the annotation info for basepair counts
    2. Recording the annotations for the extended genomic regions in another file (for later investigation in genomic regions)

    '''

    #while cgi < len(pcgenes): # modify the condition for the test runs
    while cgi < len(pcgenes): # modify the condition for the test runs
        
        gene_chr = pcgenes.iloc[cgi].chrom

        gene_start = int(pcgenes.iloc[cgi].start) - extension
        gene_end = int(pcgenes.iloc[cgi].end) + extension


        # for line in annotations:
        while (ann_start < gene_end) or not(gene_chr == ann_chr):
            ann_line_count += 1
            
            line = annotations.readline()
            fields = line.strip().split()
            
            ann_chr = fields[0]
            ann_start = int(fields[1])
            ann_end = int(fields[2])
            
            # fill up the annotation info: if label is not in the file,
            if fields[3] in labels.keys():
                
                # print('already added')
                labels[fields[3]].bp_count += int(fields[2]) - int(fields[1])
                labels[fields[3]].region_count += 1
                # (TODO - maybe?) add to distribution
                
            else:
                
                # print('not added')
                called = fields[3]
                biolabel = fields[3].split('_')[1]
                cluster = fields[3].split('_')[0]
                color = fields[8]
                bp_count = int(fields[2]) - int(fields[1])
                region_count = 1
                region_dist = 1
                labels[fields[3]] = Annotation(called, biolabel, cluster, color, bp_count, region_count, region_dist)

                '''
                # the printing block 
                allAnns = list(labels.keys())
                for ann in allAnns:
                print(str(labels[ann]))
                '''

            if previous_class == fields[3].split('_')[1]:
                classes[previous_class].bp_count +=  int(fields[2]) - int(fields[1])
            else:
                current_class = fields[3].split('_')[1]
                if current_class in classes.keys():
                    classes[current_class].bp_count +=  int(fields[2]) - int(fields[1])
                    classes[current_class].region_count += 1
                else:
                    clusters = [] # not filling it now, it can be filled later using annotations 
                    biolabel = current_class
                    color = fields[8]
                    bp_count = int(fields[2]) - int(fields[1])
                    region_count = 1
                    region_dist = 1
                    classes[biolabel] = AnnotationClass(biolabel, clusters, color, bp_count, region_count, region_dist)

            
            previous_class = fields[3].split('_')[1]
            # write the genomic region annotation
            if ann_chr == gene_chr:
                if (ann_start < gene_end and ann_start > gene_start) or (ann_end < gene_end and ann_end > gene_start) or (ann_start < gene_start and ann_end > gene_end):
                    grFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n'
                                 %(fields[0],fields[1],fields[2],fields[3],fields[4],fields[5],fields[6],fields[7],fields[8]))

        ###
        #TODO: we are probably missing something at the end of the last gene
        ###
        '''
        # just checking
        print('cgi = %d' %(cgi))
        print(pcgenes.iloc[cgi])
        print(fields)
        '''
        cgi += 1 # next gene

        
#post generation modification and saving
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# sorting classes and labels
##### TODO


# saving classes and labels
summaryAnnotation = {"classes" : classes, "labels": labels}

fileName = dataFolder + 'testBatch/' + sample + '/summaryAnnotation.pkl'
with open(fileName, 'wb') as f:
    pickle.dump(summaryAnnotation, f)

# loading classes and labels
fileName = dataFolder + 'testBatch/' + sample + '/summaryAnnotation.pkl'
with open(fileName, 'rb') as pickledFile:
    summaryAnnotation = pickle.load(pickledFile)

# some tests reading annotation file
#>>>>>>>>>>>>>>>>>>>>
with open(annFile, 'r') as f:
   for i in range(3):
      print(f.readline())
   print(f.tell())

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PLOTS for class and annotation labels
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# plot: mean length for both annotations and classes
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

classList = list(classes.keys())
classMeanLength = np.zeros(len(classList))
for i, key in enumerate(classList):
    classMeanLength[i] = int(classes[key].bp_count)/int(classes[key].region_count)

annotationList = list(labels.keys())
annotMeanLength = np.zeros(len(annotationList))
for i, key in enumerate(annotationList):
    annotMeanLength[i] = int(labels[key].bp_count)/int(labels[key].region_count)

# TODO: the plot
    
# plot: bp/all ratio for both annotations and classes
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

totalbp = 0
classList = list(classes.keys())
classbp = np.zeros(len(classList)) 
for i, key in enumerate(classList):
    classbp[i] = int(classes[key].bp_count)
    totalbp += classes[key].bp_count 

annotationList = list(labels.keys())
annotbp = np.zeros(len(annotationList))
for i, key in enumerate(annotationList):
    annotbp[i] = int(labels[key].bp_count)
    
# TODO: the plot

# plot: mean value of the histone tracks
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TODO

'''from the histone track info file'''

# plot: length distribution 
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# TODO 

'''from the length distribution file'''

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# PLOTS for class and annotation labels
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# using linecache for reading annotation file
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
import linecache
for i in range(10):
   line = linecache.getline(annFile, i)
   print(line)

with open(annFile, 'r') as annotations:
   line = linecache.getline(annotations, 1)

########################################
# 2. Getting the gene structures from GTF
########################################
# From the .bed files:
# 1. Enrichment of transcript bps in the genomic regions, and in the expressed regions: go through the file once and get this
# 2. Enrichment of promoter regions in 1000bps before the genomic regions, 1000bps before the expressed genes, and in the whole files
# (anything to do with the enhancers?)
# 3. Exploratory: for the 200bps start site and 200bp end site, enrichment of labels

'''
TODO
1. Get the gene structure from the GTF file, where I have coordinates for exons and for genes. A gene must have list of exons. 
Exons have begin, end and ID (first, second, etc) - also the intron that comes in- between.
2. Get the genomic annotations for 3k around the genes. 
3. Do the transcript and expression and all the else

Data/datastructure: plan for gene structure: object gene, object exon, I will also have the gene coordinate dataframe (sorted) and the annotation file for genomic regions (sorted)

NOTE: while filling the info, the first exon is the one with smallest start, exons should be sorted based on their start 

I will then go through the coordinates dataframe and fetch genes and exons info - from here I will know about the count of exons for each gene, first or last exon etc - exons are distinguished with their IDs

getting the enrichment: I will go through the gene dataframe and genomic annotation file using the loop in 1. But here, at each gene, I will also walk through the exons and getting label enrichment for them. 

label enrichment: we count the basepairs in that region with a specific label. So far, we do it with the "LABELS" and not "CLASSES"/BIOLABELS. This is because I highly suspect that some classes are something else. This is also a good investigation step

We do this for all the genes. 

Then we do this for expressed genes, or a specific group of genes.

'''

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> I examined gene type records in GTF
gene_type = 'protein_coding'
gene_type = 'snoRNA'
file = open(gftFile, 'r')

line = file.readline()
while not 'AP006222.1' in line:
    line = file.readline()

print(line)
line = file.readline()

while gene_type not in line:
    line = file.readline()

while gene_type in line:
    line = file.readline()

print(line)
line = file.readline()    

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Genes and Exons

class Gene(object):
    def __init__(self, name, ENS_ID, gtype, chrom, start, end, strand):
        self.name = name 
        self.ENS_ID = ENS_ID # ensemble ID
        self.gtype = gtype
        self.chrom = chrom
        self.start = start # always smaller than end, direction is defined by strand
        self.end = end
        self.length = end - start
        self.strand = strand
        self.exon_list = [] # list of exons
        self.exon_count = len(self.exon_list)
        self.intron_list = [] # list of introns
####        self._flag

    def __str__(self):
        return 'name = %s, ENS_ID = %s, gtype = %s, chrom = %s, start = %d, end = %d, strand = %s, exon_count = %d' %(self.name, self.ENS_ID, self.gtype, self.chrom, self.start, self.end, self.strand, self.exon_count)

    def printExonList(self):
        for i in range(len(self.exon_list)):
            print(self.exon_list[i])

    def getIntronList(self):
        if len(self.exon_list) > 1:
            intronList = []
            for i in range(len(self.exon_list)-1):
                in_start = self.exon_list[i].end + 1
                in_end = self.exon_list[i+1].start - 1
                intron = Intron(in_start, in_end)
                intronList.append(intron)

            self.intron_list = intronList
        else:
            print('only one exon recorded for the gene')

            
    def printIntronList(self):
        for i in range(len(self.intron_list)):
            print(self.intron_list[i])

            
    def sortExonList(self):
        
        '''sorts exon list'''
        
        exon_sorted_list = []
        exon_start_list = []
        for i in range(len(self.exon_list)):
            exon_start_list.append(self.exon_list[i].start)

        exon_sorted_list = []
        inds = sorted(range(len(exon_start_list)), key = lambda k: exon_start_list[k])
        for i in range(len(inds)):
            exon_sorted_list.append(self.exon_list[inds[i]])

        self.exon_list = exon_sorted_list

            
    def exonOverlaps(self):

        '''reports overlapping exons'''
        #TODO: define and check the exon sort flag
        # to report if there are overlapping exons
        exonOverlapIncidents = 0
        for i in range(len(self.exon_list)-1):
            if self.exon_list[i].end > self.exon_list[i+1].start:
                exonOverlapIncidents +=1

        return(exonOverlapIncidents)
    

class Exon(object):
    def __init__(self, ENS_ID, start, end):
        self.ENS_ID = ENS_ID
        self.start = start
        self.end = end

    def __str__(self):
        return 'ENS_ID = %s, start = %d, end = %d' %(self.ENS_ID, self.start, self.end)

class Intron(object): # introns don't have ID, but just identified by exon_end+1 and exon_start-1
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __str__(self):
        return 'start = %d, end = %d' %(self.start, self.end)

# NOTE: once all the genes are recorded, go through the exons and save introns
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Genes and Exons

# reading the GTF file into gene struct - recording exons and introns
###############################

fileName = dataFolder + '/geneTypesFromGTF.pkl'
with open(fileName, 'rb') as pickledFile:
    geneTypesFromGTF = pickle.load(pickledFile)

geneTypeList = geneTypesFromGTF[0]    
acceptedGeneTypes = [geneTypeList[0],
                     geneTypeList[2],
                     geneTypeList[4],
                     geneTypeList[8],
                     geneTypeList[9],
                     geneTypeList[14]]

geneList = {} # donno why, want to call each gene by its ID
geneIDList = [] # the ID ordereded as we see them
i = 0
with open(gftFile, 'r') as coors:
    line = coors.readline()
    for l in range(4): # skipping header lines
        line = coors.readline()
        print(line)

    for line in coors:
#   i = 0 # ^^^^^^^^^^ for test
#    line = coors.readline() # ^^^^^^^^^^ for test
#    while i < 20: # ^^^^^^^^^^ for test

        infos = line.strip().split()
#        print(infos) #for test

        if infos[2] == 'gene':
            index = line.index('gene_type')
            thisType = line[index:].split('"')[1]
            
            if thisType in acceptedGeneTypes:
                record = True
#                i +=1 # ^^^^^^^^^^ for test
                
                index = line.index('gene_name')
                name = line[index:].split('"')[1]
                chrom = infos[0]
                start = int(infos[3])
                end = int(infos[4])
                strand = infos[6]
                gtype = thisType
                
                ENS_ID = infos[9]
                for x in enumerate(['"', ';']): ENS_ID = ENS_ID.replace(x[1], '')

                geneList[ENS_ID] = Gene(name, ENS_ID, gtype, chrom, start, end, strand)
                geneIDList.append(ENS_ID)
                #self.exon_list = [] # list of exons
                #self.exon_count = len(self.exon_list)
                #self.intron_list = [] # list of introns
                            
            else: 
                record = False

        elif record and infos[2] == 'exon': # we are reading other elements of the gene, looking for exon, if this is a gene that we are recording
#            if record: # 
#                if infos[2] == 'exon': # if this is an exon,
            if thisType == 'transcribed_unprocessed_pseudogene': # for this type of gene, we are only recording the processed_transcript
                if 'processed_transcript' in line:
                    index = line.index('exon_id')
                    ex_ID = line[index:].split('"')[1]
                    
                    if not ex_ID in geneList[ENS_ID].exon_list:
                        ex_end = int(infos[4])
                        ex_start = int(infos[3])
                        exon = Exon(ex_ID, ex_start, ex_end)
                        geneList[ENS_ID].exon_list.append(exon)

            else:
                index = line.index('exon_id')
                ex_ID = line[index:].split('"')[1]
                
                if not ex_ID in geneList[ENS_ID].exon_list:
                    ex_end = int(infos[4])
                    ex_start = int(infos[3])
                    exon = Exon(ex_ID, ex_start, ex_end)
                    geneList[ENS_ID].exon_list.append(exon)
                        
#        line = coors.readline() # ^^^^^^^^^^ for test

geneListsAndID = [geneList, geneIDList]

fileName = dataFolder + '/geneLists.pkl'
with open(fileName, 'wb') as f:
    pickle.dump(geneListsAndID, f)

fileName = dataFolder + '/geneLists.pkl'
with open(fileName, 'rb') as pickleFile:
    geneListsAndID = pickle.load(pickleFile)

geneIDList = geneListsAndID[1]
geneList = geneListsAndID[0]

# just testing
book = geneList[list(geneList.keys())[1234]]
print(book)
book2 = geneIDList[1234] 
print(book)

for i in range(10):

    book = geneList[list(geneList.keys())[i]]
    print(book)
    book.sortExonList()
    book.printExonList()
    book.getIntronList()
    book.printIntronList()
    print(book.exonOverlaps())

##TODO: exon overlap for genes, check, modify or others - I suspect that some none coding genes would have too many overlaps. But that's ok since they are a not many. My guess is that overlap happens much less in non coding genes

#############################################################    
# 3. Comparing the annotation labels to the genomic regions and transcriptomic data
#############################################################

import linecache


# just checking to see how many incidents of gene overlaps do we have
geneSS = np.zeros(len(geneIDList))
overlaps = 0
previous = 0
for i in range(len(geneIDList)):
   geneSS[i] = gene_start = geneList[geneIDList[i]].start
   if geneSS[i] < previous:
      overlaps += 1
   previous = geneList[geneIDList[i]].end

fig = plt.figure()
ax = plt.axes()
ax.plot(geneSS)
fig.show()
# 5941 where the start of the next gene is before the end of the previous gene (overlaps)


'''
This is to examine the enrichment the labels and classes in genomic regions, and also for the expressed versus
non expressed genes. I will keep the data in a numpy matrix of labels by regions, and also classes by regions

'''

'''
# TODO: load the list of classes and labels from the file
summaryAnnotation = {"classes" : classes, "labels": labels}

fileName = dataFolder + 'testBatch/' + sample + '/summaryAnnotation.pkl'
with open(fileName, 'wb') as f:
    pickle.dump(summaryAnnotation, f)
'''


fileName = dataFolder + '/geneLists.pkl' 
with open(fileName, 'rb') as pickleFile:
    geneLists = pickle.load(pickleFile)

geneList = geneLists[0]
geneIDList = geneLists[1]
del geneLists

fileName = dataFolder + 'testBatch/' + sample + '/summaryAnnotation.pkl'
with open(fileName, 'rb') as pickledFile:
    summaryAnnotation = pickle.load(pickledFile)

print(summaryAnnotation['classes'])
print(summaryAnnotation['labels'])

classes = summaryAnnotation['classes']
labels = summaryAnnotation['labels']

#classList = list(classes.keys())
classList = segwayLabels
classListIndMap = {}
for i in range(len(classList)):
   classListIndMap[classList[i]] = i
   
labelList = list(labels.keys())
classCount = 9 # this is the default for our annotations
labelCount = len(labelList) # this could vary for each annotation file
print(labelCount)

# TODO: sort the classList and label list

# TODO: the expression{} must have been saved, load it here
# 'expression' is a dictionary mapping gene ID to the expression level read in 0.2. 

extension = 3000 # count of basepairs monitored before and after the gene coordinates

cgi = 0 # walks on the geneIDList
ann_start = 0 # the start of the current annotation
ann_end = 0 # the end of the current annotation
ann_chr = 'chr'
ann_line_count = 0 # this is just to check the progress through the annotation file
previous_class = ''

'''

- I want to record the enrichment of labels in genomic regions. For this, I only need to record labels around genomic regions. I will need two matrices for this, one for labels and one for classes. I will be looking at 100 base pair resolution before and after genomic regions (30 fields) and 100 section for genomic regions (each section gets a percentage for each label)

- Then I want to record the enrichment of labels for the expressed versus not expressed genes. Again, I will have two matrices for each category. 


I will use two sets of matrices. First will be two 100 by class and label lists. 
This I will use for 3 expression levels: zero, zero and low, medium and high. low is log10(tpm) .5-1, medium is log10(tpm) 1-3, high is log10(tpm) 3+

The other matrices I will use for the exon and stuff.

'''

# the zero, zero and low, medium and high - see above
labelMats = [np.zeros((labelCount, 160)),np.zeros((labelCount, 160)),np.zeros((labelCount, 160))]
classMats = [np.zeros((classCount, 160)),np.zeros((classCount, 160)),np.zeros((classCount, 160))]

# this is to use for the negative strand genes
tempLabelMats = [np.zeros((labelCount, 160)),np.zeros((labelCount, 160)),np.zeros((labelCount, 160))]
tempClassMats = [np.zeros((classCount, 160)),np.zeros((classCount, 160)),np.zeros((classCount, 160))]

annLineInd = 1 # this keeps the annotation line index, the first line is empty, thus the index starts at 1
genomic_region_start_annLineInd = 0
firstGenomicAnn = True

previous_gene_chr = 'chr'
previous_extension_end = 0 # doesn't matter since chr condition is never true until the first gene is processed

# TODO: add annFile definition as annotation file

'''
* we define a genomic region as its start to end, plus the extension on each side
* we define a territory of a gene with regards to annotations as all the basepairs before the gene's extended
  sequence, as well as its actual genomic region (the extended sequence)
* once in each genomic region, we keep the annotationLineInd in the genomic_region_start_annLineInd. If the genomic region starts within 
  the previous genomic region, we set the annLineInd as the genomic_region_start_annLineInd

Walking on the genomic regions from the geneList using the geneIDList filling the enrichment data matrix:

-- walk through annotations, until we reach the actual genomic region of the current gene's territory:
    
--- record the labels for the genomic region until the end of genomic region (which is the end of genomic territory as well)

NOTICE THE STRAND MATTERS HERE!

'''


# while cgi < len(geneIDList): # modify the condition for the test runs
while cgi < 1: # >>>>>>>>>> TEST

   'we are in the gene territory'
   firstGenomicAnn = True
   
   geneID = geneIDList[cgi]
   gene_chr = geneList[geneIDList[cgi]].chrom

   gene_strand = geneList[geneID].strand
   gene_start = geneList[geneID].start
   gene_end = geneList[geneID].end
   gene_length = gene_end - gene_start
   gene_length_unit = int(gene_length/100)
   gene_length_last_unit = gene_length - (99* gene_length_unit)
   # TODO: something to fix: the gene_length_last_unit for negative strand versus positive strand

   extension_start = gene_start - extension
   extension_end = gene_end + extension

   ''' 
   if this gene starts somewhere in the preivous genomic region, 
   I will go back in the annotation file to the beginning of annotation for the previous gene

   '''
   if (gene_chr == previous_gene_chr) and (previous_extension_end > extension_start):
      annLineInd = genomicRegionStartAnnLineInd
   


   ''' picking the label/class matrix based on the gene expression level'''
   #TODO catch exception for when geneID is not in expression
   gene_exp = expression[geneID]
   if gene_exp == 0:
      expMatInd = 0
   elif gene_exp > 1:
      expMatInd = 2
   else:
      expMatInd = 1
      
   geneMatWalkIndex = 0
   previous_fill = 0 # at the begining of each gene, annotations are full

   # reading the next annotation
   line = linecache.getline(annFile, annLineInd)
   annLineInd +=1
   ann_line_count += 1
   fields = line.strip().split()
            
   ann_chr = fields[0]
   ann_start = int(fields[1])
   ann_end = int(fields[2])
        
   while ((ann_start < extension_end) or not(gene_chr == ann_chr)) and geneMatWalkIndex < 160: 
           
      '''
      NOTE: the second condition is for when we are at the end of a gene's territory, and then at the end of a chromosome, so when we go to the next gene, 
      ann_start is not smaller than extension_end, and we need to read the annotations until we are on the same chromosome

      The condition to be in one gene's territory (and not the next one), and it is for reading the annotations
      while in THIS gene territory, read annotations until we reach to the genomic region (plus extension)
      
      '''

      #  print('hereIam') #>>>>>>>>>> test
      #  line = annotations.readline()
            

      if ann_chr == gene_chr: # if we are in the genomic region (with extension)
         if ((ann_start < extension_end and ann_start > extension_start) or 
             (ann_end < extension_end and ann_end > extension_start) or
             (ann_start < extension_start and ann_end > extension_end)):

                   
            ''' We are in the genonimc region (with extension)'''
            
            if firstGenomicAnn:
               ''' Keeping the line index of the first annotation of the genomic region for overlapping genes'''
               genomicRegionStartAnnLineInd = annLineInd -1
               firstGenomicAnn = False


            ''' Taking off for filling the matrices... '''
            adjusted_ann_start = max(0, ann_start - extension_start)
            adjusted_ann_end = min(ann_end - extension_start, extension_end - extension_start)
            adjusted_ann_length = adjusted_ann_end - adjusted_ann_start

            ann_label = fields[3]
            ann_class = ann_label.split('_')[1]
            classInd = classListIndMap[ann_class]
            #labelInd = #TODO: get it
            
            # expMatInd
            if gene_strand == '+':
               while(adjusted_ann_length > 0) and geneMatWalkIndex < 30:
                  if (adjusted_ann_length >= 100*(1- previous_fill)):
                     classMats[expMatInd][classInd][geneMatWalkIndex] += 1 - previous_fill
                     adjusted_ann_length -= 100*(1- previous_fill)
                     previous_fill = 0
                     geneMatWalkIndex +=1
                  else:
                     classMats[expMatInd][classInd][geneMatWalkIndex] += adjusted_ann_length/100
                     previous_fill += adjusted_ann_length/100
                     adjusted_ann_length = 0
                     if previous_fill > 1:
                        previous_fill = 1

               while(adjusted_ann_length > 0) and geneMatWalkIndex < 129:
                  if (adjusted_ann_length >= gene_length_unit*(1- previous_fill)):
                     classMats[expMatInd][classInd][geneMatWalkIndex] += 1 - previous_fill
                     adjusted_ann_length -= gene_length_unit*(1 - previous_fill)
                     previous_fill = 0
                     geneMatWalkIndex +=1
                  else:
                     classMats[expMatInd][classInd][geneMatWalkIndex] += adjusted_ann_length/gene_length_unit
                     previous_fill += adjusted_ann_length/gene_length_unit
                     adjusted_ann_length = 0
                     if previous_fill > 1:
                        previous_fill = 1

               if (adjusted_ann_length > 0) and geneMatWalkIndex == 129:
                     if (adjusted_ann_length >= gene_length_last_unit*(1- previous_fill)):
                        classMats[expMatInd][classInd][geneMatWalkIndex] += 1 - previous_fill
                        adjusted_ann_length -= gene_length_last_unit*(1 - previous_fill)
                        previous_fill = 0
                        geneMatWalkIndex +=1
                     else:
                        classMats[expMatInd][classInd][geneMatWalkIndex] += adjusted_ann_length/gene_length_last_unit
                        previous_fill += adjusted_ann_length/gene_length_last_unit
                        adjusted_ann_length = 0
                        if previous_fill > 1:
                           previous_fill = 1
                         
               while(adjusted_ann_length > 0) and (geneMatWalkIndex < 160):
                  if (adjusted_ann_length >= 100*(1- previous_fill)):
                     classMats[expMatInd][classInd][geneMatWalkIndex] += 1 - previous_fill
                     adjusted_ann_length -= 100*(1- previous_fill)
                     previous_fill = 0
                     geneMatWalkIndex +=1
                  else:
                     classMats[expMatInd][classInd][geneMatWalkIndex] += adjusted_ann_length/100
                     previous_fill += adjusted_ann_length/100
                     adjusted_ann_length = 0
                     if previous_fill > 1:
                        previous_fill = 1


            if gene_strand == '-':
               while(adjusted_ann_length > 0) and geneMatWalkIndex < 30:
                  if (adjusted_ann_length >= 100*(1- previous_fill)):
                     classMats[expMatInd][classInd][159 - geneMatWalkIndex] += 1 - previous_fill
                     adjusted_ann_length -= 100*(1- previous_fill)
                     previous_fill = 0
                     geneMatWalkIndex +=1
                  else:
                     classMats[expMatInd][classInd][159 - geneMatWalkIndex] += adjusted_ann_length/100
                     previous_fill += adjusted_ann_length/100
                     adjusted_ann_length = 0
                     if previous_fill > 1:
                        previous_fill = 1
                        
               if (adjusted_ann_length > 0) and geneMatWalkIndex == 30:
                     if (adjusted_ann_length >= gene_length_last_unit*(1- previous_fill)):
                        classMats[expMatInd][classInd][159 - geneMatWalkIndex] += 1 - previous_fill
                        adjusted_ann_length -= gene_length_last_unit*(1 - previous_fill)
                        previous_fill = 0
                        geneMatWalkIndex +=1
                     else:
                        classMats[expMatInd][classInd][159 - geneMatWalkIndex] += adjusted_ann_length/gene_length_last_unit
                        previous_fill += adjusted_ann_length/gene_length_last_unit
                        adjusted_ann_length = 0
                        if previous_fill > 1:
                           previous_fill = 1


               while(adjusted_ann_length > 0) and geneMatWalkIndex < 129:
                  if (adjusted_ann_length >= gene_length_unit*(1- previous_fill)):
                     classMats[expMatInd][classInd][159 - geneMatWalkIndex] += 1 - previous_fill
                     adjusted_ann_length -= gene_length_unit*(1 - previous_fill)
                     previous_fill = 0
                     geneMatWalkIndex +=1
                  else:
                     classMats[expMatInd][classInd][159 - geneMatWalkIndex] += adjusted_ann_length/gene_length_unit
                     previous_fill += adjusted_ann_length/gene_length_unit
                     adjusted_ann_length = 0
                     if previous_fill > 1:
                        previous_fill = 1

                         
               while(adjusted_ann_length > 0) and (geneMatWalkIndex < 160):
                  if (adjusted_ann_length >= 100*(1- previous_fill)):
                     classMats[expMatInd][classInd][159 - geneMatWalkIndex] += 1 - previous_fill
                     adjusted_ann_length -= 100*(1- previous_fill)
                     previous_fill = 0
                     geneMatWalkIndex +=1
                  else:
                     classMats[expMatInd][classInd][159 - geneMatWalkIndex] += adjusted_ann_length/100
                     previous_fill += adjusted_ann_length/100
                     adjusted_ann_length = 0
                     if previous_fill > 1:
                        previous_fill = 1


                # if gene strand is positive, or negative, flag which label list we use
                # the only difference between the two strands is that I am going to reverse the index of
                # the columns in the matrix

      if geneMatWalkIndex < 160: 
         line = linecache.getline(annFile, annLineInd)
         annLineInd +=1
         ann_line_count += 1
         fields = line.strip().split()
          
         ann_chr = fields[0]
         ann_start = int(fields[1])
         ann_end = int(fields[2])
        

   cgi += 1 # next gene
   previous_extension_end = extension_end
   previous_gene_chr = gene_chr
   ###
   #TODO: we are probably missing something at the end of the last gene
   ###


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# the following lines are in the process of redoing using linecash instead of readline for the annotation file
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
with open(annFile, 'r') as annotations:
    
    '''
    * we define a genomic region as its start to end, plus the extension on each side
    * we define a territory of a gene with regards to annotations as all the basepairs before the gene's extended
      sequence, as well as its actual genomic region (the extended sequence)

    Walking on the genomic regions from the geneList using the geneIDList filling the enrichment data matrix:

    -- walk through annotations, until we reach the actual genomic region of the current gene's territory:
    
    --- record the labels for the genomic region until the end of genomic region (which is the end of genomic territory as well)
   

    NOTICE THE STRAND MATTERS HERE!

    '''

    #while cgi < len(geneIDList): # modify the condition for the test runs
    while cig < 1:

        'we are in the gene territory'

        geneID = geneIDList[cgi]
        gene_chr = geneList[geneIDList[cgi]].chrom

        gene_start = geneList[geneID].start
        gene_end = geneList[geneID].end
        gene_length = gene_end - gene_start
        gene_length_unit = int(gene_length/100)
        gene_length_last_unit = gene_length - (99* gene_length_unit)

        extension_start = gene_start - extension
        extension_end = gene_end + extension
        
        gene_strand = geneList[geneID].strand

        #TODO catch exception for when geneID is not in expression
        gene_exp = expression[geneID]
        if gene_exp == 0:
           expMatInd = 1
        elif gene_exp > 1:
           expMatInd = 2
        else:
           expMatInd = 1

        geneMatWalkIndex = 0
           
        while (ann_start < extension_end) or not(gene_chr == ann_chr):

           
            '''
            while in THIS gene territory, read annotations until we reach to the genomic region (plus extension)

            '''

            #  print('hereIam') #>>>>>>>>>> test
            #  line = annotations.readline()
            
            line = linecache.getline(annotations)
            ann_line_count += 1
            fields = line.strip().split()
            
            ann_chr = fields[0]
            ann_start = int(fields[1])
            ann_end = int(fields[2])

            if ann_chr == gene_chr: # if we are in the genomic region (with extension)
                if ((ann_start < extension_end and ann_start > extension_start) or 
                    (ann_end < extension_end and ann_end > extension_start) or
                    (ann_start < extension_start and ann_end > extension_end)):
                   
                   'we are in the genonimc region (with extension)'

                   adjusted_ann_start = max(0, ann_start - extension_start)
                   adjusted_ann_end = ann_end - extension_start
                   adjusted_ann_length = adjusted_ann_end - adjusted_ann_start


                   ann_label = fields[3]
                   ann_class = ann_label.split('_')[1]
                   classInd = classListIndMap[ann_class]
                   #labelInd = #TODO: get it

                   # expMatInd
                   while(adjusted_ann_length > 0) and geneMatWalkIndex < 30:
                      if adjusted_ann_length >= 100 and previous_fill == 1:
                         classMats[expMatInd][classInd][geneMatWalkIndex] += 1
                         adjusted_ann_length -== 100
                         previous_fill = 1
                      else:
                         previous_fill = adjusted_ann_length/100
                         classMats[expMatInd][classInd][geneMatWalkIndex] += previous_fill
                         adjusted_ann_length = 0

                   while(adjusted_ann_length > 0) and geneMatWalkIndex < 130:
                      if adjusted_ann_length >= gene_length_unit:
                         classMats[expMatInd][classInd][geneMatWalkIndex] += 1
                         adjusted_ann_length -= gene_length_unit
                      else:
                         classMats[expMatInd][classInd][geneMatWalkIndex] += adjusted_ann_length/gene_length_unit
                         adjusted_ann_length = 0
                         
                   while(adjusted_ann_length > 0) and geneMatWalkIndex < 160:
                      if adjusted_ann_length >= 100:
                         classMats[expMatInd][classInd][geneMatWalkIndex] += 1
                         adjusted_ann_length -== 100
                      else:
                         classMats[expMatInd][classInd][geneMatWalkIndex] += adjusted_ann_length/100
                         adjusted_ann_length = 0
                      
                   
                  # pick the matrix: zero, zero and low, medium and high based on the expression
                  # expMatInd has the index of the matrix

                  # fill the matrix based on the label and class: get the index of the class and label
                  ann_label = fields[3]
                  ann_class = ann_label.split('_')[1]
                  classInd = classListIndMap[ann_class]
                  labelInd = #TODO: get it

                # if gene strand is positive, or negative, flag which label list we use
                # the only difference between the two strands is that I am going to reverse the index of
                # the columns in the matrix
                if 'gene strand is positive':
                   
                   
                else: 
            
            ## idnetify the ratio for each label in the index
            ## += the matrix with the label ratios :

            ### if exp is zero 
            ### if exp is zero + low
            ### if exp is medium + high


   cgi += 1 # next gene

   ###
   #TODO: we are probably missing something at the end of the last gene
   ###




#########################################
# CODE DRAFT
#########################################

with open(gftFile, 'r') as coors:
    line = coors.readline()

gene_type = 'protein_coding'
gene_type = 'snoRNA'
file = open(gftFile, 'r')
while gene_type not in line:
    line = file.readline()

while gene_type in line:
    line = file.readline()

print(line)
line = file.readline()    


# extending the genomic region by <extension> for the start, based on the direction
if pcgenes.iloc[cgi].strand == '+':  
    start = int(pcgenes.iloc[cgi].start) - extension
    end = int(pcgenes.iloc[cgi].end) 
else:
    start = int(pcgenes.iloc[cgi].start)
    end = int(pcgenes.iloc[cgi].end + extension)

