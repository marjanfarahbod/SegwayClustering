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
## 1. compare to the transcript data

###############################
# 0. initials
###############################

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# data folder
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'

# the GTF file
gftFile = dataFolder + 'testBatch/gencode.v29.primary_assembly.annotation_UCSC_names.gtf'

# sample list (by folder name)
sample = 'H1'

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
    for line in coors: # reading
        index = line.index('gene_type')
        thisType = line[index:].split('"')[1]
        if not(thisType in geneTypeList):
            geneTypeList.append(thisType)

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
    Walking on the genomic coordinates and the annotation file doing:
    1. Filling up the annotation info
    2. Recording the annotations for the extended genomic regions

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

########################################
# 2. compare to the transcript data
########################################
# From the .bed files:
# 1. Enrichment of transcript bps in the genomic regions, and in the expressed regions: go through the file once and get this
# 2. Enrichment of promoter regions in 1000bps before the genomic regions, 1000bps before the expressed genes, and in the whole files
# (anything to do with the enhancers?)
# 3. Exploratory: for the 200bps start site and 200bp end site, enrichment of labels








#########################################
# CODE DRAFT
#########################################

# extending the genomic region by <extension> for the start, based on the direction
if pcgenes.iloc[cgi].strand == '+':  
    start = int(pcgenes.iloc[cgi].start) - extension
    end = int(pcgenes.iloc[cgi].end) 
else:
    start = int(pcgenes.iloc[cgi].start)
    end = int(pcgenes.iloc[cgi].end + extension)

