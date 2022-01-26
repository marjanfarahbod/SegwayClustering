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
# 1. compare to the transcript data
########################################
# From the .bed files:
# 1. Enrichment of transcript bps in the genomic regions, and in the expressed regions: go through the file once and get this
# 2. Enrichment of promoter regions in 1000bps before the genomic regions, 1000bps before the expressed genes, and in the whole files
# (anything to do with the enhancers?)
# 3. Exploratory: for the 200bps start site and 200bp end site, enrichment of labels 

coors_frame
exp_frame
annFile = dataFolder + 'testBatch/' + sample + '/ccSegway.bed'

# get all the protein coding genes from the coors
pcgenes = coors_frame[coors_frame['geneType'].str.contains('protein_coding')]

# subset the file for genomic regions and 1500bps before and after
label = namedtuple('label', ['called', 'category', 'color', 'bp_count', 'region_count'])
ann_labels = [] # unique combination of cluster + classes

cgi = 0
with open(annFile, 'r') as annotations: # two things I do in this following loop: fill up the label info, subset annotations for genomic regions

    if pcgenes.iloc[cgi].strand == '+': # the current genomic region with start and end, regions are sorted
        start = pcgenes.iloc[cgi].start - 1500
        end = pcgenes.iloc[cgi].end
    else:
        start = pcgenes.iloc[cgi].start
        end = pcgenes.iloc[cgi].end + 1500

    for line in annotations:

        #fill up the annotation info
        fields = line.strip().split()
        
        
    

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 1.1 get the labels 
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>



