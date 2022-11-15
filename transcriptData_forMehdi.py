

from QC_transcriptionComparison_util import Gene, Exon, Annotation, AnnotationClass

# 0. the gene coordinate file exploration
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# the pickled gene coordinate data (same for all samples)
fileName = dataFolder + '/geneLists.pkl'
with open(fileName, 'rb') as pickleFile:
    geneListsAndIDs = pickle.load(pickleFile)
    
geneList = geneListsAndIDs[0] # this is a dictionary of gene ID to gene objects (each object has gene info)
geneIDList = geneListsAndIDs[1] # this is list of gene IDs
del geneListsAndIDs

print(geneIDList[0]) # what is printed is the gene's ensemble ID
print(geneList[geneIDList[0]]) # what is printed is gene info, see below

'''
name = DDX11L1, ENS_ID = ENSG00000223972.5, gtype = transcribed_unprocessed_pseudogene, chrom = chr1, start = 11869, end = 14409, strand = +, exon_count = 0
'''

print(geneList[geneIDList[0]].name)
print(geneList[geneIDList[0]].ENS_ID)
print(geneList[geneIDList[0]].chrom)
print(geneList[geneIDList[0]].start)
print(geneList[geneIDList[0]].end)
print(geneList[geneIDList[0]].strand) # direction matters for placement of promoter


# 1. the transcript file 
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# the pickled transcript file - you have multiple transcript files
# NOTE: transcript levels are log10 transformed 
transcriptFile = dataFolder + 'Mehdi_testRun/tomarjan/RNA_seq/K562/geneExp_dict_ENCFF840UYD.pkl'
with open(transcriptFile, 'rb') as f:
    explevels = pickle.load(f)

print(len(explevels)) # this is the expression level of all genes in that sample

# 2. get the expression level and coordinates of all genes in a sample
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

for geneID in geneIDList[0:5]:

    # gene coordinate info:
    print('______')
    print(geneList[geneID].name)
    print(geneList[geneID].ENS_ID)
    print(geneList[geneID].chrom)
    print(geneList[geneID].start)
    print(geneList[geneID].end)
    print(geneList[geneID].strand)

    # gene expression level
    print(explevels[geneID])
