# The goal here is to just do structure for the phantom5 enhancer regions so that I can adopt the expression code for it.
#

# report on enhancer coverage throughout the genome.

# can we predict enhancer bay enrichment of labels around them? This would solve the starting site segway.


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

inputFile = dataFolder + dataSubFolder + 'metaInfo.pkl'
with open(inputFile, 'rb') as f:
    annMeta = pickle.load(f)

sample_count = len(annMeta)

annAccessionList = list(annMeta.keys())
annAccession = annAccessionList[104]
print(annAccession)

# 0. process one of the files (phanotm5 is not until March31)
#########################################

