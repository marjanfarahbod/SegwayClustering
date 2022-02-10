# This I am getting from QC_transcriptionComparison.py, but just the transcription file running
# This code is different than the original in that it works for a limited region of the genome, and for cluster labels only
#
# Segments:
# 0. Initials
# 1. Specifics
# 2. Main


########################################
# 0. Initials 
########################################

import linecache

# General data folder
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'

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

# .bed file


# expression file


########################################
# 2. Main
########################################



