# this is for the 105 runs, to get the transcript files (for those with transcript data) from the protal. Here I am starting from where I have the medatadata JSON file for the annotations
# starting from a list of folders, each has a json file in them with the metadata of the annotations. Notice that at this stage I have all the segway and classifier output in these folders, I am only checking the samples to find the ones with RNAseq data so I can run the RNA-seq QC part for them.
# 
#
#
###################################################
# 0. Initials
###################################################

from urllib.request import urlopen
import json

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'

inputFolder = 'testBatch105/'

###################################################
# 1. In the JSON of the main page, get the annoation ID - we need to look into the page for each annotation ID.
###################################################


JSONFile = dataFolder + inputFolder + 'metadata.json'
f = open(JSONFile)
metadata = json.load(f)

# now extract the biosample and donor ID - then search the website for RNAseq sampls for these two






