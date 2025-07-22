# Remains unfinished until I figure how can I get the RNAseq data from the two batch of samples, for the two batches of May and the other one. The other one has the samples I guess. I will go over the command for the other ones. 
# To get the RNAseq data from portal, we need to search with two keys, the donorID and the term_name for the sample. These I can get from any of the track data for the sample.
# 0. Initial
# 1. For a list of annotations, fetch the trackIDs
# 2. For each of the trackIDs, fetch the donorID and the term_mod
# 3. Get the annotation
# 4. while you are at it, get the tissue info as well from the track data


###################################################
# 0. Initials
###################################################

from urllib.request import urlopen
import urllib
import json
import http.cookiejar
import requests
import re
import pickle

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'

# I already have the data for the May runs. I need the data for the other two batches, and I want to get them using the track data only. 
dataSubFolder = '/testBatch_May112022/'

###################################################
# 1. get the RNAseq and tissue info
###################################################

# load that file and map the runID to the donorID
metaInfoFile = dataFolder + dataSubFolder + '92segwaymay2022.tsv'
runID_donorID_map = {}
with open(metaInfoFile, 'r') as f:
    for line in f:
        content = line.strip().split()
        runID = content[0]
        donorID = content[3].split('_')[0]
        runID_donorID_map[runID] = donorID

# TODO: get the accession for these files

# load the runID for these samples
inputFile = dataFolder + dataSubFolder + 'runID_list.pkl'
with open(inputFile, 'rb') as f:
    runID_list = pickle.load(f)

# Just checking if the two runID lists match
for runID in runID_list:
    print(runID_donorID_map[runID])

# get list of samples based on the folder
index = 0
runID_tissue = {}
for runID in runID_list:
    print(index)

    sampleFolder = dataFolder + dataSubFolder + runID + '/'
    trackFile = sampleFolder + 'trackname_assay.txt'

    with open(trackFile, 'r') as tf:
        trackAccession = tf.readline().split('\t')[0]

    url = 'https://www.encodeproject.org/experiments/%s/?format=json' %(trackAccession)

    response = requests.get(url)
    jsonData = response.json()
    
    termName_id_list = [jsonData['biosample_ontology']['term_name'], jsonData['biosample_ontology']['@id'], jsonData['biosample_ontology']['classification'], runID_donorID_map[runID]]

    term_name = jsonData['biosample_ontology']['term_name']
    runID_tissue[runID] = termName_id_list

    term_name_mod = term_name.replace(' ', '+')

    donorID = runID_donorID_map[runID]
    
    search_url = 'https://www.encodeproject.org/search/?type=Experiment&biosample_ontology.term_name=%s&status=released&assay_title=total+RNA-seq&replicates.library.biosample.donor.accession=%s' %(term_name_mod ,donorID)
    print(search_url)
    
    response = requests.get(search_url)
    #search_results = response.json()

    if response.status_code == '404':
        continue
    else:
        for result in search_results['@graph']:
            RNA_accession = result['accession']
            print('annotation accession' + accession)
            print(RNA_accession)
        
            # go to the RNA page and download the gene quantity files in the folder
#            rna_url = 'https://www.encodeproject.org/experiments/ENCSR052FJA/'
            rna_url = 'https://www.encodeproject.org/experiments/%s' %(RNA_accession)
            response = requests.get(rna_url)
            rna_json = response.json()

            rna_fileList = rna_json['files']

            for file in rna_fileList:
                if file['output_type'] == 'gene quantifications':
                    
                    file_accession = file['accession']
                    downloadURL = 'https://www.encodeproject.org/files/%s/@@download/%s.tsv' %(file_accession, file_accession)
                    if 'preferred_default' in file:
                        if file['preferred_default'] == True:
                            downloadPath = dataFolder + dataSubFolder + accession + '/preferred_default_' + file_accession + '.tsv'
                    else:
                        downloadPath = dataFolder + dataSubFolder + accession + '/' + file_accession + '.tsv'
                        
                    print('doing the download')
                    urllib.request.urlretrieve(downloadURL, downloadPath)

# just recording the IDs and stuff.

file = dataFolder + dataSubFolder + 'runID_tissue_donor.pkl'
with open(file, 'wb') as f:
    pickle.dump(runID_tissue, f)




