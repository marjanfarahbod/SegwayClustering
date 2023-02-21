# for each run ID, search the portal,

# 1. Get the annotation list and make folders
# 2. Get the tissue type, name and ID from ENCODE based on the track accession, runID
# 3. Get the tissue type, name and ID from ENCODE based on the track accession, Segway accession

###################################################
# 0. Initials
###################################################

from urllib.request import urlopen
import urllib
import json
import http.cookiejar
import requests
import re

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'

dataSubFolder = 'classifier_data' # likely to be obsolete
dataSubFolder = 'testBatch105/fromAPI/' # the 105Jan run
dataSubFolder = 'testBatch_May112022/' # the 92May run

###################################################
# 1. Get the annotation list and make folders
###################################################

# 1.1 get results of the annotation list
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

wholeList_url = 'https://www.encodeproject.org/search/?type=Annotation&annotation_type=chromatin+state&lab.title=Maxwell+Libbrecht%2C+SFU&organism.scientific_name=Homo+sapiens&limit=all'
cookies = {'required_cookie': 'values'}
headers = {'accept' : 'application/json'}
cookies = {'_gid': 'GA1.1.1312437877.1676396824', '_ga': 'GA1.1.1807347774.1675726254', 'session': 'iQbtInJ0gMtIfk3wvDuznxURjP2XyOBrhSztxHmXjEwLjaO-WAWZdInJ1vANxpjAMs728SG8__llxKZIbJi2iVsxNjc1NzI2MjYwLCAxNjc1NzI2MjU2LjU3NDQxMzMsIHsiX2NzcmZ0XyI6ICI0ZGJiZGIxOTEwMmY0MTBjOGE2NTcwZTc4ZjYzYjIzZWY3NDY4YzRlIiwgImF1dGgudXNlcmlkIjogIm1hcmphbi5mYXJhaGJvZEBnbWFpbC5jb20ifV0'} 
#response = requests.get(wholeList_url, cookies=cookies, headers=headers)
#wholeList_jdata = response.json()

# get the list of IDs

book = os.listdir(dataFolder + dataSubFolder)

runID_accession_map = {}
for item in book:
    identifier = item.split('_')[0] # I need it only if I have extra items in the runID
    print(identifier)

    url = 'https://www.encodeproject.org/search/?type=Annotation&analyses.aliases=encode-processing-pipeline:' + identifier
    response = requests.get(url, cookies = cookies, headers=headers)

    search_result = response.json()

    if len(search_result) < 1:
        continue
    else:
        result = search_result['@graph'][0]

    accession = result['accession']
    print(accession)

    runID_accession_map[identifier] = accession

    outputFile = dataFolder + 'runID_accession_map_105run.pkl'
    with open(outputFile, 'wb') as output:
        pickle.dump(runID_accession_map, output)


# 2. For the runs that are not on ENCODE, based on the track data, get the tissue info
########################################

# here, there is going to be no accession, but only the sample tissue and tissue ID stuff

dataSubFolder = 'testBatch_May112022/' # the 92May run
book = os.listdir(dataFolder + dataSubFolder) # this needs filtering if there are other files, I need to modify

runID_tissueInfo_map = {}
for identifier in book:
    print(identifier)

    # open the track assay file and read the line with the H3K4me1 ID
    sampleFolderAdd = dataFolder + dataSubFolder + identifier + '/'

    try:
        mapping_file = sampleFolderAdd + 'trackname_assay.txt'
        with open(mapping_file) as inputFile:
            for line in inputFile:
                fields = line.strip().split()
                if fields[1] == 'H3K4me1':
                    trackID = fields[0]
                    break

        url = 'https://www.encodeproject.org/experiments/%s/?format=json' %(trackID)

        response = requests.get(url, cookies=cookies, headers=headers)
        jsonData = response.json()
        
        term_name = jsonData['biosample_ontology']['term_name']
        termID = jsonData['biosample_ontology']['term_id']
        term_class = jsonData['biosample_ontology']['classification']
        term_name_mod = term_name.replace(' ', '+')

        tissue_info = (term_name, termID, term_class)
        runID_tissueInfo_map[identifier] = tissue_info

    except NotADirectoryError: 
        print(mapping_file)
        print('not a dir')


outputFile = dataFolder + dataSubFolder + 'runID_tissueInfo.pkl'
with open(outputFile, 'wb') as output:
    pickle.dump(runID_tissueInfo_map, output)


# 3. Get the tissue type, name and ID from ENCODE based on the track accession, Segway accession
########################################

dataSubFolder = 'testBatch105/fromAPI/' # the 105Jan run
inputFileName = 'all_annInfo_list.pkl'
inputFile = dataFolder + dataSubFolder + inputFileName
with open(inputFile, 'rb') as f:
    annInfo_list = pickle.load(f)

accession_tissueInfo_map = {}
tissue_categories = []
for ann in annInfo_list:

    index = ann['index']
    print(index)
    sampleFolder =  ann['accession']

    sampleFolderAdd = dataFolder + dataSubFolder + sampleFolder + '/'

    mapping_file = sampleFolderAdd + 'trackname_assay.txt'
    with open(mapping_file) as inputFile:
        for line in inputFile:
            fields = line.strip().split()
            if fields[1] == 'H3K4me1':
                trackID = fields[0]
                break

    url = 'https://www.encodeproject.org/experiments/%s/?format=json' %(trackID)

    response = requests.get(url, cookies=cookies, headers=headers)
    jsonData = response.json()
        
    term_name = jsonData['biosample_ontology']['term_name']
    termID = jsonData['biosample_ontology']['term_id']
    term_class = jsonData['biosample_ontology']['classification']
    term_name_mod = term_name.replace(' ', '+')

    tissue_info = (term_name, termID, term_class)
    accession_tissueInfo_map[sampleFolder] = tissue_info

outputFile = dataFolder + dataSubFolder + 'segwayAccession_tissueInfo.pkl'
with open(outputFile, 'wb') as output:
    pickle.dump(accession_tissueInfo_map, output)

tissue_categories = []
for sample in accession_tissueInfo_map.values():
    term_class = sample[2]
    if not(term_class in tissue_categories):
        tissue_categories.append(term_class)
    
