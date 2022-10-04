# for each run ID, search the portal,

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

dataSubFolder = 'classifier_data'


###################################################
# 1. Get the annotation list and make folders
###################################################

# 1.1 get results of the annotation list
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

wholeList_url = 'https://www.encodeproject.org/search/?type=Annotation&annotation_type=chromatin+state&lab.title=Maxwell+Libbrecht%2C+SFU&organism.scientific_name=Homo+sapiens&limit=all'
cookies = {'required_cookie': 'values'}
headers = {'accept' : 'application/json'}
cookies = {'_gid': 'GA1.1.1650294807.1655355136', '_ga': 'GA1.2.1811198692.1634853775', 'session': 'czQEFJeze3kbLDf0VUaY_gmPM4bVDfo3F4ZSTOFJkhAeTqwHY5tqZn-oaVtqpOp1R886l5FXCYuT2vi4HP9g71sxNjU1MzU1MTQzLCAxNjU1MzU1MTQwLjA2ODI3NjIsIHsiX2NzcmZ0XyI6ICI3NzRhMGJkZWMwMWU4NjUxOTBlN2Q4ZjlhMjc4Y2EwOTNjZmZmN2MyIiwgImF1dGgudXNlcmlkIjogIm1hcmphbi5mYXJhaGJvZEBnbWFpbC5jb20ifV0'} 
#response = requests.get(wholeList_url, cookies=cookies, headers=headers)
#wholeList_jdata = response.json()

# get the list of IDs

book = os.listdir(dataFolder + dataSubFolder)

runID_accession_map = {}
for item in book:
    identifier = item.split('_')[0]
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



