# to get for each sample from the sample ID, name of the track list


from urllib.request import urlopen
import urllib
import json
import http.cookiejar
import requests
import re

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
dataSubFolder = 'testBatch105/fromAPI/'

# I don't need this url and the search response for this code, just copy pasted from before
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
wholeList_url = 'https://www.encodeproject.org/search/?type=Annotation&annotation_type=chromatin+state&lab.title=Maxwell+Libbrecht%2C+SFU&organism.scientific_name=Homo+sapiens&limit=all'
#cookies = {'required_cookie': 'values'}
headers = {'accept' : 'application/json'}
cookies = {'_gid': 'GA1.1.466526283.1652932878', '_ga': 'GA1.1.1847928511.1604510440', 'session': 'e5pxO0f-5mChVJq7fu7_UMvxSse3sVnbUuxeiQQOkQnuZklLJ5x-vZGxbAKIj9A5Z9oSRwRAKQIkJnX3Lm9TlVsxNjUyOTMyODgzLCAxNjUyNzIzNjUzLjIyNzgyNTYsIHsiX2NzcmZ0XyI6ICIwOWVkZmEyMmU5ZDAzYTBhNzcxYmUxYjdmNGFlMDYyYTc2MGNjODllIiwgImF1dGgudXNlcmlkIjogIm1hcmphbi5mYXJhaGJvZEBnbWFpbC5jb20ifV0'} 
response = requests.get(wholeList_url, cookies=cookies, headers=headers)
wholeList_jdata = response.json()
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# list of annotation folders
inputFile = dataFolder + dataSubFolder + 'metaInfo.pkl'
with open(inputFile, 'rb') as pickledFile:
    annMeta = pickle.load(pickledFile)

accession_list = list(annMeta.keys())

c = 0
for accession in sampleFolder_list:
    print(c)
    print(accession)
    c += 1
    
    url = 'https://www.encodeproject.org/annotations/%s/?format=json' %(accession)
    response = requests.get(url, cookies=cookies, headers=headers)
    jsonData = response.json()

    trackFile_list = jsonData['contributing_files']

    track_assay_map = {}
    print('tracks')
    for file in trackFile_list:

        if file.split('/')[2].endswith('sizes'):
            continue
        
        file_accession = file.split('/')[2]

        file_url = 'https://www.encodeproject.org/experiments/%s' %(file_accession)
        file_response = requests.get(file_url, cookies=cookies, headers=headers)
        file_json = file_response.json()

        print(file_accession)

        if 'target' in file_json:
            assay = file_json['target']['label']

        else:
            assay = file_json['assay_title']

        track_assay_map[file_accession] = assay

    sampleFolder = dataFolder + dataSubFolder + accession + '/'
    outputFile = sampleFolder + 'trackname_assay.txt'
    with open(outputFile, 'w') as f:
        for track in list(track_assay_map.keys()):
            f.write('%s\t%s\n' %(track, track_assay_map[track]))
        


