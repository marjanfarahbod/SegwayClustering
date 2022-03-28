# this is the code to get the necessary data for the 105 test samples
#
# pre-code stuff:
# 1. loging to the page for the samples. DONE
# 2. get JSON view. DONE
#
# I think I will need for each annotation, a folder. But then I am going to need a key file, where I have the donor, tissue and the annotation key - so I can look at annotations by tissue. 
#
# coding: 
# 1. In the JSON of the main page, get the annoation ID - we need to look into the page for each annotation ID. - get the files and the meta from GCP files. 
# 2. In each page of the annotation, get the biosample_ontology.term_name (I think you need to get the donor ID as well)
# 3. Via code, get the serach result for the donor + ontology term.
# 4. In the search result, look for RNA-seq
# 5. Go to the page of the RNA-seq
# 6. Get the RNAseq data
#
#
###################################################
# 0. Initials
###################################################

from urllib.request import urlopen
import urllib
import json
import http.cookiejar
import requests

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'

inputFolder = 'testBatch105/'

###################################################
# 1. In the JSON of the main page, get the annoation ID - we need to look into the page for each annotation ID.
###################################################

# get the JSON file: since this is a login 

url = 'https://www.encodeproject.org/search/?type=Annotation&annotation_type=chromatin+state&lab.title=Maxwell+Libbrecht%2C+SFU&organism.scientific_name=Homo+sapiens'

url = 'https://www.encodeproject.org/annotations/ENCSR596QDJ/?format=json'
response = urlopen(url)


cookies = {'required_cookie': 'values'}
headers = {'accept' : 'application/json'}
cookies = {'_gid': 'GA1.1.827781495.1645033198', '_ga': 'GA1.2.1811198692.1634853775', 'session': 'R08tk2RYfGA-eClOXepXPUohq40r62Q-t5AqIOj78ZG1pavyA_isX4zDkoa8Z1DK0l1dgsh3YqjxSYNvFW-EhVsxNjQ1MDUwOTI0LCAxNjQ1MDMzMjAzLjIwOTE3NjMsIHsiX2NzcmZ0XyI6ICJiN2Q3NDZiYTdiNmMzMzUwMGE0Y2EwMDE5Yjg3YzBlNDRkYWFiMWQ3IiwgImF1dGgudXNlcmlkIjogIm1hcmphbi5mYXJhaGJvZEBnbWFpbC5jb20ifV0'} # I need to be logged in, get it from the cookie thing

response = requests.get(url, cookies=cookies, headers=headers)
meta = response.json()
print(meta['aliases'])
print(meta['biosample_ontology']['term_name'])


wholeList_url = 'https://www.encodeproject.org/search/?type=Annotation&annotation_type=chromatin+state&lab.title=Maxwell+Libbrecht%2C+SFU&organism.scientific_name=Homo+sapiens'

response = requests.get(wholeList_url, cookies=cookies, headers=headers)
jsonData = response.json()

# I can't do the above cause these pages need me logged in - and I don't want to pass my googel password through. I will work with the data from GCP for this part for now, until I find a solution for that.




