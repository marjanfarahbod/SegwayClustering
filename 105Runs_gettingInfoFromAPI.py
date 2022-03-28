# this is the code to get the necessary data for the 105 test samples
#
# pre-code stuff:
# 1. loging to the page for the samples. DONE
# 2. get JSON view. DONE
#
# I think I will need for each annotation, a folder. But then I am going to need a key file, where I have the donor, tissue and the annotation key - so I can look at annotations by tissue.
#
# 1. Get the annotation lists and make folders
# 2. Get the RNA-seq files 
# 3. Get the rest of files
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
import re

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'

dataSubFolder = 'testBatch105/'


###################################################
# 1. Get the annotation list and make folders
###################################################

# 1.1 get results of the annotation list
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

wholeList_url = 'https://www.encodeproject.org/search/?type=Annotation&annotation_type=chromatin+state&lab.title=Maxwell+Libbrecht%2C+SFU&organism.scientific_name=Homo+sapiens&limit=all'
cookies = {'required_cookie': 'values'}
headers = {'accept' : 'application/json'}
cookies = {'_gid': 'GA1.1.621565880.1648437766', '_ga': 'GA1.2.1811198692.1634853775', 'session': 'r3Uv1ma0Q01rs-Q3brkRSwkJhyUMdFZtDQioXhmgUgCNKBGpJYyAOHL4dQ3IgYgKCnOdTmmNknmK7sTwDqtIiFsxNjQ4NDM3Nzc2LCAxNjQ4MTU2NTE3LjY0ODY0NTIsIHsiX2NzcmZ0XyI6ICJkYjE4M2MzNTllMmY5OWYzNzg0MWFjOWVmNWM0NGM0ODAyNDk0MTFhIiwgImF1dGgudXNlcmlkIjogIm1hcmphbi5mYXJhaGJvZEBnbWFpbC5jb20ifV0'} 
response = requests.get(wholeList_url, cookies=cookies, headers=headers)
wholeList_jdata = response.json()

# 1.2 get the annotation IDs and make folders
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
annotations_IDlist = []
for annotation in wholeList_jdata['@graph']:
    accession = annotation['accession']
    annotations_IDlist.append(accession)

    path = dataFolder + dataSubFolder + accession
    os.mkdir(path)


###################################################
# 2. Get the RNA-seq files 
###################################################

# coding: 
# 1. In the JSON of the main page, get the annoation ID - we need to look into the page for each annotation ID. - get the files and the meta from GCP files. 
# 2. In each page of the annotation, get the biosample_ontology.term_name (I think you need to get the donor ID as well)
# 3. Via code, get the serach result for the donor + ontology term.
# 4. In the search result, look for RNA-seq
# 5. Go to the page of the RNA-seq
# 6. Get the RNAseq data

# get the JSON file: passing cookies since the access requires login
# get the cookies from the browser


# from the jsonData, get the list of IDs
i = 0
for accession in annotations_IDlist:

    print('accession count %d' %(i))
    i +=1

    # get the json meta for the annotation
    url = 'https://www.encodeproject.org/annotations/%s/?format=json' %(accession)
    response = requests.get(url, cookies=cookies, headers=headers)
    jsonData = response.json()

    # TODO download all the files you need for running the classifier and QC


    # get the RNA-seq data if present
    aliases = jsonData['aliases'][0]
    index = aliases.index(':')
    
    donorID = aliases[index+1:].split('_')[0]

    term_name = jsonData['biosample_ontology']['term_name']
    term_name_mod = term_name.replace(' ', '+')

    search_url = 'https://www.encodeproject.org/search/?type=Experiment&biosample_ontology.term_name=%s&status=released&assay_title=total+RNA-seq&replicates.library.biosample.donor.accession=%s' %(term_name_mod ,donorID)

    response = requests.get(search_url, cookies=cookies, headers=headers)
    search_results = response.json()

    # TODO: what if len('@graph')>1 ?

    if len(search_results['@graph']) > 3:
        print('more than 3 RNAseq')
    else:
        for result in search_results['@graph']:
            RNA_accession = result['accession']
            print('annotation accession' + accession)
            print(RNA_accession)
        
            # go to the RNA page and download the gene quantity files in the folder
#            rna_url = 'https://www.encodeproject.org/experiments/ENCSR052FJA/'
            rna_url = 'https://www.encodeproject.org/experiments/%s' %(RNA_accession)
            response = requests.get(rna_url, cookies=cookies, headers=headers)
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

                    
##################################################                    
# 3. Get the rest of files
##################################################

# the default bedbed9 has the classifier/cluster label: download it
# the non default bedbed9 page has the extra files, download them! 

i = 0
for accession in annotations_IDlist:

    print('accession count %d' %(i))
    i +=1

    # get the json meta for the annotation
    url = 'https://www.encodeproject.org/annotations/%s/?format=json' %(accession)
    response = requests.get(url, cookies=cookies, headers=headers)
    jsonData = response.json()

    #
    fileList = jsonData['files']
    for file in fileList:
        if file['file_format'] == 'bed':
            print('is bed')
            file_accession = file['accession']
            if 'preferred_default' in file:
                print('downloading')
                downloadURL = 'https://www.encodeproject.org/files/%s/@@download/%s.bed.gz' %(file_accession, file_accession)
                bedFileName = downloadURL.split('/')[-1]
                downloadPath = dataFolder + dataSubFolder + accession + '/' + bedFileName
                
                # urllib.request.urlretrieve(downloadURL, downloadPath) # doesnt work since I need to pass cookies

                r = requests.get(downloadURL, cookies=cookies, headers=headers)
                with open(downloadPath, 'wb') as f:
                    f.write(r.content)

            else: # go to the page and download all the other files
                url = 'https://www.encodeproject.org/annotations/%s/?format=json' %(file_accession)
                response = requests.get(url, cookies=cookies, headers=headers)
                bedFile_json = response.json()

                id = bedFile_json['quality_metrics'][0]['@id']
                encode_section = 'https://www.encodeproject.org'

                
                file_section = '/signal_distribution_tab/signal_distribution.tab.txt'
                file_name = file_section.split('/')[-1]

                downloadURL = '%s%s@@download%s' %(encode_section, id, file_section)
                downloadPath = dataFolder + dataSubFolder + accession + '/' + file_name
                
                r = requests.get(downloadURL, cookies=cookies, headers=headers)
                with open(downloadPath, 'wb') as f:
                    f.write(r.content)
                    

                file_section = '/feature_aggregation_tab/feature_aggregation.tab.txt'
                file_name = file_section.split('/')[-1]

                downloadURL = '%s%s@@download%s' %(encode_section, id, file_section)
                downloadPath = dataFolder + dataSubFolder + accession + '/' + file_name
                
                r = requests.get(downloadURL, cookies=cookies, headers=headers)
                with open(downloadPath, 'wb') as f:
                    f.write(r.content)

                    
                file_section = '/length_distribution_tab/length_distribution.tab.txt'
                file_name = file_section.split('/')[-1]

                downloadURL = '%s%s@@download%s' %(encode_section, id, file_section)
                downloadPath = dataFolder + dataSubFolder + accession + '/' + file_name
                
                r = requests.get(downloadURL, cookies=cookies, headers=headers)
                with open(downloadPath, 'wb') as f:
                    f.write(r.content)

                    
                file_section = '/segment_sizes_tab/segment_sizes.tab.txt'
                file_name = file_section.split('/')[-1]

                downloadURL = '%s%s@@download%s' %(encode_section, id, file_section)
                downloadPath = dataFolder + dataSubFolder + accession + '/' + file_name
                
                r = requests.get(downloadURL, cookies=cookies, headers=headers)
                with open(downloadPath, 'wb') as f:
                    f.write(r.content)


##################################################
# DRAFT
##################################################
                    
        for item in book:
            accession = item['accession']
            if accession == 'ENCFF330DCY':
                print(i)
            i += 1

        # item 29 is the one I need
        myfile = book[29]
        print(myfile['output_type'])

        downloadURL = 'https://www.encodeproject.org/files/ENCFF330DCY/@@download/ENCFF330DCY.tsv'
        print('doing the download')
        urllib.request.urlretrieve(downloadURL, '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/testTS.tsv')

    
# for each file, search to see if you have an RNAseq data, for the search we need the
'''
       get ENCOD (donor) from a.aliases
        get term_name from biosample_ontology.term_name
        res = requests.get(https://www.encodeproject.org/search/?type=Experiment&biosample_ontology.term_name=spleen&status=released&assay_title=total+RNA-seq&replicates.library.biosample.donor.accession=ENCDO451RUA)
        # list of experiment datasets (c.a. 1) -> ENCSRxxx123 (@graph[0].accession)
        res2 = requests.get(https://www.encodeproject.org/experiments/{accession}
        res.json experiment.files list (output_type = 'gene quantification' and preferred_default=true ) // 1 file
'''

url = 'https://www.encodeproject.org/search/?type=Annotation&annotation_type=chromatin+state&lab.title=Maxwell+Libbrecht%2C+SFU&organism.scientific_name=Homo+sapiens'

wholeList_url = 'https://www.encodeproject.org/search/?type=Annotation&annotation_type=chromatin+state&lab.title=Maxwell+Libbrecht%2C+SFU&organism.scientific_name=Homo+sapiens'
wholeList_url = 'https://www.encodeproject.org/search/?type=Annotation&annotation_type=chromatin+state&lab.title=Maxwell+Libbrecht%2C+SFU&organism.scientific_name=Homo+sapiens&limit=all'
cookies = {'required_cookie': 'values'}
headers = {'accept' : 'application/json'}
cookies = {'_gid': 'GA1.1.827781495.1645033198', '_ga': 'GA1.2.1811198692.1634853775', 'session': 'R08tk2RYfGA-eClOXepXPUohq40r62Q-t5AqIOj78ZG1pavyA_isX4zDkoa8Z1DK0l1dgsh3YqjxSYNvFW-EhVsxNjQ1MDUwOTI0LCAxNjQ1MDMzMjAzLjIwOTE3NjMsIHsiX2NzcmZ0XyI6ICJiN2Q3NDZiYTdiNmMzMzUwMGE0Y2EwMDE5Yjg3YzBlNDRkYWFiMWQ3IiwgImF1dGgudXNlcmlkIjogIm1hcmphbi5mYXJhaGJvZEBnbWFpbC5jb20ifV0'} 
response = requests.get(wholeList_url, cookies=cookies, headers=headers)
wholeList_jdata = response.json()

# from the jsonData, get the list of IDs
annotations_IDlist = []
for annotation in wholeList_jdata['@graph']:
    accession = annotation['accession']
    annotations_IDlist.append(accession)

    # TODO make a folder for the annotation, if it doesn't exist already
    
    url = 'https://www.encodeproject.org/annotations/%s/?format=json' %(accession)
    response = requests.get(url, cookies=cookies, headers=headers)
    jsonData = response.json()

    # TODO download all the files you need for running the classifier and QC

    aliases = jsonData['aliases'][0]
    index = aliases.index(':')
    
    donorID = aliases[index+1:].split('_')[0]

    term_name = jsonData['biosample_ontology']['term_name']
    term_name_mod = term_name.replace(' ', '+')

    search_url = 'https://www.encodeproject.org/search/?type=Experiment&biosample_ontology.term_name=%s&status=released&assay_title=total+RNA-seq&replicates.library.biosample.donor.accession=%s' %(term_name_mod ,donorID)

    response = requests.get(search_url, cookies=cookies, headers=headers)
    search_results = response.json()

    # TODO: what if len('@graph')>1 ?

    for result in search_results['@graph']:
        RNA_accession = result['accession']
        print(RNA_accession)
        
        # go to the RNA page and download the gene quantity files in the folder

        rna_url = 'https://www.encodeproject.org/experiments/ENCSR052FJA/'
        response = requests.get(rna_url, cookies=cookies, headers=headers)
        rna_json = response.json()

        book = rna_json['files']

        i = 0
        for item in book:
            accession = item['accession']
            if accession == 'ENCFF330DCY':
                print(i)
            i += 1

        # item 29 is the one I need
        myfile = book[29]
        print(myfile['output_type'])

        downloadURL = 'https://www.encodeproject.org/files/ENCFF330DCY/@@download/ENCFF330DCY.tsv'
        print('doing the download')
        urllib.request.urlretrieve(downloadURL, '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/testTS.tsv')

    
# for each file, search to see if you have an RNAseq data, for the search we need the
'''
       get ENCOD (donor) from a.aliases
        get term_name from biosample_ontology.term_name
        res = requests.get(https://www.encodeproject.org/search/?type=Experiment&biosample_ontology.term_name=spleen&status=released&assay_title=total+RNA-seq&replicates.library.biosample.donor.accession=ENCDO451RUA)
        # list of experiment datasets (c.a. 1) -> ENCSRxxx123 (@graph[0].accession)
        res2 = requests.get(https://www.encodeproject.org/experiments/{accession}
        res.json experiment.files list (output_type = 'gene quantification' and preferred_default=true ) // 1 file
'''

# use the url to get the json for each .bed file



url = 'https://www.encodeproject.org/annotations/ENCSR596QDJ/?format=json'
response = urlopen(url)


cookies = {'required_cookie': 'values'}
headers = {'accept' : 'application/json'}
cookies = {'_gid': 'GA1.1.827781495.1645033198', '_ga': 'GA1.2.1811198692.1634853775', 'session': 'R08tk2RYfGA-eClOXepXPUohq40r62Q-t5AqIOj78ZG1pavyA_isX4zDkoa8Z1DK0l1dgsh3YqjxSYNvFW-EhVsxNjQ1MDUwOTI0LCAxNjQ1MDMzMjAzLjIwOTE3NjMsIHsiX2NzcmZ0XyI6ICJiN2Q3NDZiYTdiNmMzMzUwMGE0Y2EwMDE5Yjg3YzBlNDRkYWFiMWQ3IiwgImF1dGgudXNlcmlkIjogIm1hcmphbi5mYXJhaGJvZEBnbWFpbC5jb20ifV0'} # I need to be logged in, get it from the cookie thing

response = requests.get(url, cookies=cookies, headers=headers)
meta = response.json()
print(meta['aliases'])
print(meta['biosample_ontology']['term_name'])




response = requests.get(wholeList_url, cookies=cookies, headers=headers)
jsonData = response.json()

# I can't do the above cause these pages need me logged in - and I don't want to pass my googel password through. I will work with the data from GCP for this part for now, until I find a solution for that.




