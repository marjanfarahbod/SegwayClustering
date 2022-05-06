# this is to get the CCRe and chromhmm annotations from the API

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

dataSubFolder = 'testBatch105/fromAPI/'

###################################################
# 1. Get the annotation list (folders exist)
###################################################

# loging to the API and get the cookies

wholeList_url = 'https://www.encodeproject.org/search/?type=Annotation&annotation_type=chromatin+state&lab.title=Maxwell+Libbrecht%2C+SFU&organism.scientific_name=Homo+sapiens&limit=all'
cookies = {'required_cookie': 'values'}
headers = {'accept' : 'application/json'}
cookies = {'_gid': 'GA1.1.1961915786.1649266470', '_ga': 'GA1.1.1847928511.1604510440', 'session': 'JxV9F3PJu9u4DMhNS4ggUf6eACelfke_PP9SiA_gl88FuYq4WbkRiBDdOI7gWbOHKaZqa_jQdejBioirTyEnI1sxNjQ5NDM1NDAyLCAxNjQ5NDM1Mzk4LjY3OTQ1MDMsIHsiX2NzcmZ0XyI6ICJjYWNiMmIxMjlkZDc0YTlkZDBiMDNmNmRlOTBmMjQ0YjQwYzc3ZDRhIiwgImF1dGgudXNlcmlkIjogIm1hcmphbi5mYXJhaGJvZEBnbWFpbC5jb20ifV0'} 
response = requests.get(wholeList_url, cookies=cookies, headers=headers)
wholeList_jdata = response.json()

# 1.2 get the annotation IDs and make folders
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

annotations_IDlist = []
for annotation in wholeList_jdata['@graph']:
    accession = annotation['accession']
    annotations_IDlist.append(accession)

# to figure out the URL for CCRE and ChromHMM, I searched through the annotations page, using the tissue. I don't have the donor ID there (which I need to check at some point)


# from the jsonData, get the list of IDs
i = 0
for accession in annotations_IDlist:

    print('accession count %d' %(i))
    i +=1

    # get the json meta for the annotation
    url = 'https://www.encodeproject.org/annotations/%s/?format=json' %(accession)
    response = requests.get(url, cookies=cookies, headers=headers)
    jsonData = response.json()

    # get the RNA-seq data if present
    aliases = jsonData['aliases'][0]
    index = aliases.index(':')
    
    donorID = aliases[index+1:].split('_')[0]

    term_name = jsonData['biosample_ontology']['term_name']
    term_name_mod = term_name.replace(' ', '+')

    if 'relevant_timepoint' in jsonData.keys():
        donor_age = jsonData['relevant_timepoint']
    else:
        donor_age = 'none'
    
    # I can't search by the donorID , I can only search for the biosample, which is tissue or cell type. We do not have the donor ID for chromhmm and CCRe - ChromHMM has the donor ID in description, but I don't seem to be able to search with that, it is in the description, so pick it from the description - the actual description. Then download and get the file. phew. 

    search_url = 'https://www.encodeproject.org/search/?type=Annotation&biosample_ontology.term_name=%s&status=released&limit=all&annotation_type=candidate+Cis-Regulatory+Elements&annotation_type=chromatin+state' %(term_name_mod)

    response = requests.get(search_url, cookies=cookies, headers=headers)
    search_results = response.json()

    aliases = jsonData['aliases'][0]
    index = aliases.index(':')
    
    donorID = aliases[index+1:].split('_')[0]

    print(len(search_results['@graph']))

    chmm_description = 'ChromHMM'
    CCRe_description = 'Cis-Regulatory'
    reg_description = 'regulatory'

    for result in search_results['@graph']:

        description = result['description']
        if (chmm_description in description):

            if (donorID in description):
                print(description)

                # go to the annotation page
                annot_accession = result['accession']

                annot_url = 'https://www.encodeproject.org/annotations/%s' %(annot_accession)
                response = requests.get(annot_url, cookies=cookies, headers=headers)
                annot_json = response.json()

                annot_fileList = annot_json['files']
                
                donor_age = segway_donor_age

                # get the file with the right version

                for file in annot_fileList:
                    if file['assembly'] == 'GRCh38' and file['file_type'] == 'bed bed9':
                        file_accession = file['accession']

                    
                        downloadURL = 'https://www.encodeproject.org/files/%s/@@download/%s.bed.gz' %(file_accession, file_accession)

                        bedFileName = downloadURL.split('/')[-1]
                        downloadFileName = '/ChromHMM_SegwayDonor_age%s_%s' %(donor_age, bedFileName)
                        downloadPath = dataFolder + dataSubFolder + accession + downloadFileName

                        r = requests.get(downloadURL, cookies=cookies, headers=headers)
                        with open(downloadPath, 'wb') as f:
                            f.write(r.content)

            else:

                # go to the annotation page
                annot_accession = result['accession']
                annot_url = 'https://www.encodeproject.org/annotations/%s' %(annot_accession)
                response = requests.get(annot_url, cookies=cookies, headers=headers)
                annot_json = response.json()

                annot_fileList = annot_json['files']

                if 'relevant_timepoint' in annot_json.keys():
                    donor_age = annot_json['relevant_timepoint']
                else:
                    donor_age = 'none'

                for file in annot_fileList:
                    if file['assembly'] == 'GRCh38' and file['file_type'] == 'bed bed9':
                        file_accession = file['accession']

                        downloadURL = 'https://www.encodeproject.org/files/%s/@@download/%s.bed.gz' %(file_accession, file_accession)

                        bedFileName = downloadURL.split('/')[-1]
                        downloadFileName = '/ChromHMM_age%s_%s' %(donor_age, bedFileName)
                        downloadPath = dataFolder + dataSubFolder + accession + downloadFileName

                        r = requests.get(downloadURL, cookies=cookies, headers=headers)
                        with open(downloadPath, 'wb') as f:
                            f.write(r.content)

                
 

        if (CCRe_description in description):
            print(description)
            annot_accession = result['accession']

            annot_url = 'https://www.encodeproject.org/annotations/%s' %(annot_accession)
            response = requests.get(annot_url, cookies=cookies, headers=headers)
            annot_json = response.json()

            annot_fileList = annot_json['files']

            # get the file with the right version

            if 'relevant_timepoint' in annot_json.keys():
                donor_age = annot_json['relevant_timepoint']
            else:
                donor_age = 'none'


            for file in annot_fileList:
                if file['assembly'] == 'GRCh38' and file['file_format'] == 'bed':
                    file_accession = file['accession']
                    downloadURL = 'https://www.encodeproject.org/files/%s/@@download/%s.bed.gz' %(file_accession, file_accession)

                    bedFileName = downloadURL.split('/')[-1]
                    downloadFileName = '/ccre_age%s_%s' %(donor_age, bedFileName)
                    downloadPath = dataFolder + dataSubFolder + accession + downloadFileName

                    r = requests.get(downloadURL, cookies=cookies, headers=headers)
                    with open(downloadPath, 'wb') as f:
                        f.write(r.content)
                    

        if (reg_description in description):
            print(description)

            annot_accession = result['accession']

            annot_url = 'https://www.encodeproject.org/annotations/%s' %(annot_accession)
            response = requests.get(annot_url, cookies=cookies, headers=headers)
            annot_json = response.json()

            annot_fileList = annot_json['files']

            # get the file with the right version
            if 'relevant_timepoint' in annot_json.keys():
                donor_age = annot_json['relevant_timepoint']
            else:
                donor_age = 'none'

            for file in annot_fileList:
                if file['assembly'] == 'GRCh38' and file['file_format'] == 'bed':
                    file_accession = file['accession']
                    downloadURL = 'https://www.encodeproject.org/files/%s/@@download/%s.bed.gz' %(file_accession, file_accession)

                    bedFileName = downloadURL.split('/')[-1]
                    downloadFileName = '/creg_age%s_%s' %(donor_age, bedFileName)
                    downloadPath = dataFolder + dataSubFolder + accession + downloadFileName

                    r = requests.get(downloadURL, cookies=cookies, headers=headers)
                    with open(downloadPath, 'wb') as f:
                        f.write(r.content)


        # we want the .bed file (not big bed) for the GRCh38 version -

        # it is all in the description:

        
        # if it is chrohmm
        # check the donor ID
        # if the same, get it

        # if it is 

    # from the output, if it is chromhmm, I can see if it is the same donor.

    # if it is the CCRE, there is no way to see the donor, just get it 
    

    #search_url = 'https://www.encodeproject.org/search/?type=Experiment&biosample_ontology.term_name=%s&status=released&assay_title=total+RNA-seq&replicates.library.biosample.donor.accession=%s' %(term_name_mod ,donorID)

    #search_url = 'https://www.encodeproject.org/search/?type=Annotation&biosample_ontology.term_name=%s&status=released&replicates.library.biosample.donor.accession=%s' %(term_name_mod ,donorID)




