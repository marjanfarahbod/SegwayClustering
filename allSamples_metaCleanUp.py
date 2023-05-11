# Collect list of samples based on accession. For each sample, record the address of the transcriptomic file, if exists. Record the address of chromhmm file, if exists.

import util # this is for features_from_segtools_dir
import gzip
import pickle
import pandas as pd

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
# the file with most metadata
file = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/MetaSheet310Samples.tsv'
accession_meta = {} # just mapping the accession to the folder

# I have three batches of Segway samples. The 105JanBatch, the 112HomeBatch (actually 38 samples), the 95MayBatch. Total of human samples: 238

# for each file, I want to get the folder address for now, and which address it belongs to. 105, 95, or 38
df = pd.read_table(file)
for i,header in enumerate(list(df)):
    print('%d %s' %(i, header))

#  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# the 38 runs
dataSubFolder = 'the38batch/'

# load the address file
inputFile = dataFolder + dataSubFolder + 'hg_accessionList.pkl'
with open(inputFile, 'rb') as f:
    accessionList = pickle.load(f)

for accession in accessionList:
    sampleFolder = dataFolder + dataSubFolder + accession + '/'
    accession_meta[accession] = sampleFolder 

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# the 105 runs
dataSubFolder = 'testBatch105/fromAPI/'

inputFile = dataFolder + dataSubFolder + 'metaInfo.pkl'
with open(inputFile, 'rb') as f:
    annMeta = pickle.load(f)

accessionList = list(annMeta.keys())
for accession in accessionList:
    sampleFolder = dataFolder + dataSubFolder + accession + '/'
    accession_meta[accession] = sampleFolder

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# the 92 runs
dataSubFolder = 'testBatch_May112022/'

# get the runID with the accession (cause that's where the folder is)
metaInfoFile = dataFolder + dataSubFolder + 'runID_accession.tsv'

df = pd.read_table(metaInfoFile)
for i,header in enumerate(list(df)):
    print('%d %s' %(i, header))

for i, runID in enumerate(df['cromwell ID']):
    sampleFolder = dataFolder + dataSubFolder + runID + '/'
    accession = df['Annotation on portal'][i].split('/')[2]
    accession_meta[accession] = sampleFolder + '/'

########################################
# building the annotation meta list (with API)
########################################


from urllib.request import urlopen
import urllib
import json
import http.cookiejar
import requests
import re

cookies = {'required_cookie': 'values'}
headers = {'accept' : 'application/json'}
cookies = {'_gid': 'GA1.1.455974259.1683650576', '_ga': 'GA1.1.1807347774.1675726254', 'session': 'AQGy9-buf509Y1nLBAsy4FnrBWkltjU8v2DNcR8P9C4VXpIlnrKC81EBW2nFv3KOzxINqj8BcNjr1PN-Bn0T8VsxNjgzNjUwNTgyLCAxNjgzNjUwNTc4LjI4MjY3NzIsIHsiX2NzcmZ0XyI6ICIyYTBiZjA3N2YzYjgxYjNmMWY0OTk5ZGNkNTM3MTdlODA0ZjIyYWMyIiwgImF1dGgudXNlcmlkIjogIm1hcmphbi5mYXJhaGJvZEBnbWFpbC5jb20ifV0'} 

allAccession = list(accession_meta.keys())

# for the 38
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# get the tissue
# load the list and folder address, get the transcript and chromhmm, save the data

#>>>>>>>>>> TODO: add the address of the rnaseq files to the meta info. add the address of the tissue to the meta info
accession38_92RNAseq = {}
accession38_92TissueInfo = {}
accession38_92Chmm = {}

inds = list(range(0,38)) + list( range(143, len(allAccession)))

for i in inds:

    accession = allAccession[i]
    print(accession)
    print(i)
    rnaSeqFiles = []
    annFolder = accession_meta[accession]

    # get the tissue from API
    url = 'https://www.encodeproject.org/annotations/%s/?format=json' %(accession)
    response = requests.get(url, cookies=cookies, headers=headers)
    jsonData = response.json()

    # get the donor info
    aliases = jsonData['aliases'][0]
    index = aliases.index(':')
    donorID = aliases[index+1:].split('_')[0]

    # get the tissue info: name, whole name and ID, type
    termName_id_list = [jsonData['biosample_ontology']['term_name'], jsonData['biosample_ontology']['@id'], jsonData['biosample_ontology']['classification']]

    accession38_92TissueInfo[accession] = [termName_id_list, donorID]

    # get the transcriptomic data
    term_name = jsonData['biosample_ontology']['term_name']
    term_name_mod = term_name.replace(' ', '+')

    search_url = 'https://www.encodeproject.org/search/?type=Experiment&biosample_ontology.term_name=%s&status=released&assay_title=total+RNA-seq&replicates.library.biosample.donor.accession=%s' %(term_name_mod ,donorID)

    response = requests.get(search_url, cookies=cookies, headers=headers)
    search_results = response.json()

    if response.status_code == '404':
        continue
    else:
        c = np.min([len(search_results['@graph']), 2])
        for result in search_results['@graph'][0:c]:
            RNA_accession = result['accession']
            print('annotation accession ' + accession)
            print(RNA_accession)
        
            # go to the RNA page and download the gene quantity files in the folder
#            rna_url = 'https://www.encodeproject.org/experiments/ENCSR052FJA/'
            rna_url = 'https://www.encodeproject.org/experiments/%s' %(RNA_accession)
            response = requests.get(rna_url, cookies=cookies, headers=headers)
            rna_json = response.json()

            rna_fileList = rna_json['files']

            for file in rna_fileList:
                # if it is gene quantifaction, preferred default and .tsv
                if file['output_type'] == 'gene quantifications' and file['file_format'] == 'tsv' and 'preferred_default' in file:
                    print(file['accession'])

                    file_accession = file['accession']
                    downloadURL = 'https://www.encodeproject.org/files/%s/@@download/%s.tsv' %(file_accession, file_accession)
                    
                    downloadPath = annFolder + 'preferred_default_' + file_accession + '.tsv'
                        
                    print('doing the download')
                    urllib.request.urlretrieve(downloadURL, downloadPath)

                    rnaSeqFiles.append(downloadPath)

    accession38_92RNAseq[accession] = rnaSeqFiles

# get the chrom data


for i in inds:

    accession = allAccession[i]
    print(accession)
    print(i)

    annotation_folder = accession_meta[accession]
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

    gotFile = False
    for result in search_results['@graph']:

        description = result['description']
        if (chmm_description in description and gotFile == False):# and not('treated' in description)):

            if (donorID in description and not('treated' in description)):
                print(description)

                # go to the annotation page
                annot_accession = result['accession']

                annot_url = 'https://www.encodeproject.org/annotations/%s' %(annot_accession)
                response = requests.get(annot_url, cookies=cookies, headers=headers)
                annot_json = response.json()

                annot_fileList = annot_json['files']
                
                # get the file with the right version

                for file in annot_fileList:
                    if file['assembly'] == 'GRCh38' and file['file_type'] == 'bed bed9':
                        file_accession = file['accession']

                    
                        downloadURL = 'https://www.encodeproject.org/files/%s/@@download/%s.bed.gz' %(file_accession, file_accession)

                        bedFileName = downloadURL.split('/')[-1]
                        downloadFileName = '/ChromHMM_SegwayDonor_age%s_%s' %(donor_age, bedFileName)
                        downloadPath = annotation_folder + downloadFileName

                        r = requests.get(downloadURL, cookies=cookies, headers=headers)
                        with open(downloadPath, 'wb') as f:
                            f.write(r.content)

                        accession38_92Chmm[accession] = downloadPath
                        gotFile = True

            else:

                print('not same donor')
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
                        downloadPath = annotation_folder + downloadFileName

                        r = requests.get(downloadURL, cookies=cookies, headers=headers)
                        with open(downloadPath, 'wb') as f:
                            f.write(r.content)

                        accession38_92Chmm[accession] = downloadPath
                        gotFile = True


accession38_92RNAseq # where the RNAseq file is
accession38_92TissueInfo # all tissue info
accession38_92Chmm # where the chrom file is
accession_meta # the address of the annotation info

rawDataMerge = [accession38_92RNAseq, accession38_92TissueInfo, accession38_92Chmm, accession_meta]
outputFileName = 'rawData_meta_38_92.pkl'
outputFile = dataFolder +  outputFileName
with open(outputFile, 'wb') as f:
    pickle.dump(rawDataMerge, f)

                        
# for the 105
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# General data folder
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
dataSubFolder = 'testBatch105/fromAPI/'

# load the meta and accession, fill the info based on that
inputFile = dataFolder + dataSubFolder + 'all_annInfo.pkl'
with open(inputFile, 'rb') as f:
    annMeta = pickle.load(f)

inputFile = dataFolder + dataSubFolder + 'metaInfo.pkl'
with open(inputFile, 'rb') as f:
    annMeta105BedRNA = pickle.load(f)

inputFile = dataFolder + dataSubFolder + 'segwayAccession_tissueInfo.pkl'
with open(inputFile, 'rb') as f:
    tissueInfo = pickle.load(f)

inputFile = dataFolder + dataSubFolder + 'chmmFile_list_dict.pkl'
with open(inputFile, 'rb') as f:
    chmmFile_dict = pickle.load(f)
    
accession105 = list(annMeta.keys())
accessionOthers = list(accession_meta.keys())

rnaseqAccession = list(accession38_92RNAseq.keys()) # where the RNAseq file is
accession38_92TissueInfo # all tissue info
chromAccession = list(accession38_92Chmm.keys()) # where the chrom file is

# for accession in accession_meta, if it is in 105, fetch data based on that.
allAccessionMeta = {}
for i, accession in enumerate(accession_meta):
    print(i)
    print(accession)
    thisAccession = {}
    if accession in accession105:
        
        annFolder = dataFolder + 'testBatch105/fromAPI/' + accession + '/'
        
        thisAccession['accession'] = accession
        thisAccession['folder'] = annFolder
        
        originalBedFile = annMeta105BedRNA[accession]['bedFile']
        bedAccession = originalBedFile.split('.')[0]
        thisAccession['bedFile'] = annFolder + '%s_filteredSorted.bed' %(bedAccession)
        
        thisAccession['tissueInfo'] = tissueInfo[accession]
        thisAccession['donorInfo'] = 'none'
        thisAccession['chromFile'] = chmmFile_dict[accession]
        
        if annMeta105BedRNA[accession]['expressionFile'] == 'none':
            thisAccession['RNAseqFile'] = 'none'
        else:
            thisAccession['RNAseqFile'] = annFolder + annMeta105BedRNA[accession]['expressionFile']

    else:

        annFolder = accession_meta[accession]
        thisAccession['accession'] = accession
        thisAccession['folder'] = annFolder

        thisAccession['bedFile'] = annFolder + 'bedFile_filteredSorted.bed' # to be made
        
        thisAccession['tissueInfo'] = accession38_92TissueInfo[accession][0]
        thisAccession['donorInfo'] = accession38_92TissueInfo[accession][1]

        if accession in chromAccession:
            # fill up
            thisAccession['chromFile'] = accession38_92Chmm[accession]
        else:
            thisAccession['chromFile'] = 'none'

        if accession in rnaseqAccession:
            if len(accession38_92RNAseq[accession]) == 0:
                thisAccession['RNAseqFile'] = 'none'
            else:
                thisAccession['RNAseqFile'] = accession38_92RNAseq[accession]
        else:
            thisAccession['RNAseqFile'] = 'none'

    allAccessionMeta[accession] = thisAccession

outputFileName = 'all235Annot_meta.pkl'
outputFile = dataFolder +  outputFileName
with open(outputFile, 'wb') as f:
    pickle.dump(allAccessionMeta, f)

########################################
# how many samples has transcriptomic

rnaseqCount = 0
for annot in allAccessionMeta:
    if not(allAccessionMeta[annot]['RNAseqFile'] == 'none'): 
        rnaseqCount+=1

# we have 94 samples with RNAseq
