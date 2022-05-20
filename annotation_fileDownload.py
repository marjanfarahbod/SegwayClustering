# This is the code to get the annotation files
# Here, we have a file with a list of annotation accession and .tab files where the download link for each of the files exists.

########################################
# 0. Initials 
########################################

import urllib
import linecache
import pickle
import re
import numpy as np
import pandas as pd
import subprocess as sp
import seaborn as sns
import matplotlib.pyplot as plt
from QC_transcriptionComparison_util import Gene, Exon, Annotation, AnnotationClass

# General data folder
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'

dataSubFolder = 'testBatch105/fromAPI/'

inputFile = dataFolder + dataSubFolder + 'metaInfo.pkl'
with open(inputFile, 'rb') as pickledFile:
    annMeta = pickle.load(pickledFile)

sampleFolder_list = list(annMeta.keys())

########################################
# get the probs.txt files from the old files
########################################

fileList = list(glob.iglob(dataFolder + dataSubFolder + 'croo_may_16_2022/*.tsv'))
# get all the .tsv files


# from the html files, get the download links and the workflow_id

# the file with list of samples to identifiers
accession_ID_file = dataFolder + dataSubFolder + 'Jan_segway_runs_Paul - Detailed table.tsv'

identifier_sample_map = {}
with open(accession_ID_file, 'r') as inputFile:
    line = inputFile.readline()
    for line in inputFile:
        fields  = line.strip().split()
        identifier_sample_map[fields[0]] = fields[1]


c = 0
for file in fileList:
    identifier = file.split('.')[-2]
    if identifier in list(identifier_sample_map.keys()):
        accession = identifier_sample_map[identifier]
        downloadFile = dataFolder + dataSubFolder + accession + '/probs.txt'

        # open the file and pars the file for the line with the probs.txt

        with open(file, 'r') as addressFile:
            for line in addressFile:
                gs_link = line.strip().split('\t')[1]
                #print(gs_link[-4:])
                if gs_link.endswith('probs.txt'):
                    url = line.strip().split('\t')[2]
                    # download the probs file to its folder
                    urllib.request.urlretrieve(url, downloadFile)
                    c += 1
                    print(c)
                    print(accession)
                    

########################################
# get the new batch of files 
########################################

# the input folder:
runFolder = dataFolder + 'testBatch_May112022/'

# download file folder
fileList = list(glob.iglob(runFolder + 'croo_may_11_22/*.tsv'))

# Note: these samples are not on the portal yet so they don't have Segway ID, they have the run ID

# make the folders based on the file list:
runID_list = []
id_file_map = {}
for file in fileList:
    run_id = file.split('.')[-2]
    runID_list.append(run_id)
    os.mkdir(runFolder + run_id)
    id_file_map[run_id] = file

file = runFolder + 'runID_list.pkl'
with open(file, 'wb') as outputFile:
    pickle.dump(runID_list, outputFile)


# for each id in the runID_list, download the files listed in the files from id_file_map[id]
c = 0
for run_id in runID_list:

    sampleFolder = runFolder + run_id + '/' # the folder where files are going to be downloaded
    addressesFile = id_file_map[run_id] # the .tab file with addresses of files to be downloaded

    # going to the .tab file with addresses of all files
    with open(addressesFile, 'r') as inputFile:
    #line = inputFile.readline() # the first line is genome data, we don't want that
        print(inputFile)
        print(c)
        c += 1
        downloadFolder = sampleFolder
    
        # make segtools folder
        os.mkdir(downloadFolder + 'call-segtools')
        segtoolsFolder = downloadFolder + 'call-segtools/'
    
        # make interpretation folder
        os.mkdir(downloadFolder + 'call-interpretation')
        interpretationFolder = downloadFolder + 'call-interpretation/'

        # make segway folder
        os.mkdir(downloadFolder + 'call-segway')
        segwayFolder = downloadFolder + 'call-segway/'

        for line in inputFile:
            download_url = line.split('\t')[2]
            fileName = line.split('\t')[1]

            downloadFolder = sampleFolder

            if 'genomedata' in fileName:
                print('gdata - skipping')
                continue
        
            if 'call-segway' in fileName:
                downloadFolder = segwayFolder
                
            if 'call-segtools' in fileName:
                downloadFolder = segtoolsFolder

            if 'call-interpretation' in fileName:
                downloadFolder = interpretationFolder

            downloadFileName = fileName.split('/')[-1]
        
            downloadFile = downloadFolder + downloadFileName
            print(downloadFile)
        
            urllib.request.urlretrieve(download_url, downloadFile)

            print('done download')




