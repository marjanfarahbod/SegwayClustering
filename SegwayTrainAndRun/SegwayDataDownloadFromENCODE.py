# This is a main file. Given is the input file, it usese the function fileDownloadFromENCODE.py

import ast
import urllib
import json
import http.cookiejar
import requests
import re
import os


### local add
#dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/the112batch/'

### cedar add
dataFolder = '/home/mfarahbo/scratch/'

sampleFileLists = 'SegwayFileSets.tsv'
# for each sample, make the folder with the accession name
inputFile = dataFolder + sampleFileLists
with open(inputFile, 'r') as f:
    header = f.readline()
    for line in f:
        fields = line.strip().split('\t')
        print(line)

        # with fields 0 make the folder
        segwayAccession = fields[0]
        print('now in sample %s' %(segwayAccession))
        
        # make the folder
        newFolder = dataFolder + segwayAccession
        os.mkdir(newFolder)

        # with fields 1 make the trackname_assay.txt and download
        tracks = ast.literal_eval(fields[1])
        trackList = list(tracks.keys())

        # download the files
        for track in tracks:
            dladd = tracks[track][1]
            fileAccession = dladd.split('/')[-1]
            downloadURL = 'https://www.encodeproject.org/' + dladd
            downloadPath = dataFolder + segwayAccession + '/' + fileAccession
            urllib.request.urlretrieve(downloadURL, downloadPath)
            print('file %s downloaded' %(track))

            
        # while at it, make the track file
        tnaFile = dataFolder + segwayAccession + '/trackname_assay.txt'
        with open(tnaFile, 'w') as tna:
            for track in tracks:
                fileAccession = tracks[track][1].split('/')[-1][0:-7]
                print(fileAccession)
                tna.write('%s\t%s\n' %(fileAccession, track))
                
            
