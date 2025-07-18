# checking from the ENCODE api to see which samples are human which are mouse.


from urllib.request import urlopen
import urllib
import json
import http.cookiejar
import requests
import re
import pickle

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
dataSubFolder = 'the112batch/'

##### IMPORTANT: classifier training data
feature_names = ['(09) initial exon', '(01) H3K9me3', '(10) initial intron', '(02) H3K27me3', '(11) internal exons', '(04) H3K4me3', "(16) 3' flanking (1000-10000 bp)", '(12) internal introns', '(03) H3K36me3', '(13) terminal exon', '(06) H3K4me1', '(14) terminal intron', "(07) 5' flanking (1000-10000 bp)", '(05) H3K27ac', "(15) 3' flanking (1-1000 bp)", "(08) 5' flanking (1-1000 bp)"]

# load the address file
inputFile = dataFolder + dataSubFolder + 'accessionList.pkl'
with open(inputFile, 'rb') as f:
    accessionList = pickle.load(f)

for i,accession in enumerate(accessionList[2:15]):

    index = i
    print(accession)

    


