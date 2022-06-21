# To get the transcript evaluation for Mehdi's work
###################################################
# 0. Initials
###################################################

from urllib.request import urlopen
import urllib
import json
import http.cookiejar
import requests
import re

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/Mehdi_testRun/tomarjan/'

###################################################
# 0. Initials
###################################################


# from each folder, get the tissue name, get the transcript data from the portal for that tissue

# get the list of folders:

os.listdir(dataFolder + 'segway')
