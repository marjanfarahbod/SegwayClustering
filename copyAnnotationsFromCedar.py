# file copy from and to cedar

# just running the scp command for getting annotations from cedar
 
import os
import pickle

# 0. getting segway output runs from cedar
########################################
cedarDataFolder = '/home/mfarahbo/projects/def-maxwl/mfarahbo/segway112Samples/'
localDataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/the38batch/'

# load the address file
inputFile = localDataFolder + 'hg_accessionList.pkl'
with open(inputFile, 'rb') as f:
    accessionList = pickle.load(f)

for accession in accessionList:
    print(accession)
    
    remoteDir = cedarDataFolder + accession + '/segOutput/'
    localDir = localDataFolder + accession + '/'

    # make the folder
    os.mkdir(localDir)

    command = 'scp -r mfarahbo@cedar.computecanada.ca:%s %s' %(remoteDir, localDir)
    print(command)
    os.system(command)

# 1. get the accession track list
########################################

# load the address file
inputFile = localDataFolder + 'hg_accessionList.pkl'
with open(inputFile, 'rb') as f:
    accessionList = pickle.load(f)

for accession in accessionList[15:]:
    print(accession)
    remoteFile = cedarDataFolder + accession + '/trackname_assay.txt'
    localDir = localDataFolder + accession + '/'

    command = 'scp -r mfarahbo@cedar.computecanada.ca:%s %s' %(remoteFile, localDir)
    print(command)
    os.system(command)


# 2. Sending the mnemonics files to cedar
########################################

cedarDataFolder = '/home/mfarahbo/projects/def-maxwl/mfarahbo/segway112Samples/'
localDataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/the38batch/'

# load the address file
inputFile = localDataFolder + 'hg_accessionList.pkl'

inputFile = localDataFolder + 'hg_accessionList.pkl'
with open(inputFile, 'rb') as f:
    accessionList = pickle.load(f)

for accession in accessionList:
    print(accession)
    
    remoteDir = cedarDataFolder + accession + '/'
    localFile = localDataFolder + accession + '/mnemonics_v04.txt'

    command = 'scp -r %s mfarahbo@cedar.computecanada.ca:%s' %(localFile, remoteDir)
    print(command)
    os.system(command)
    



