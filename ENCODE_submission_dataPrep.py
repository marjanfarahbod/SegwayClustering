# This is for preparing the data for submission to ENDODE. The first step is to have for each sample, the address of the mnemonics file and the address of the .bed file on Cedar. There are three folders in Cedar that contain segway and segtools outputs for the three sets of runs. I also have a local file with meta data for each of the samples that distinguishes which set each sample belongs too.
#
# Output 1: a dictionary of accession to 1. address of the .bed file on cedar, 2. address of the mnemonics file on cedar
# Output 2: a list of label to color mapping
# The above outputs didn't work as Cedar does not accept some copies. I am going to make mnemonics with the accession and load them somwhere else so Abe can access. 
#
# TODO: just add the accession to each mnemoincs and put them on cedar home folder. 
#
# 0. Initials
# 1. Make the two outputs and copy files to cedar (this didn't work cause cedar does not allow copying of some files)
# 2. cp the mnemnonics files to the folder nemonics with the accession name attached. 
# 3. the speradsheet for the files
########################################
# 0. initials
########################################

import pickle
import os

segwayStates = ['Enhancer', 'EnhancerLow', 'Promoter', 'PromoterFlanking', 'Transcribed', 'CTCF', 'K9K36', 'Bivalent', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent']
print(len(segwayStates))

# General data folder
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'

# # load the metadata local file
inputFileName = 'all235Annot_meta_corrected.pkl'
#inputFileName = 'all235Annot_meta.pkl'
inputFile = dataFolder +  inputFileName
with open(inputFile, 'rb') as f:
    allMeta = pickle.load(f)

accessionList = list(allMeta.keys())

########################################
# 1. make the two outputs and copy files to cedar
########################################

cedarFiles = {}

# for the 105 batch
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
idAccessionFile = dataFolder + 'testBatch105/runID_accession_map_105run.pkl'
with open(idAccessionFile, 'rb') as file:
    idAccessionMap105 = pickle.load(file)

count = 0
for id in idAccessionMap105:
    accession = idAccessionMap105[id]
    print(count)
    count+=1
    
    # the address of the .bed file on cedar
    bedFile = 'projects/rrg-maxwl/mfarahbo/encode-segway/Jan2022/' + id + '/call-segway_annotate/segway.bed.gz'

    # the address of the mnemonics file on cedar
    mnemFile = 'projects/rrg-maxwl/mfarahbo/encode-segway/Jan2022/' + id + '/mnemonics_v04.txt'

    # add the addresses to the dictionary
    cedarFiles[accession] = {'bedFile': bedFile, 'mnemFile':mnemFile}
    
    # >>>>>>>>>>>>>>>>>>>>
    # move the mnemonics file to cedar folder
    localFile = dataFolder + 'testBatch105/fromAPI/' + accession + '/mnemonics_v04.txt'

    remoteDir = '/home/mfarahbo/projects/rrg-maxwl/mfarahbo/encode-segway/Jan2022/' + id 

    command = 'scp -r %s mfarahbo@cedar.computecanada.ca:%s' %(localFile, remoteDir)
    print(command)
    os.system(command)


# for the 95 batch
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

count = 0
for accession in accessionList:
    annotation = allMeta[accession]
    annFolder = annotation['folder']
    print(count)
    count +=1

    # if it is from the 95 runs
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if ('May112' in annFolder) and (count-1 in remainingCounts):
        id = annFolder.split('/')[-2]

        bedFile = 'projects/def-maxwl/mfarahbo/May2022_ENCODEBatch/testBatch_May112022/' + id + '/call-segway/segway.bed.gz'

        mnemFile = 'projects/def-maxwl/mfarahbo/May2022_ENCODEBatch/testBatch_May112022/' + id + '/mnemonics_v04.txt.'

        #cedarFiles[accession] = {'bedFile': bedFile, 'mnemFile':mnemFile}

        # >>>>>>>>>>>>>>>>>>>>
        # move the mnemonics file to cedar folder
        localFile = annFolder + 'mnemonics_v04.txt'

        remoteDir = 'projects/def-maxwl/mfarahbo/May2022_ENCODEBatch/testBatch_May112022/' + id

        #remoteFile = 'mfarahbo/remainingStuff/' + accession +' _mnemonics_v04.txt'
        command = 'scp -r %s mfarahbo@cedar.computecanada.ca:%s' %(localFile, remoteDir)
        print(command)
        os.system(command)

# counts that didn't work for this one:
remainingCounts = [154, 155, 156, 161, 166, 167, 171,172, 176, 179, 181,182, 184, 185,186, 189,190, 194, 195196, 198, 201, 202, 203, 204, 209, 210, 211, 212, 214, 215, 217, 220, 224, 225, 226, 227, 230, 231, 232, 233, 234]

# for the 38 batch
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

count = 0
for accession in accessionList:
    annotation = allMeta[accession]
    annFolder = annotation['folder']
    print(count)
    count +=1

    # if it is from the 95 runs
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if ('38batch' in annFolder):
        bedFile = '/home/mfarahbo/projects/def-maxwl/mfarahbo/segway112Samples/' + accession + '/segOutput/segway.bed'

        mnemFile = '/home/mfarahbo/projects/def-maxwl/mfarahbo/segway112Samples/' + accession + '/mnemonics_v04.txt.'

        cedarFiles[accession] = {'bedFile': bedFile, 'mnemFile':mnemFile}

        # >>>>>>>>>>>>>>>>>>>>
        # move the mnemonics file to cedar folder
        localFile = annFolder + 'mnemonics_v04.txt'

        remoteDir = '/home/mfarahbo/projects/def-maxwl/mfarahbo/segway112Samples/' + accession

        command = 'scp -r %s mfarahbo@cedar.computecanada.ca:%s' %(localFile, remoteDir)
        print(command)
        os.system(command)


outputFile = dataFolder + 'cedar_bed_mnem_addresses_accessionDict.pkl'
with open(outputFile, 'wb') as f:
    pickle.dump(cedarFiles, f)

        
########################################
# 2. cp the mnemnonics files to the folder nemonics with the accession name attached.
########################################

destinationFolder = dataFolder + 'allMnemonics_v04/'
for accession in accessionList:
    annotation = allMeta[accession]
    annFolder = annotation['folder']

    mnemFile = annFolder + 'mnemonics_v04.txt'
    destMnemFile = destinationFolder + accession + '_mnemonics_v04.txt'

    command = 'cp %s %s' %(mnemFile, destMnemFile)
    os.system(command)

########################################
# 3. the spreadsheet for the files, sample accession to the list of track accessions
########################################

count = 0
trackAccessions = {}
for accession in accessionList:
    annotation = allMeta[accession]
    annFolder = annotation['folder']
    count+=1

    trackFile = annFolder + 'trackname_assay.txt'

    # read file by line, save the accession list with comma
    accessionStr = ''
    with open(trackFile, 'r') as file:
        for line in file:
            trackAccession = line.split()[0]
            accessionStr = accessionStr + trackAccession + ', '

        accessionStr = accessionStr[0:-2]

    trackAccessions[accession] = accessionStr


outputFile = dataFolder + 'accession_trackAccession.tsv'
import csv
with open(outputFile, 'w') as f:
    for accession in accessionList:
        f.write('%s\t%s\n' %(accession, trackAccessions[accession]))
            
        



