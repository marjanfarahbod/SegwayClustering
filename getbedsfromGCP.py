# from a file with addresses, download all the .bed files.
# 0. General initiation
# 1. Open the file and read it into a dictionary with fields
# 2. Filtering and download setting # 3. Download the files
# 4. Unzip and merge
# #. DRAFT

# 0. general initiation
#####################
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/testBatch/'
gsc = 'gs://segway-testing/segway/'

# 1. Open the file and read it into a dictionary with fields
#####################
# input: the address of the ID/name/GCP link table file
# output: the dictionary where I have info for each sample
# the input is the common file format I get from Paul for different runs, however for the October runs I just got list of IDs and addresses and had to match IDs to previous runs. The dir changes for different runs while the name and ID (sample/experiment ID) remains the same

# the ID table (input)
metaTableFileName = 'Segway_runs_imputation_challenge_samples_6_1runs.tsv'
metaFile = dataFolder + metaTableFileName

# reading lines and adding items to a list
metaData = []
with open(metaFile, 'r') as input:
    next(input) # skip the header
    for line in input:
        fields = line.strip().split('\t')
        if len(fields) > 1: # last lines are blank and a note
            exp = {"name": fields[0], "ID": fields[1], "dir": fields[2]}
            metaData.append(exp)

# reading the lines and getting a name/id dictionary    
nameIDDict = {}
with open(metaFile, 'r') as input:
    next(input)
    for line in input:
        fields = line.strip().split('\t')
        if len(fields) > 1: # last lines are blank and a note
            nameIDDict[fields[1]] = fields[0]


# 2. Filtering and download setting # 3. Download the files
#####################
#####################
# 2. This depends on the project stage. For the Oct2021 runs, I had to match the experiment ID to dir and name
# 3. for each itme, make the folder and download the file
# TODO: 

# oct2021 runs
filterList = os.listdir(dataFolder)
infoFile = dataFolder + 'oct2021_IDs.tsv'
experimentFolderNames = [] # list of experiment names (we might have some folders from before or the folder names might change)
with open(infoFile, 'r') as file:
    next(file)
    for line in file:
        fields = line.strip().split('\t')
        name = nameIDDict[fields[1]]

        # downloading the files
        if name in filterList:
            experimentFolderNames.append(name)
            command1 = '/Users/marjanfarahbod/./google-cloud-sdk/bin/gsutil cp -r ' + fields[0][0:-1] + '/call-recolor_bed/' + 'recolored.bed.gz ' + dataFolder + name + '/'
            command2 = '/Users/marjanfarahbod/./google-cloud-sdk/bin/gsutil cp -r ' + fields[0][0:-1] + '/call-segway_annotate/' + 'segway.bed.gz ' + dataFolder + name + '/'
            os.system(command1)
            os.system(command2)
        else: # else we need to make the folder for it, after we remove the space and adjust the name for the command
            name = name.replace(' ', '')
            name = name.replace('(', r"\(")
            name = name.replace(')', r"\)")
            experimentFolderNames.append(name)
            os.mkdir(dataFolder + name)
            command1 = '/Users/marjanfarahbod/./google-cloud-sdk/bin/gsutil cp -r ' + fields[0][0:-1] + '/call-recolor_bed/' + 'recolored.bed.gz ' + dataFolder + name + '/'
            command2 = '/Users/marjanfarahbod/./google-cloud-sdk/bin/gsutil cp -r ' + fields[0][0:-1] + '/call-segway_annotate/' + 'segway.bed.gz ' + dataFolder + name + '/'
            os.system(command1)
            os.system(command2)

# 4. unzip and merge
####################
# for all the folders in the data directory do this
# TODO: add catch exception for cases where files are already unzipped or do not exist

import os.path
for item in experimentFolderNames:
    path = dataFolder + item
    if os.path.isdir(path):
        print(path)
        name = item
        #unzipping the files
#        name = name.replace('(', r"\(")
#       name = name.replace(')', r"\)")
        command = 'gunzip ' + dataFolder + name + '/recolored.bed.gz'
        os.system(command)
        command = 'gunzip ' + dataFolder + name + '/segway.bed.gz'
        os.system(command)

# merging
# remove the header and merge the column. The two files have same count of lines and the lines match
for expFolder in experimentFolderNames:
    file1 =  dataFolder + expFolder + '/recolored.bed' 
    file2 =  dataFolder + expFolder + '/segway.bed' 
    output = dataFolder + expFolder + '/ccSegway.bed'
    file1 = file1.replace('\\', '')
    file2 = file2.replace('\\', '')
    output = output.replace('\\', '')
    print(output)
    with open(file1, 'r') as classFile, open(file2, 'r') as clusterFile, open(output, 'w') as mergedFile:
        next(classFile)
        next(clusterFile)
        for line in classFile:
            classLine = line.strip().split()
            clusterLine = clusterFile.readline().strip().split()
            mergedFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(classLine[0], classLine[1], classLine[2], clusterLine[3] + '_' + classLine[3], classLine[4], classLine[5], classLine[6], classLine[7], classLine[8]))


# >>>>>>>>>>>>>>>>>>>> DRAFT
# D2: draft from task 2
#  below is not used currently, it is fetch the bed files from info in the meta file 
# get the list of folders in the data file, for folders in the meta file do the task
filterList = os.listdir(dataFolder)
for item in metaData:
    if item['name'] in filterList: # filtering 
        command1 = '/Users/marjanfarahbod/./google-cloud-sdk/bin/gsutil cp -r gs://segway-outputs/segway/' + item['dir'] + '/call-recolor_bed/' + 'recolored.bed.gz ' + dataFolder + item['name'] + '/'
        os.system(command1)


