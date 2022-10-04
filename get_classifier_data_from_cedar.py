# get the list of files
# cp ones that start with [] and end with [] to the scratch with the folders.
# scp from scratch to the desktop
# manage on the desktop

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'

# whole file list
allFiles_file = dataFolder + 'fileList.txt'

# classifier data files
class_files = dataFolder + 'classifier_data_file_list.txt'

c = 0
with open(allFiles_file, 'r') as inputFile, open(class_files, 'w') as outputFile:
    for line in inputFile:
        if line.strip().endswith('classifier_data.tab'):
            local_add = line[30:]
            outputFile.write(local_add)


