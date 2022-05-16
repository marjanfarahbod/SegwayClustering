# for comparing the annotations with CCRE
#
# ITEMS IN THE CODE:
# ########################################
# 0. Initials
# 1. pre-process the CCRE
# 2. Get the summary info for each ccre file
#
########################################
# 0. Initials 
########################################

import linecache
import pickle
import re
import numpy as np
import pandas as pd
import subprocess as sp
import seaborn as sns
import matplotlib.pyplot as plt
import glob
from QC_transcriptionComparison_util import Gene, Exon, Annotation, AnnotationClass

# General data folder
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'

dataSubFolder = 'testBatch105/fromAPI/'

# list of annotation folders
inputFile = dataFolder + dataSubFolder + 'metaInfo.pkl'
with open(inputFile, 'rb') as pickledFile:
    annMeta = pickle.load(pickledFile)

sampleFolder_list = list(annMeta.keys())

########################################
# 1. pre-process the CCRE 
########################################

count = 0
ccreFile_dict = {}
for sampleFolder in sampleFolder_list:

    print(count)
    print(sampleFolder)
    count +=1

    sampleFolder_add = dataFolder + dataSubFolder + sampleFolder + '/'

    '''
    get the list of files in the folder, pick ones that starts with ccre_*

    '''

    fileList = list(glob.iglob(sampleFolder_add + 'ccre_*', recursive = True))
    if len(fileList) > 0:
        ccreFile = fileList[0]

        # unzip the file
        if ccreFile.endswith('.gz'):
            os.system('gunzip %s' %(ccreFile))
            ccreFile = ccreFile[:-3]

        # sort the chromosome column
        fileName = ccreFile.split('/')[-1]
        sortedCcreFile = sampleFolder_add + 'sorted_' + fileName
        os.system('sort -V -k1,1 -k2,2 %s > %s' %(ccreFile, sortedCcreFile))
        ccreFile_dict[sampleFolder] = sortedCcreFile
    else:
        print('none')
        ccreFile_dict[sampleFolder] = 'none'
        
outputFile = dataFolder + dataSubFolder + 'ccre_list_dict.pkl'
with open(outputFile, 'wb') as f:
    pickle.dump(ccreFile_dict, f)

########################################
# 2. Get general info about CCRE files
########################################

inputFile = dataFolder + dataSubFolder + 'ccre_list_dict.pkl'
with open(outputFile, 'rb') as f:
    ccreFile_dict = pickle.load(f)

sampleFolder_list = list(ccreFile_dict.keys())

# get an idea of the count of ccre or ccreg files:
for sample in sampleFolder_list:
    sampleFolder_add = dataFolder + dataSubFolder + sample + '/'

    fileList = list(glob.iglob(sampleFolder_add + 'ccre_*', recursive = True))
    print(len(fileList))
    fileList = list(glob.iglob(sampleFolder_add + 'creg_*', recursive = True))
    print(len(fileList))

#make the plot folders
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/testBatch105/'
for sample in sampleFolder_list:
    plotFolder_add = plotFolder + sample + '/'
    os.mkdir(plotFolder_add)

# get the boxplot of lenght distribution for each class in ccre for each plot
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/testBatch105/'
for sample in sampleFolder_list:
    sampleFolder_add = dataFolder + dataSubFolder + sample + '/'

    plotFolder_add = plotFolder + sample + '/'

    # for ccre files
    fileList = list(glob.iglob(sampleFolder_add + 'ccre_*', recursive = True))
    for ccreFile in fileList:
        if ccreFile.endswith('.gz'):
            os.system('gunzip %s' %(ccreFile))
            ccreFile = ccreFile[:-3]

        ccreAccession = re.search('ENCFF.*\.bed', ccreFile)[0][:-4]
        print(ccreAccession)

        # do the plot
        
        labels_dict = {}
        max_bp_count = 0
        min_bp_count = 1000
        with open(ccreFile, 'r') as annotations:

            for line in annotations:
                fields = line.strip().split()
                length = int(fields[2]) - int(fields[1])
                label = fields[9]
                if label in list(labels_dict.keys()):
                    labels_dict[label].append(length)
                else:
                    labels_dict[label] = list()
                    labels_dict[label].append(length)

                max_bp_count = max(length, max_bp_count)
                min_bp_count = min(length, min_bp_count)

        labels_list = list(labels_dict.keys())

        # for each label, append the list
        data = list()
        for label in labels_list:
            data.append(labels_dict[label])

        label_bpCount = np.zeros(len(labels_list))
        for i in range(len(labels_list)):
            label_bpCount[i] = sum(data[i])

        #ccreFile_name = ccreFile.split('/')[-1][0:-4]
        #plotFile_name = 'boxplot_%s.pdf' %(ccreFile_name)
        #figFile = plotFolder_add + plotFile_name
        #print(figFile)
        
        # just the boxplot

        label_bpCount_frac = label_bpCount/(3e9)
        label_bpCount_frac_rounded = [round(x, 4) for x in label_bpCount_frac]

        numBoxes = len(labels_list)
        fig, ax1 = plt.subplots(figsize=(10,7))
        # fig.canvas.set_window_title('ccre label distribution')
        plt.subplots_adjust(left=.075, right=.95, top=.9, bottom=.25)

        ax1.set_ylabel('base pair count')
        ax1.set_title('ccre labels in genome - fraction of the base pairs on top')

        bp = plt.boxplot(data, positions = np.arange(1, len(labels_list) + 1, step=1))
        plt.setp(bp['boxes'], color='black')

        # set axes ranges and labels
        ax1.set_xlim(.5, numBoxes + .5)
        xtickNames = plt.setp(ax1, xticklabels=labels_list)
        plt.setp(xtickNames, rotation=45, fontsize=10)
        top = max_bp_count + 10
        bottom = min_bp_count - 10
        ax1.set_ylim(bottom, top)

        # adding the top xtick labels
        pos = np.arange(numBoxes) + 1
        upperLabels = label_bpCount_frac_rounded

        for tick, label in zip(range(numBoxes), ax1.get_xticklabels()):
            k = tick % 2
            ax1.text(pos[tick], top -5, upperLabels[tick],
                     horizontalalignment='center', size='small')

        #plt.show()

        # save the plot
        ccreFile_name = ccreFile.split('/')[-1][0:-4]
        plotFile_name = 'boxplot_%s.pdf' %(ccreFile_name)
        figFile = plotFolder_add + plotFile_name
        print(figFile)
        plt.savefig(figFile)
        plt.close('all')

#########################################
# 3. The comparison
#########################################




#########################################
# DRAFT
#########################################

        plt.savefig(figFile)
        plt.close('all')
        
        # >>>>> plot draft
        ax1.set_xlabel('xlabel')
        ax1.set_ylabel('base pair count')


        ax1.set_xticklabels(labels_list, rotation = 45, ha = 'right')
        #ax1.set_xticks(np.arange(1, len(labels_list) + 1, step =1))

        ax2 = ax1.twiny() # instantiate a second axes that shares the same x-axis

        top_axis = np.arange(1, len(labels_list) + 1, step=1)
        top_axis = np.append(top_axis, len(roundedBook))

        ax2.set_xticks(top_axis)
        ax2.set_xticklabels(label_bpCount_frac_rounded)

        
        plt.show()
        
        top_axis = np.arange(.5, len(labels_list) + .5, step=1)
        top_axis = np.append(top_axis, len(roundedBook))
        ax1.set_xticks(top_axis)
        ax1.set_xticklabels(label_bpCount_frac_rounded)

        ax1.tick_params(axis = 'y')

        
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        fig = plt.figure(figsize =(7, 10))
        
        ax = fig.add_subplot(111)
        bp = ax.boxplot(data, positions = np.arange(1, len(labels_list) + 1, step=1))
        #bp = ax.boxplot(data)
        ax.set_xticklabels(labels_list, rotation = 45, ha = 'right')
        ax.set_xticks(np.arange(1, len(labels_list) + 1, step =1))

        axT = ax.twinx()
        top_axis = np.arange(.5, len(labels_list) + .5, step=1)
        top_axis = np.append(top_axis, len(roundedBook))
        axT.set_xticks(top_axis)
        axT.set_xticklabels(label_bpCount_frac_rounded)
        plt.ylabel('base pair count')
        plt.xlabel('fraction of the genome')

        title = 'ccre accession %s' %(ccreAccession)
        plt.title(title)
        plt.gcf().subplots_adjust(bottom=0.30)
        plt.show()

    # for creg files
    fileList = list(glob.iglob(sampleFolder_add + 'creg_*', recursive = True))
    for cregFile in fileList:

        # if zip, unzip

        # get the accession

        

    ccreFile = ccreFile_dict[sample]

    labels_dict = {}
    with open(ccreFile, 'r') as annotations:

        # annotations have no header
        # header = annotations.readline()

#        f = open(bedFileAdd, 'r') # >>>> test
#       line = f.readline() # >>>> test

        for line in annotations:
            fields = line.strip().split()
            length = int(fields[2]) - int(fields[1])
            label = fields[9]
            if label in list(labels_dict.keys()):
                labels_dict[label].append(length)
            else:
                labels_dict[label] = list()
                labels_dict[label].append(length)

        labels_list = list(labels_dict.keys())
        data = [labels_dict[labels_list[0]], labels_dict[labels_list[1]], labels_dict[labels_list[2]], labels_dict[labels_list[3]], labels_dict[labels_list[4]], labels_dict[labels_list[5]], labels_dict[labels_list[6]]]
        data = 

        fig = plt.figure(figsize =(7, 7))
        
        #ax = fig.add_axes([0, 0, 1, 1])
        ax = fig.add_subplot(111)
    
        
        bp = ax.boxplot(data, positions = np.arange(1, len(roundedBook) + 1, step=1))
        #bp = ax.boxplot(data)
        # hp = ax.hist(labels_dict[labels_list[0]])
        ax.set_xticklabels(labels_list, rotation = 45, ha = 'right')
        ax.set_xticks(np.arange(1, len(roundedBook) + 1, step =1))

        axT = ax.twiny()
        book = label_bpCount/(3e9)
        roundedBook = [round(x, 4) for x in book]
        kado = np.arange(.5, len(roundedBook) + .5, step=1)
        kado = np.append(kado, len(roundedBook))
        axT.set_xticks(kado)
        axT.set_xticklabels(roundedBook)
        plt.ylabel('base pair count')
        plt.xlabel('fraction of the genome')


        title = 'ccre accession %s' %(ccreAccession)
        plt.title(title)

        plt.gcf().subplots_adjust(bottom=0.30)

        plt.show()

# TODO: put the code in, scale and save the files!
# TODO: comparison to the chrom hmm and do the plots! 

# just add the bp count for each section
# do the plot perhaps for all of them and also the creg
        
########################################
# Get the summary info for each ccre file
########################################

inputFile = dataFolder + dataSubFolder + 'ccre_list_dict.pkl'
with open(outputFile, 'rb') as f:
    ccreFile_dict = pickle.load(f)

sampleFolder_list = list(ccreFile_dict.keys())

for sample in sampleFolder_list:

    ccreFile = ccreFile_dict[sample]

    if ccreFile != 'none':
        print(ccreFile)
        f = open(ccreFile, 'r') # >>>> test
        line = f.readline() # >>>> test

        fields = line.strip().split()
        print(fields[-1])
    


