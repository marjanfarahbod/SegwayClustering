# Plots for paper Figure 1
# 0. Initials
# 1. The heatmap for assay presence DONE
# 2. The actual annotation for samples sorted - I donno what this is
# 3. The long heatmap with all things sorted for all the signal. For genes without transcriptomic data, just the genomic regions(supplement). (so I want each "interpretation term" to be clustered, then them presented in a long heatmap)
# 4. The prediction probabilities DONE
# *. the big fat genome regional plot
# *. The prediction some sort of filtering?
# *. The other tracks bar thing plot
#
#
#########################################
# 0. Initials
#########################################
from numpy import random
import pickle
import re
import numpy as np
import pandas as pd
import subprocess as sp
import seaborn as sns
import matplotlib.pyplot as plt

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/'

# load the meta file 
inputFileName = 'all235Annot_meta_corrected.pkl'
#inputFileName = 'all235Annot_meta.pkl'
inputFile = dataFolder +  inputFileName
with open(inputFile, 'rb') as f:
    allMeta = pickle.load(f)

accessionList = list(allMeta.keys())
accessionList.remove(accessionList[205])
# Segway states:
segwayStates = ['Promoter', 'PromoterFlanking', 'Enhancer', 'EnhancerLow', 'Bivalent', 'CTCF', 'Transcribed', 'K9K36', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent']
print(len(segwayStates))
segwayStateCount = len(segwayStates)

########################################
# 1. get the sample first quartile probabilities
########################################

#pMetricList = np.zeros(len(accessionList[39:145]))
#for a,accession in enumerate(accessionList[39:145]):

pMetricList = np.zeros(len(accessionList))
for a,accession in enumerate(accessionList):
    print(accession)
    annotation = allMeta[accession]
    annotationFolder = annotation['folder']

    # get the track file
    prob_file = annotationFolder + 'probs_v04.csv'
    df = pd.read_csv(prob_file)

    label_term_mapping = {}
    mnemonics_file = annotationFolder + 'mnemonics_v04.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term
            
    labelCount = len(label_term_mapping)

    # fetch all the probs from the df
    allProbs = np.zeros(labelCount)
    for i in range(labelCount): # for each label
        l = str(i)
        term = label_term_mapping[l] # to find the column of the df
        # >>>>
        
        if 'Enhancer' in term:
            allProbs[i] = df._get_value(i, 'Enhancer') + df._get_value(i, 'EnhancerLow')
        else:
            if 'Promoter' in term:
                allProbs[i] = df._get_value(i, 'Promoter') + df._get_value(i, 'PromoterFlanking')
            else:
                allProbs[i] = df._get_value(i, term)

    
        #allProbs[i] = df._get_value(i, term)

    sortedPs = np.sort(allProbs)

    #sib = np.mean(sortedPs[0:int(labelCount*.6)])
    #sib = np.mean(sortedPs[0:5])
    sib = np.quantile(allProbs, .5)
    pMetricList[a] = sib

sortedMP = np.sort(pMetricList)
sortedMPA = np.argsort(pMetricList)
plt.plot(sortedMP)
plt.grid()
plt.plot([16, 16], [.35, .7])
#plt.boxplot(pMetricList)
#plt.show()
figFile = plotFolder + 'median_probs_dist.pdf'
plt.savefig(figFile)

outputFile = dataFolder + 'medianProbFile_segway.txt'
with open(outputFile, 'w') as f:
    for index in sortedMPA:
        accession = accessionList[index]
        medProb = pMetricList[index]
        f.write('%s\t%.3f\n' %(accession, medProb))
    

########################################
# 2. get the samples chromHMM coverage of quiescent
########################################

book =  'EnhA1 EnhA2 EnhBiv EnhG1 EnhG2 EnhWk Het Quies ReprPC ReprPCWk TssA TssBiv TssFlnk TssFlnkD TssFlnkU Tx TxWk ZNF/Rpts'
chromLabels = book.split()
print(len(chromLabels))

    # read each line of chromFile, get Quiescent coverage
count = 0
coverageList = {}
for accession in accessionList[68:]:
    annotation = allMeta[accession]
    print(annotation['chromFile'])

    if accession == 'ENCSR273LUT': # the file is empty, skipping it for now
        continue

    #if (('38batch' in annotation['folder']) or ('May11' in annotation['folder'])) and not(annotation['RNAseqFile'] == 'none') and not(annotation['chromFile'] == 'none'):
    if  not(annotation['chromFile'] == 'none'):
        print(accession)
        count+=1
        print(count)
        coverage = np.zeros(len(chromLabels))

        annotationFolder = annotation['folder']
        print(annotationFolder)

        chmmFileName = annotation['chromFile'].split('/')[-1]

        if (('38batch' in annotation['folder']) or ('May11' in annotation['folder'])):
            annFile = annotationFolder + 'sorted_' + chmmFileName
        else:
            annFile = annotationFolder + chmmFileName

            
        # prepare ann
        if annFile.endswith('.gz'):
            os.system('gunzip %s' %(annFile))
            annFile = annFile[0:-3]
        else:
            os.system('gunzip %s.gz' %(annFile))

        with open(annFile, 'r') as anns:
            for line in anns:
                #print(line)
                fields = line.strip().split()
                if fields[0] == 'chrX':
                    print('break')
                    break
                ann_start = int(fields[1])
                ann_end = int(fields[2])
                ann_cluster = fields[3]
                labelInd = chromLabels.index(ann_cluster)

                coverage[labelInd] += ann_end - ann_start
            
        coverageList[accession] = coverage
        os.system('gzip %s' %(annFile))

coverageList.pop('ENCSR273LUT')

accessionsWithChrom = list(coverageList.keys())
chromQs = {}
segPVals = {}
for accessionInd in sortedMPA:
    accession = accessionList[accessionInd]
    if accession in accessionsWithChrom:
        chromQs[accession] = np.sum((coverageList[accession][7])) / np.sum(coverageList[accession])
        segPVals[accession] = pMetricList[accessionInd]

chromQsArr = np.zeros(len(chromQs))
segPsArr = np.zeros(len(chromQs))
for i,q in enumerate(chromQs):
    chromQsArr[i] = chromQs[q]
    segPsArr[i] = segPVals[q]

plt.scatter(chromQsArr, segPsArr)
plt.plot(chromQsArr)
plt.scatter(range(len(chromQsArr)), chromQsArr)
plt.grid()
plt.ylim([.2,1])
plt.show()

print(np.corrcoef(chromQsArr, segPsArr))

# save the plot, report the observation. 
        
