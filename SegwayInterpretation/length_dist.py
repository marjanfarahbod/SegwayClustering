# length file:

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
from util import *

# mapping Segway accession to a number - since the accession is too big or long, this I can do with that accession tissue thing.

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
dataSubFolder = 'testBatch105/fromAPI/'
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/testBatch105/'

inputFile = dataFolder + dataSubFolder + 'metaInfo.pkl'
with open(inputFile, 'rb') as f:
    annMeta = pickle.load(f)

inputFile = dataFolder + dataSubFolder + 'biosample_tissue_info.pkl'
with open(inputFile, 'rb') as f:
    tissue_info = pickle.load( f)

    
inputFile = dataFolder + 'testBatch105/' + 'runID_accession_map_105run.pkl'
with open(inputFile, 'rb') as f:
    runID_map = pickle.load(f)

sampleFolder_list = list(tissue_info.keys())

accession_IDmap = {y: x for x, y in runID_map.items()}

segtoolFolder = dataFolder + 'testBatch105/all_segtools/'

runID_list = list(runID_map.keys())

for runID in runID_list:
    
    accession = runID_map[runID]
    
    lengthFile = segtoolFolder + runID + '/length_distribution/length_distribution.tab'

    #mydata = np.loadtxt(lengthFile, dtype = int, skiprows=1)
    #print('length file read')
    mydata = pd.read_csv(lengthFile, sep='\t')
    print('length file read')

    totalbpCount = mydata['length'].sum(axis=0)

    segmentSizesFile = segtoolFolder + runID + '/length_distribution/segment_sizes.tab'
    segment_stats = pd.read_csv(segmentSizesFile, sep='\t')
    
    labelCount = segment_stats.shape[0]-1
    mybins = np.linspace(1.9,6.1, 42)
    fig, axs = plt.subplots(labelCount, 1, figsize=(4,10))
    for i in range(labelCount):

        print(i)

        sib = mydata.loc[mydata['label'] ==i]
        
        if len(sib) == 0:
            break
        print(i)

        d = sib['length'] #[1:1000]
        labelGenomeCoverage = sum(d) / totalbpCount
        labelMedian = np.median(d)
        plotData = np.log10(d.to_numpy())
    
        count = np.count_nonzero(plotData > 4)
        print(round(count/len(sib), 4))

        #plotData = numpy.copy(dnp)
        #plotData = np.log10(plotData)

        n, bins, patches = axs[i].hist(x=plotData, bins=mybins)
        #plt.xticks(ticks=[2, 2.27, 2.44, 2.56, 2.66, 2.76, 2.87, 2.97, 4, 5, 6], labels=['100', '200', '300', '400', '500','600','700-800','900-1000', '1e4', '1e5', '1e6'], rotation=60)
        #axs[0].text(3,3, 'here')
        #plt.close('all')

        yticks = [0, max(n)]
        ylabels = [0, round(max(n)/len(sib), 2)]
        axs[i].set_yticks(yticks)
        axs[i].set_yticklabels(ylabels, fontsize=8)
        axs[i].set_ylim([0, max(n)*1.2])
        axs[i].set_ylabel(i, rotation=90)
        
        ticks = [2, 2.27, 2.44, 2.56, 2.66, 2.76, 2.87, 2.97, 4, 5]
        axs[i].set_xticks(ticks)
        axs[i].set_xticklabels([])
        axs[i].set_xlim([1.9, 4.1])

        axs[i].text(3.5, max(n)*.7, 'med=%d' %(labelMedian), fontsize=8)
        axs[i].text(3.5, max(n)*.3, 'cov=%.3f' %(labelGenomeCoverage), fontsize=8)

    ticks = [2, 2.27, 2.44, 2.56, 2.66, 2.76, 2.87, 2.97, 4, 5]
    labels=['100', '200', '300', '400', '500','600','700-800','900-1000', '1e4', '1e5']
    #axs[15].set_xticks(ticks)
    axs[labelCount-1].set_xticklabels(labels, rotation=90)
    plotFile = plotFolder + accession + '/length_dist_hist.pdf'
    plt.savefig(plotFile)
    #plt.tight_layout()
    plt.close('all')
    
    plt.show()


# get that median coverage plot
# General data folder
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
dataSubFolder = 'testBatch105/fromAPI/'

inputFileName = 'all_annInfo_list.pkl'
inputFile = dataFolder + dataSubFolder + inputFileName
with open(inputFile, 'rb') as f:
    annInfo_list = pickle.load(f)

# getting the data for plot
segwayLabels = ['Enhancer_low', 'Enhancer', 'Promoter_flanking', 'Promoter', 'Transcribed', 'CTCF', 'K9K36', 'Bivalent', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent', 'Unclassified']
coverage_data = {}
for term in segwayLabels:
    coverage_data[term] = []

group_coverage_data = {}
for term in segwayLabels:
    group_coverage_data[term] = np.zeros(len(annInfo_list))

median_data = {}
for term in segwayLabels:
    median_data[term] = []

for ann in annInfo_list:

    index = ann['index']

    clusters = ann['segway_anns']['clusters']
    totalbp = 0
    bp_counts = {}
#    segment_counts = {}
    for cluster in clusters:
        totalbp += clusters[cluster].bp_count
        cluster_label = cluster.split('_')[0]
        bp_counts[cluster_label] = clusters[cluster].bp_count
#        segment_counts[cluster_label] = clusters[cluster].region_count
    
    # get the coverage: a dictionary from labels to coverage
    bp_coverage = {}
    for label in bp_counts:
        bp_coverage[label] = np.round(bp_counts[label]/totalbp, 4)

    # get the label from label from mnemonics: a dictionary from labels to terms
    sampleFolderAdd = dataFolder + dataSubFolder + ann['accession'] + '/'
    label_term_mapping = {}
    mnemonics_file = sampleFolderAdd + 'mnemonics_v02.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term


    # get the median from the length dist: ann array indexed by label
    accession = ann['accession']
    runID = accession_IDmap[accession]
    segment_file = segtoolFolder + runID + '/length_distribution/segment_sizes.tab'
    mydata = np.loadtxt(segment_file, dtype = float, skiprows=2)
    medians = mydata[:,3]


    # for each label, have a list of items and values
    for label in label_term_mapping.keys():
        clusterID = str(index) + '___' + label
        term = label_term_mapping[label]
        cov = {clusterID:bp_coverage[label]}
        coverage_data[term].append(cov)
        med = (clusterID, medians[int(label)])
        median_data[term].append(med)

        term_index = segwayLabels.index(term)
        group_coverage_data[term][index] += cov[clusterID]


    # plot the individual label coverage
    x = range(len(bp_coverage))
    y = np.zeros(len(x))
    xticklabels = []
    xtick_ind = 0
    total_term_coverage = np.zeros(len(segwayLabels))
    total_xticks = []
    for i, term in enumerate(segwayLabels):
        print(term)
        for label in label_term_mapping.keys():
            if label_term_mapping[label] == term:
                total_term_coverage[i] += bp_coverage[label]
                y[xtick_ind] = bp_coverage[label]
                ticklabel = label + '_' + label_term_mapping[label]
                xticklabels.append(ticklabel)
                xtick_ind += 1

    # plot the total label coverage
    fig, axs = plt.subplots(1, 2, figsize=(8,6))
    #plt.scatter(x, y)
    axs[0].bar(x,y)
    axs[0].set_xticks(range(xtick_ind))
    axs[0].set_xticklabels(xticklabels, rotation=90)
    axs[0].set_ylim([0,max(total_term_coverage)+.05])
    axs[0].set_ylabel('genome coverage')

    axs[1].bar(range(len(segwayLabels)), total_term_coverage)
    axs[1].set_xticks(range(len(segwayLabels)))
    axs[1].set_xticklabels(segwayLabels, rotation=90)
    axs[1].set_ylim([0,max(total_term_coverage)+.05])
    plt.tight_layout()

    figFile = plotFolder + accession + '/genome_coverage.pdf'
    plt.savefig(figFile)

    plt.show()
    


# do the two plots
sampleCount = len(annInfo_list)

# group coverage plot
for i,term in enumerate(segwayLabels[0:11]):
    y = group_coverage_data[term]
    y = y[y>0]
    x = np.random.uniform(-.25, .25, len(y)) + i
    plt.scatter(x,y , c = 'black', alpha = .3)
    

plt.xticks(range(len(segwayLabels)-1), segwayLabels[0:11], rotation = 90)
plt.tight_layout()
plotFile = plotFolder + 'group_coverage_v02.pdf'
plt.savefig(plotFile)
plt.show()

# median plot
for i,term in enumerate(segwayLabels[0:11]):
    y = np.zeros(len(median_data[term]))
    for j, item in enumerate(median_data[term]):
        y[j] = median_data[term][j][1]
    x = np.random.uniform(-.25, .25, len(y)) + i
    plt.scatter(x,y , c = 'black', alpha = .05)
    
plt.xticks(range(len(segwayLabels)-1), segwayLabels[0:11], rotation = 90)
plt.tight_layout()
plt.show()

# checking the large promoter thing:
counter = 0
prom_data = median_data['Promoter']
for item in prom_data:
    if item[1]>= 800:
        print(item[0])
        counter+=1


# coverage plot 

########################################
# DRAFT    
########################################

fix, axs = plt.subplots(16, 1, figsize=(6,12))
fix, axs = plt.subplots(2, 1, figsize=(6,12)) 
n, bins, patches = plt.hist(x=plotData, bins=book)#, density=True)
#axs[4] = plt.hist(x=plotData, bins=book, density=True)
axs[1].hist(x=plotData, bins=book, density=True)
ticks = [2, 2.27, 2.44, 2.56, 2.66, 2.76, 2.87, 2.97, 4, 5, 6]
labels=['100', '200', '300', '400', '500','600','700-800','900-1000', '1e4', '1e5', '1e6']
axs[1].set_xticks(ticks)
axs[1].set_xticklabels(labels, rotation=60)
plt.xticks(ticks=[2, 2.27, 2.44, 2.56, 2.66, 2.76, 2.87, 2.97, 4, 5, 6], labels=['100', '200', '300', '400', '500','600','700-800','900-1000', '1e4', '1e5', '1e6'], rotation=60)
#axs[0].text(3,3, 'here')
#plt.close('all')
plt.tight_layout()

plt.show()
plt.close('all')

count = np.count_nonzero(plotData > 6)
plotData[plotData>6] = 6
sns.set_style('darkgrid')
bins = np.linspace(1.9,6, 50)
#book = np.log10(np.linspace(90,1090, 11))
bookplus = np.linspace(3.09, 6.1, 30)
bookall = np.concatenate((book, bookplus), axis=0)
book = np.linspace(1.9,6.1, 42)
myhist = sns.histplot(plotData, stat='percent', bins=bins)
plt.xticks(range(2,7), [100, 1000, 10000, 100000, 1000000])
plt.show()
 
n, bins, patches = plt.hist(x=plotData, bins=book, density=True)
plt.xticks(ticks=[2, 2.27, 2.44, 2.56, 2.66, 2.76, 2.87, 2.97, 4, 5, 6], labels=['100', '200', '300', '400', '500','600','700-800','900-1000', '1e4', '1e5', '1e6'], rotation=60)
ax.text(3,3, 'here')
#plt.close('all')
plt.tight_layout()

plt.show()
plt.grid(axis='y', alpha=0.75)
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('My Very Own Histogram')
plt.text(23, 45, r'$\mu=15, b=3$')
maxfreq = n.max()
# Set a clean upper y-axis limit.
plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)


Q1 = np.percentile(dnp, 25)
Q3 = np.percentile(dnp, 75)
IQ = Q3 - Q1

uof = Q3 + 3*IQ

#print(dnp[book])
#count = np.count_nonzero(dnp > uof)


