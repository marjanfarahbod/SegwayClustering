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

    mydata = np.loadtxt(lengthFile, dtype = int, skiprows=1)
    mydata = pd.read_csv(lengthFile, sep='\t')

    print('length file read')

    totalbpCount = mydata['length'].sum(axis=0)

    mybins = np.linspace(1.9,6.1, 42)
    fix, axs = plt.subplots(16, 1, figsize=(4,10)) 
    for i in range(16):


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
    axs[15].set_xticklabels(labels, rotation=90)
    plotFile = plotFolder + accession + '/length_dist_hist.pdf'
    plt.savefig(plotFile)
    #plt.tight_layout()
    plt.close('all')
    
    plt.show()



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


