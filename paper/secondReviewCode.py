# 1. Getting the correlation of AUCs and posterior probabilities
# 2. Getting the new labels for the old classifier data
# 3. Figure 2 Max requested predictivity thing
# 4. Figure 2 another one 

########################################
# 0. Initials
########################################


import util # this is for features_from_segtools_dir
import gzip
import pickle
import pandas as pd
import numpy as mp

import seaborn as sns
import matplotlib.pyplot as plt

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'

# review plot folder
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/review02Plots/'
########################################
# 1. getting the correlation of AUCs and posterior probabilities.
########################################

# 1. find the AUCs for the samples.


file = dataFolder + 'Segway_logisticRegression_output_v04_expFile_v02.pkl'
with open(file, 'rb') as f:
    regOutput = pickle.load(f)

# which key we are going for? both gene and promoter, unbalanced, and the 1st choice: 1

key = ('genes', 1, 'unbalanced')
key = ('promoter', 1, 'unbalanced')

book = regOutput[key]
myAUCs = [kado[1] for kado in book]  # got the AUCs
myAccessions = [kado[0] for kado in book] # got the accessions

# 2. find the posterior from the samples.
meanProbs = []
medianProbs = []
quantProbs = []
quantProbsLow = []
highmeanProbs = []
for accession in myAccessions:

    annotation = allMeta[accession]
    annotationFolder = annotation['folder']
    
    # the mnemonics file
    label_term_mapping = {}
    mnemonics_file = annotationFolder + 'mnemonics_v04.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term

    labelCount = len(label_term_mapping)

    # get the track file
    prob_file = annotationFolder + 'probs_v04.csv'
    df = pd.read_csv(prob_file)
    df = pd.read_csv(prob_file, usecols=list(range(1, df.shape[1])))

    highProbs = (df.max(axis=1))
    '''
    meanProbs.append(highProbs.mean())
    medianProbs.append(highProbs.median())
    quantProbs.append(highProbs.quantile(.75))
    quantProbsLow.append(highProbs.quantile(.25))
    '''
    myArr = np.array(highProbs)

    meanProbs.append(highProbs.mean())    
    highmeanProbs.append(np.mean(np.sort(myArr)[-3:]))
    medianProbs.append(highProbs.median())
    quantProbs.append(highProbs.quantile(.75))
    quantProbsLow.append(highProbs.quantile(.25))


# 3. do the plots

highmeanProbsArr  = np.array(highmeanProbs)
meanProbsArr  = np.array(meanProbs)
medianProbsArr = np.array(medianProbs)
myAUCsArr = np.array(myAUCs)
quantileArr = np.array(quantProbs)
quantileLowArr = np.array(quantProbsLow)
print(np.corrcoef(quantileLowArr, myAUCsArr))
print(np.corrcoef(meanProbsArr, myAUCsArr))
print(np.corrcoef(highmeanProbsArr, myAUCsArr))
print(np.corrcoef(quantileArr, myAUCsArr))
print(np.corrcoef(medianProbsArr, myAUCsArr))

sortInds = np.argsort(myAUCsArr)
sortInds = np.argsort(quantileArr)
sortInds = np.argsort(meanProbsArr)
sortInds = np.argsort(medianProbsArr)

plt.scatter(meanProbsArr[sortInds], myAUCsArr[sortInds])
#title('mean probs, prediction AUC')
plt.xlabel('meanProbs')
plt.ylabel('AUCs')

plt.scatter(medianProbsArr[sortInds], myAUCsArr[sortInds])
#title('median probs, prediction AUC')
plt.xlabel('medianProbs')
plt.ylabel('AUCs')

plt.scatter(quantileArr[sortInds], myAUCsArr[sortInds])
#title('0.75 percentile probs, prediction AUC')
plt.xlabel('0.75percentile Probs')
plt.ylabel('AUCs')

plt.title('gene predictions')
plt.title('promoter predictions')

plt.show()
plt.close('all')

plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/'
figFile = plotFolder + 'AUC_75Probs_gene.pdf'
figFile = plotFolder + 'AUC_meanProbs_gene.pdf'
figFile = plotFolder + 'AUC_medianProbs_gene.pdf'

figFile = plotFolder + 'AUC_75Probs_promoter.pdf'
figFile = plotFolder + 'AUC_meanProbs_promoter.pdf'
figFile = plotFolder + 'AUC_medianProbs_promoter.pdf'

plt.savefig(figFile)
plt.close('all')


# 2. Getting the new labels for the old classifier data
#########################################
# (see code classifier_train.py)

runFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/runs/'

# get the file: classifier_data.tab 
classifier_tab_fname = dataFolder + 'classifier_data.tab' # this is the old classifier data

# the old labels for the old data (2019, 2019)
# oldDataLabels_file = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/backup_allinterpretation/interpretation/exp/02_2020-06-03_interpretation/16_2020-08-06_run/label_mappings.txt'


# list of labels
# list of features - classifier input
feature_names = ['(09) initial exon', '(01) H3K9me3', '(10) initial intron', '(02) H3K27me3', '(11) internal exons', '(04) H3K4me3', "(16) 3' flanking (1000-10000 bp)", '(12) internal introns', '(03) H3K36me3', '(13) terminal exon', '(06) H3K4me1', '(14) terminal intron', "(07) 5' flanking (1000-10000 bp)", '(05) H3K27ac', "(15) 3' flanking (1-1000 bp)", "(08) 5' flanking (1-1000 bp)"]


feature_names_plotOrder = ['(09) initial exon',  "(08) 5' flanking (1-1000 bp)",'(10) initial intron', '(11) internal exons','(12) internal introns', '(13) terminal exon',  '(14) terminal intron', "(15) 3' flanking (1-1000 bp)" ,"(07) 5' flanking (1000-10000 bp)","(16) 3' flanking (1000-10000 bp)",'(04) H3K4me3','(05) H3K27ac','(06) H3K4me1', '(03) H3K36me3','(02) H3K27me3','(01) H3K9me3']


# 2.1 Get the old mean-feature-label from the old training data, based on the labels and mean features
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

classifier_data_frame = pd.read_csv(classifier_tab_fname, sep="\t")
old_classifier_labels = new_label_archive # this is from the file classifier_train.py, section 1.0 - these are the old labels. new_label_archive[concatenation_key][orig_label] = new_label

all_ref_bio_labels = set.union(*map(set, map(lambda x: x.values(), old_classifier_labels.values())))

# removing rows with 'Unclassified label' from the dataframe
archive_label_list = []
dropList = []
for i in range(len(classifier_data_frame)):
    if (old_classifier_labels[classifier_data_frame.iat[i,0]][classifier_data_frame.iat[i,1]]) == 'Unclassified':
        dropList.append(i)
    else:
        archive_label_list.append(old_classifier_labels[classifier_data_frame.iat[i,0]][classifier_data_frame.iat[i,1]])

classifier_data_frame_dropped = classifier_data_frame.drop(dropList)

# get the model
runFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/runs/'
model_file = runFolder + "run03/model_300_reg.020_auc0.89V04.pickle.gz"
with gzip.open(model_file, "r") as f:
    the_model = pickle.load(f)

# prepare the input data 
feature_names = ['(09) initial exon', '(01) H3K9me3', '(10) initial intron', '(02) H3K27me3', '(11) internal exons', '(04) H3K4me3', "(16) 3' flanking (1000-10000 bp)", '(12) internal introns', '(03) H3K36me3', '(13) terminal exon', '(06) H3K4me1', '(14) terminal intron', "(07) 5' flanking (1000-10000 bp)", '(05) H3K27ac', "(15) 3' flanking (1-1000 bp)", "(08) 5' flanking (1-1000 bp)"]

mydata = classifier_data_frame_dropped[feature_names]
classifier_input_features = mydata.to_numpy()

# get the output
all_labels = the_model.predict(classifier_input_features) # classifier output - new
probs = the_model.predict_proba(example_features) # probs

# now we have all_labels and archive_label_list, we need to make two heatmaps with them

# first, the old_classifier_labels plot
meanMat = np.zeros((8, 16)) 
wholeMat = classifier_input_features
from scipy.stats import zscore
plot_data_z = zscore(wholeMat, axis = 0)
# plot_data_z_thr = np.where(plot_data_z > 1, 1.1, plot_data_z)

segwayOldLabel_list = ['Promoter', 'Enhancer', 'RegPermissive', 'Bivalent', 'Transcribed', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent']

for i in range(len(segwayOldLabel_list)):
    book = np.where(np.array(archive_label_list) == np.array(segwayOldLabel_list[i]))[0]
    subMat = plot_data_z[book, :]
    meanMat[i, :] = np.mean(subMat, axis=0)

# change meanMat to dataframe
import matplotlib.colors as colors
meanMat_thr = np.where(meanMat >2, 2, meanMat)
plot_old_data = pd.DataFrame(meanMat_thr, columns=feature_names, index=segwayOldLabel_list)
cmap = plt.cm.coolwarm
norm = colors.BoundaryNorm(np.arange(-1, 2, .5), cmap.N)
book = plot_old_data[feature_names_plotOrder]
sns.heatmap(book, center=0, cmap=cmap, norm=norm, linewidth=.1, linecolor='white')



plt.title('average z-score of the feature values\n among the 210 training labels from Segway2019, \n for 2019 interpretation terms')
plt.tight_layout()
plt.show()


figFile = plotFolder + 'oldSegway_trainData_plot_oldLabels_review02.pdf'
print(figFile)
plt.savefig(figFile)
plt.close('all')

# getting the new plot
segwayStates = ['Promoter', 'PromoterFlanking', 'Enhancer', 'EnhancerLow', 'Bivalent', 'CTCF', 'Transcribed', 'K9K36', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent']

meanMat = np.zeros((11, 16)) 
for i in range(len(segwayStates)):
    book = np.where(np.array(all_labels) == np.array(segwayStates[i]))[0]
    subMat = plot_data_z[book, :]
    meanMat[i, :] = np.mean(subMat, axis=0)

# change meanMat to dataframe
import matplotlib.colors as colors
meanMat_thr = np.where(meanMat >2, 2, meanMat)
plot_old_data = pd.DataFrame(meanMat_thr, columns=feature_names, index=segwayStates)
cmap = plt.cm.coolwarm
norm = colors.BoundaryNorm(np.arange(-1, 2, .5), cmap.N)
book = plot_old_data[feature_names_plotOrder]
sns.heatmap(book, center=0, cmap=cmap, norm=norm, linewidth=.1, linecolor='white')

plt.title('average z-score of the feature values\n among the 210 training labels from Segway2019, \n for new interpretation terms')
plt.tight_layout()
#plt.show()


figFile = plotFolder + 'oldSegway_trainData_plot_newLabels_review02.pdf'
print(figFile)
plt.savefig(figFile)
plt.close('all')

########################################
# 3. Figure 2 Max requested predictivity thing
########################################

# for this section see the loop in section 1 paperFigure02.py "get the intermediate data for the matrix"
# we need to get the raw fraction of the occurrences for each bar, all the rest stay the same. 
# code below is copied from that file, so mind all the imports and variable definitions
            
rnaAccessionList = []
for accession in accessionList:

    annotation = allMeta[accession]
    #if ((('38batch' in annotation['folder']) or ('May11' in annotation['folder'])) and not(annotation['RNAseqFile'] == 'none')):
    if((annotation['RNAseqFile'] != 'none')):
        rnaAccessionList.append(accession)

#allBarValues = np.zeros((94, 3, 11, 160)) # a, accession; s, segwayStates; m, expMat; b, bars: log(obs/expected) [note:fractions]
allObsMats = np.zeros((94, 3, 11, 160)) # a, accession; s, segwayStates; m, expMat; b, bars: obs [note:fractions]
barPresence =  np.zeros((94, 11)) # some samples don't have the some labels
geneCountBars = np.zeros((94, 3, 11, 160)) # a, accession; s, segwayStates; m, expMat; b, bars: obs [just the count]
for a, accession in enumerate(rnaAccessionList):
    print(a)

    annotation = allMeta[accession]

    annotationFolder = annotation['folder']
    expFile = annotationFolder + 'defaultExp_5kg_expSummary_newSort_Q30_toTerms.pkl'
    print(expFile)
    with open(expFile, 'rb') as file:
        expMats = pickle.load(file)

    # 1.1 get the mnemonics - we need this fraction coverage
    # >>>>>>>>>>>>>>>>>>>>>>>>>>
    label_term_mapping = {}
    mnemonics_file = annotationFolder + 'mnemonics_v04.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term
            
    labelCount = len(label_term_mapping)


    # 1.2 get the fraction of coverage for each basepair
    # >>>>>>>>>>>>>>>>>>>>>>>>>>
    lfile = length_files[accession]
    fractions = np.zeros(segwayStateCount)
    with open(lfile, 'r') as inputFile:
        header = inputFile.readline()
        header = inputFile.readline()
        for line in inputFile:
            fields = line.strip().split()
            label = fields[0]
            stateInd = segwayStates.index(label_term_mapping[str(label)])
            coverage = fields[-1]
            #print(label, coverage)
            fractions[stateInd] += float(coverage)

    barPresence[a, fractions>0] = 1  # filling the presence
    fractions[fractions == 0] = .0001 # for the calculations only

    # fig, axs = plt.subplots(segwayStateCount, 3, figsize=(12,8))
                
    # xticks = [30, 130]
    # xticksLabels = ['TSS', 'TTS']

    indexList = np.array(list(range(160)))
    for m in range(3):

        thisMat = expMats[m]

        geneCountBars[a, m, :, :] = np.copy(thisMat)

        # make it the ratio
        thisMat = thisMat / np.sum(thisMat, axis=0)[np.newaxis,:]
        # just saving the observed here (thisMat so far is the observed)
        observedMat = np.copy(thisMat)
        allObsMats[a, m, :, :] = observedMat

        # we don't need the following part, and just the observed mats
        '''
        # versus expected
        thisMat = thisMat / fractions[:, None]
        logMat = np.log10(thisMat, out=np.zeros_like(thisMat), where=(thisMat!=0))
        allBarValues[a, m, :, :] = logMat
        '''
# we just got the observed mats, now we need to plot the mean from the observed mat for the three categories. 

meanBars = np.zeros((3, 11, 160)) # m, expMat; s, segwayState; b, bars

for m in range(3): # for each matrix
    for s in range(11): # for each label
        thisBarPresence = barPresence[:, s]
        for b in range(160): # for each bar
            barValues = allObsMats[:, m, s, b]
            values = barValues[thisBarPresence ==1]
            meanBars[m,s,b] = np.mean(values)

#########################################

fig, axs = plt.subplots(segwayStateCount, 3, figsize=(12,8))

xticks = [30, 130]
xticksLabels = ['TSS', 'TTS']

indexList = np.array(list(range(160)))

# make it the ratio
#thisMat = thisMat / np.sum(thisMat, axis=0)[np.newaxis,:]
# versus expected
#thisMat = thisMat / fractions[:, None]
#logMat = np.log10(thisMat, out=np.zeros_like(thisMat), where=(thisMat!=0))
for m in range(3):
    thisMeanBars = meanBars[m, :, :]
    for s in range(segwayStateCount):
        positiveInds = indexList[thisMeanBars[s,:] >= 0]
        negativeInds = indexList[thisMeanBars[s,:] < 0]
        posBar = np.copy(thisMeanBars[s, :])
        posBar[negativeInds]=0
        axs[s,m].bar(range(160), posBar, color=[249/256, 80/256, 97/256,1], width=1)
        negBar = np.copy(thisMeanBars[s, :])
        negBar[positiveInds]=0
        axs[s,m].bar(range(160), negBar, color=[50/256,164/256,249/256,1], width=1)
        axs[s,m].plot([0, 160], [.5, .5], linestyle='--', linewidth=.5)
        axs[s,m].set_ylim((-.1,1.1))
        ylabel = segwayStates[s]
        #axs[i,0].text(60, .55, ylabel, fontsize=8)
        axs[s,m].set_xticks(xticks)
        axs[s,m].set_yticks([0,.5, 1])

        axs[s,m].set_yticklabels([0,.5, 1], fontsize=8)
        if m == 0:
            axs[s,m].set_ylabel(ylabel, rotation=0, fontsize=8, labelpad=3, ha='right', va='center')

    axs[segwayStateCount-1,m].set_xticklabels(xticksLabels)
    axs[segwayStateCount-1,m].set_xlabel('Position relative to gene')


axs[0,0].set_title('For genes with zero expression levels, \nratio of the genes with each label\n at different genomic positions', fontsize = 9)
axs[0,1].set_title('For genes with bottom 30% expression levels, \nratio of the genes with each label\n at different genomic positions', fontsize = 9)
axs[0,2].set_title('For genes with top 70% expression levels, \nratio of the genes with each label\n at different genomic positions', fontsize = 9)

plt.show()


#figFile = plotFolder + figureFolder + 'transcriptomicPlot_median_QartileLines.pdf'
figFile = plotFolder + 'fig2_a_predictive_ratio.pdf'
plt.savefig(figFile)
plt.close('all')

########################################
# 4. Figure 2 another one
########################################

# the raw gene counts in each bar are in the geneCountBars, this is the number that Max wants for each bar b at the 3rd column, row i: (M1/M)/(N1/N)
# M1: value at that bar (number of genes with label i and top 70% expression) - topLabelBar
# M: sum of values at bar b, label i for column 1, 2, and 3 (number of genes with label i) - labelBarTotal
# N: sum of column 2, for all 3 matrices - totalGeneCount
# N1: sum of column 2, for matrix 2, top 70 geneCount - topGeneCount
# we are caclculating (M1/M)/(N1/N)
# (topLabelBar/labelBarTotal)/(topGeneCount/totalGeneCount)

# the result is a one matrix of [11, 160] for each of the samples.
# expMats [3, 11, 160]
result = np.zeros((94, 11, 160))
# I already have the barPresence so I skipped it for this part
# I also have the geneCountBars [94, 3, 11, 160], which has everything I need

for i in range(94): # sample count, for each sample
    zeroGeneCount = np.sum(geneCountBars[i, 0, :, :], 0)[1]
    lowGeneCount = np.sum(geneCountBars[i, 1, :, :], 0)[1]
    topGeneCount = np.sum(geneCountBars[i, 2, :, :], 0)[1]
    totalGeneCount = zeroGeneCount + lowGeneCount + topGeneCount
    for l in range(11): # for each label
        for b in range(160): # for each bar    
            labelBarTotal = np.max(((geneCountBars[i, 0, l, b] + geneCountBars[i, 1, l, b] + geneCountBars[i, 2, l, b]), 1))
            topLabelBar = geneCountBars[i, 2, l, b]
            result[i, l, b] = (topLabelBar/labelBarTotal)/(topGeneCount/totalGeneCount)


meanBars = np.zeros((11, 160)) # m, expMat; s, segwayState; b, bars

for s in range(11): # for each label
    thisBarPresence = barPresence[:, s]
    for b in range(160): # for each bar
        barValues = result[:, s, b]
        values = barValues[thisBarPresence ==1]
        meanBars[s,b] = np.mean(values)

# plotting the results
fig, axs = plt.subplots(segwayStateCount, 1, figsize=(5, 10))

xticks = [30, 130]
xticksLabels = ['TSS', 'TTS']

indexList = np.array(list(range(159)))

thisMeanBars = meanBars[:, 1:]
thisMeanBarsLog = np.log10(thisMeanBars, out=np.zeros_like(thisMeanBars), where=(thisMeanBars!=0))
# thisMeanBars = result[0, :, :]
for s in range(segwayStateCount):
    positiveInds = indexList[thisMeanBarsLog[s,:] >= 0]
    negativeInds = indexList[thisMeanBarsLog[s,:] < 0]
    posBar = np.copy(thisMeanBarsLog[s, :])
    posBar[negativeInds]=0
    axs[s].bar(range(159), posBar, color=[249/256, 80/256, 97/256,1], width=1)
    negBar = np.copy(thisMeanBarsLog[s, :])
    negBar[positiveInds]=0
    axs[s].bar(range(159), negBar, color=[50/256,164/256,249/256,1], width=1)
    #axs[s].plot([0, 160], [.5, .5], linestyle='--', linewidth=.5)
    axs[s].set_ylim((-.8,.6))
    ylabel = segwayStates[s]
    #axs[i,0].text(60, .55, ylabel, fontsize=8)
    axs[s].set_xticks(xticks)
    axs[s].set_yticks([-.5, 0, .5])
    axs[s].set_ylabel(ylabel, rotation=0, fontsize=9, labelpad=3, ha='right', va='center')
    if s < 10:
        axs[s].tick_params(axis='x', labelbottom=False)


    #axs[s].set_yticklabels([-1, 0, .5], fontsize=8)
    #axs[s].set_ylabel(ylabel, rotation=0, fontsize=8, labelpad=3, ha='right', va='center')

axs[segwayStateCount-1].set_xticklabels(xticksLabels)
axs[segwayStateCount-1].set_xlabel('Position relative to gene')
axs[0].set_title('title goes here', fontsize = 9)

#
plt.tight_layout()
plt.show()

#figFile = plotFolder + figureFolder + 'transcriptomicPlot_median_QartileLines.pdf'
figFile = plotFolder + 'M1_M_N1_N_log.pdf'
figFile = plotFolder + 'M1_M.pdf'
plt.savefig(figFile)
plt.close('all')



########################################
# please see paperFigure02
########################################

fig, axs = plt.subplots(segwayStateCount, 1, figsize=(4,8))

xticks = [30, 130]
xticksLabels = ['TSS', 'TTS']

indexList = np.array(list(range(160)))

thisMeanBars = np.log(meanFRQ)
for s in range(segwayStateCount):
    positiveInds = indexList[thisMeanBars[s,:] >= 0]
    negativeInds = indexList[thisMeanBars[s,:] < 0]
    posBar = np.copy(thisMeanBars[s, :])
    posBar[negativeInds]=0
    axs[s].bar(range(160), posBar, color=[249/256, 80/256, 97/256,1], width=1)
    negBar = np.copy(thisMeanBars[s, :])
    negBar[positiveInds]=0
    axs[s].bar(range(160), negBar, color=[50/256,164/256,249/256,1], width=1)
    axs[s].set_ylim((-3,5))
    ylabel = segwayStates[s]
    #axs[i,0].text(60, .55, ylabel, fontsize=8)
    axs[s].set_xticks(xticks)
    axs[s].set_yticks([-2, 4])
    axs[s].set_ylabel(ylabel, rotation=0, fontsize=9, labelpad=3, ha='right', va='center')
    if s < 10:
        axs[s].tick_params(axis='x', labelbottom=False)

#plt.show()
axs[segwayStateCount-1].set_xticklabels(xticksLabels)
axs[segwayStateCount-1].set_xlabel('Position relative to gene')
#plt.figure(constrained_layout=True)
plt.tight_layout()
#plt.show()


#figFile = plotFolder + figureFolder + 'transcriptomicPlot_median_QartileLines.pdf'
figFile = plotFolder + 'thatFreqPlot.pdf'
plt.savefig(figFile)
plt.close('all')



########################################
# DRAFT
########################################
'''
# read the file into a dataframe
 dataFeatures = pd.read_csv(oldDataFeatures_file, sep='\t')

# reorder featureNames
reorderFeatures = list(['concatenation_key', 'orig_label'])
reorderFeatures.extend(feature_names_plotOrder)

dataFeatures_reorder = dataFeatures[reorderFeatures]

for i in range(len(dataFeatures)):
    print(dataFeatures['concatenation_key'][i], dataFeatures['orig_label'][i])

# read the labels 
oldDataLabels = pd.read_csv(oldDataLabels_file, sep='\t')

for i in range(len(oldDataLabels)):
    print(oldDataLabels['concatenation_key'][i], oldDataLabels['orig_label'][i])

# note: we have 223 labels with features, and 294 labels with "labels" - all features have labels, but features are missing for some of the labels. It's fine.

'''

# Add the new_label (which is the old classifier label) to the feature set.

# get the mean value of each feature 


# 2.2 Get the new mean-feature-label from the the old training data, using the new classifier 


'''
# reordering the features 
classifier_input_features = features_frame_to_matrix(classifier_data_frame, feature_names)
'''


# getting the right order of features for the plot

feature_names_plotOrder = ['(09) initial exon' , "(08) 5' flanking (1-1000 bp)" , '(10) initial intron', '(11) internal exons', '(12) internal introns', '(13) terminal exon', '(14) terminal intron', "(15) 3' flanking (1-1000 bp)", "(07) 5' flanking (1000-10000 bp)",  "(16) 3' flanking (1000-10000 bp)",'(04) H3K4me3', '(05) H3K27ac', '(06) H3K4me1','(03) H3K36me3', '(02) H3K27me3', '(01) H3K9me3']

feature_names_plotOrder_renamed = ['initial exon',  "5' flanking (1-1000 bp)",'initial intron', 'internal exons','internal introns', 'terminal exon',  'terminal intron', "3' flanking (1-1000 bp)" , "5' flanking (1000-10000 bp)","3' flanking (1000-10000 bp)",'H3K4me3','H3K27ac','H3K4me1', 'H3K36me3','H3K27me3','H3K9me3']

# removing rows with unclassified from the dataframe






plot_features = features_frame_to_matrix(classifier_data_frame, feature_names_plotOrder)


# TODO: find the heatmap plot and do the heatmap for old and new labels for the 210 and be done




