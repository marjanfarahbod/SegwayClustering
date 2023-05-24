# how do they compare regarding the transcription data.
# How well do they do regarding the transcription data? active prmoter
# How well do they do regarding the transcription data? transcribed region etc.

# expressed genes act very much the same between the low and high expression.

# for each gene, get the ratio of its promoter region and its region covered with each of the labels
# examine the distribution of these labels between the expressed and non expressed genes.

# do this for both chromhmm and segway. hahaha.

# 0. Initials
# 0.0 The function for extracting data - Segway
# 0.1 The function for extracting data - chromhmm
# 1. Segway data extraction
# 2. CHromHMM data extraction
# 3. The AUC curve for genes which were included
# 4. get the AUC plot

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
from QC_transcriptionComparison_util import Gene, Exon, Annotation, AnnotationClass

# General data folder
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
dataSubFolder = 'testBatch105/fromAPI/'

# GTF data structure
fileName = dataFolder + '/geneLists.pkl'
with open(fileName, 'rb') as pickleFile:
    geneListsAndIDs = pickle.load(pickleFile)
    
geneList = geneListsAndIDs[0]
geneIDList = geneListsAndIDs[1]
del geneListsAndIDs

inputFile = dataFolder + dataSubFolder + 'metaInfo.pkl'
with open(inputFile, 'rb') as f:
    annMeta = pickle.load(f)

sample_count = len(annMeta)

annAccessionList = list(annMeta.keys())
annAccession = annAccessionList[104]
print(annAccession)


# obsolete segwayLabels = ['Quiescent', 'ConstitutiveHet', 'FacultativeHet', 'Transcribed', 'Promoter', 'Enhancer', 'RegPermissive', 'Bivalent', 'LowConfidence']

# How well do they do regarding the transcription data? active prmoter
# How well do they do regarding the transcription data? transcribed region etc.

# expressed genes act very much the same between the low and high expression.

# for each gene, get the ratio of its promoter region and its region covered with each of the labels
# examine the distribution of these labels between the expressed and non expressed genes.

# do this for both chromhmm and segway. hahaha.
# 0. The function

# 1. Segway data extraction
# 2. CHromHMM data extraction

########################################
# 0.0 The function - Segway
########################################
from transcription_overlap import EPLabelCover_Segway

'''
Seway bed files are already sorted and filtered. Please see the QC_transcriptionComparison_03.py section 5
'''

# I need to do it for the 21 remaining samples with the transcriptomic data
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'

# load the meta file 
inputFileName = 'all235Annot_meta_corrected.pkl'
#inputFileName = 'all235Annot_meta.pkl'
inputFile = dataFolder +  inputFileName
with open(inputFile, 'rb') as f:
    allMeta = pickle.load(f)

accessionList = list(allMeta.keys())

count = 0
for accession in accessionList:
    annotation = allMeta[accession]
    count+=1

    if ((('38batch' in annotation['folder']) or ('May11' in annotation['folder'])) and not(annotation['RNAseqFile'] == 'none')):
        RNAFile = annotation['RNAseqFile'][0]
        print(count)
        print(accession)
        print(RNAFile)

        annotationFolder = annotation['folder']
        print(annotationFolder)

        expAccession = RNAFile[-15:-4]
        expFile = annotationFolder +  '/geneExp_dict_' + expAccession + '.pkl'
        annFile = annotation['bedFile']
        mnemFile = annotationFolder + 'mnemonics_v04.txt'
        #geneList
        #geneIDList =
        promLength = 3000

        EPLabelCover_Segway(annotationFolder, expFile, annFile, mnemFile, geneList, geneIDList, promLength)

########################################
# 0.1 The function - chromhmm
########################################

# sort the chromHMM file, and then the annotation file will be the sorted_ChromHMM_[chromFile]
count = 0
for accession in accessionList:
    annotation = allMeta[accession]

    if ((('38batch' in annotation['folder']) or ('May11' in annotation['folder'])) and not(annotation['chromFile'] == 'none')):
        count+=1
        print(count)
        print(accession)
        annotationFolder = annotation['folder']
        print(annotationFolder)
        chromFile = annotation['chromFile']
        print(chromFile)

        fileName = chromFile.split('/')[-1][0:-3]
        print(fileName)

        os.system('gunzip %s' %(chromFile))
        chromFile = chromFile[0:-3]
        print(chromFile)

        sortedBedFile = annotationFolder + 'sorted_%s' %(fileName)
        print(sortedBedFile)
        os.system('sort -V -k1,1 -k2,2 %s > %s' %(chromFile, sortedBedFile))

        os.system('gzip %s' %(chromFile))
        os.system('gzip %s' %(sortedBedFile))

'''
Active TSS
Flanking TSS
Flanking TSS upstream
Flanking TSS downstream
Strong transcription
Weak transcription
Genic enhancer 1
Genic enhancer 2
Active enhancer 1
Active enhancer 2
Weak enhancer
ZNF genes & repeats
Heterochromatin
Bivalent/Poised TSS
Bivalent enhancer
Repressed PolyComb
Weak Repressed PolyComb
Quiescent/Low
'''

book =  'EnhA1 EnhA2 EnhBiv EnhG1 EnhG2 EnhWk Het Quies ReprPC ReprPCWk TssA TssBiv TssFlnk TssFlnkD TssFlnkU Tx TxWk ZNF/Rpts'
chromLabels = book.split()

count = 0
for accession in accessionList[6:]:
    annotation = allMeta[accession]


    if (('38batch' in annotation['folder']) or ('May11' in annotation['folder'])) and not(annotation['RNAseqFile'] == 'none') and not(annotation['chromFile'] == 'none'):
        count+=1
        print(count)
        RNAFile = annotation['RNAseqFile'][0]
        print(RNAFile)
        print(accession)

        annotationFolder = annotation['folder']
        print(annotationFolder)

        expAccession = RNAFile[-15:-4]
        expFile = annotationFolder +  'geneExp_dict_' + expAccession + '.pkl'
        chmmFileName = annotation['chromFile'].split('/')[-1]
        annFile = annotationFolder + 'sorted_' + chmmFileName
        #geneList
        #geneIDList =
        promLength = 3000

        print(expAccession)
        print(expFile)
        print(annFile)

        EPLabelCover_chromHMM(annotationFolder, expFile, annFile, chromLabels, geneList, geneIDList, promLength)


        
# 3. The AUC curve for genes which were included
########################################
# get the AUC for Segway runs

inputFile = dataFolder + dataSubFolder + 'biosample_tissue_info.pkl'
with open(inputFile, 'rb') as f:
    tissue_info = pickle.load( f)

plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/testBatch105/'


aucs = {} # gene, promoter
for annAccession in annAccessionList:
    
    if (annMeta[annAccession]['expressionFile'] != 'none'):

        expFileName = annMeta[annAccession]['expressionFile']
        expAccession = re.split('_|\.', expFileName)[2]
        inputFile = dataFolder + dataSubFolder + annAccession + '/geneExp_dict_' + expAccession + '.pkl'
        print(inputFile)
        with open(inputFile, 'rb') as pickledFile:
            expression = pickle.load(pickledFile)

        label_term_mapping = {}
        mnemonics_file = dataFolder + dataSubFolder + annAccession + '/' + 'mnemonics_v04.txt'
        with open(mnemonics_file, 'r') as mnemonics:
            for line in mnemonics:
                print(line)
                label = line.strip().split()[0]
                term = line.strip().split()[1]
                label_term_mapping[label] = term

        clusterCount = len(label_term_mapping)
        
        inputFile = dataFolder + dataSubFolder + annAccession + '/' + 'exp_promoter_labelCover.pkl'
        with open(inputFile, 'rb') as f:
            transPromoMat = pickle.load(f)  # promoter, genes

        labelExpMat = transPromoMat['genes']
        sib = labelExpMat.sum(axis = 1)
        filterGene = sib > 0

        expArray = np.asarray([expression[x] for x in geneIDList])
        filterExp = expArray[filterGene,]
        filterExpMat = labelExpMat[filterGene, :]
        sumExp = filterExpMat.sum(axis = 1)
        filterExpMat = filterExpMat/sumExp[:,None]
    
        notExp = filterExp == 0

        #fig = plt.figure(figsize =(17, 7))
        #ax = fig.add_subplot(111)

        fig, axs = plt.subplots(nrows = 4, ncols=4, figsize =(10,9))

        geneAUC = []

        for i in range(clusterCount):
            print(i)
            filteredExpArray = filterExpMat[:, i]
            sib = np.argsort(filteredExpArray)[::-1]

            AUCCurve = np.zeros((len(sib), 2))
            xwalk = 0
            ywalk = 0
            area = 0 
            for j in range(len(sib)):
                if notExp[sib[j]] == True:
                    xwalk +=1
                    area+= ywalk
                else:
                    ywalk +=1
                    
                AUCCurve[j, 0] = xwalk
                AUCCurve[j, 1] = ywalk

            auc = area / (xwalk*ywalk)
            geneAUC.append(auc)
            
            pi = int(i /4)
            pj = i%4
            axs[pi, pj].plot(AUCCurve[:, 0], AUCCurve[:,1])
            if pj > 0:
                axs[pi, pj].set_yticks([])

            if pi < 3:
                axs[pi, pj].set_xticks([])

            clusterLabel = '%d' %(i)
            axs[pi, pj].text(xwalk*.48,ywalk/17, 'AUC: %.3f' %(auc))
            axs[pi, pj].text(xwalk*.012,ywalk*.8, '%s\n%d' %(label_term_mapping[clusterLabel], i))
        #plt.scatter(AUCCurve[:, 0], AUCCurve[:,1], s=4, ax=axs[1])
        #plt.boxplot(data)
        
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        tissue = tissue_info[annAccession]
        fig.suptitle('prediction for gene expression \n %s - %s' %(annAccession, tissue))
        fig.text(0.5, 0.04, 'zero expression', ha='center')
        fig.text(0.04, 0.5, 'expressed', va='center', rotation='vertical')

        plotFolder_add = plotFolder + annAccession + '/'
        figFile = plotFolder_add + 'geneAUC.pdf'
        print(figFile)
        plt.savefig(figFile, bbox_inches='tight')
        plt.close('all')

        #plt.show()


        labelExpMat = transPromoMat['promoter']
        sib = labelExpMat.sum(axis = 1)
        filterGene = sib > 0

        expArray = np.asarray([expression[x] for x in geneIDList])
        filterExp = expArray[filterGene,]
        filterExpMat = labelExpMat[filterGene, :]
        sumExp = filterExpMat.sum(axis = 1)
        filterExpMat = filterExpMat/sumExp[:,None]
    
        notExp = filterExp == 0

        #fig = plt.figure(figsize =(17, 7))
        #ax = fig.add_subplot(111)

        fig, axs = plt.subplots(nrows = 4, ncols=4, figsize =(10,9))

        promoterAUC = []

        for i in range(clusterCount):
            print(i)
            filteredExpArray = filterExpMat[:, i]
            sib = np.argsort(filteredExpArray)[::-1]

            AUCCurve = np.zeros((len(sib), 2))
            xwalk = 0
            ywalk = 0
            area = 0 
            for j in range(len(sib)):
                if notExp[sib[j]] == True:
                    xwalk +=1
                    area+= ywalk
                else:
                    ywalk +=1
                    
                AUCCurve[j, 0] = xwalk
                AUCCurve[j, 1] = ywalk

            auc = area / (xwalk*ywalk)

            promoterAUC.append(auc)
            
            pi = int(i /4)
            pj = i%4
            axs[pi, pj].plot(AUCCurve[:, 0], AUCCurve[:,1])
            if pj > 0:
                axs[pi, pj].set_yticks([])

            if pi < 3:
                axs[pi, pj].set_xticks([])

            clusterLabel = '%d' %(i)
            axs[pi, pj].text(xwalk*.48,ywalk/17, 'AUC: %.3f' %(auc))
            axs[pi, pj].text(xwalk*.012,ywalk*.8, '%s\n%d' %(label_term_mapping[clusterLabel], i))
        #plt.scatter(AUCCurve[:, 0], AUCCurve[:,1], s=4, ax=axs[1])
        #plt.boxplot(data)
        
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        tissue = tissue_info[annAccession]
        fig.suptitle('prediction for active promoter \n %s - %s' %(annAccession, tissue))
        fig.text(0.5, 0.04, 'zero exprssion', ha='center')
        fig.text(0.04, 0.5, 'expressed genes', va='center', rotation='vertical')

        plotFolder_add = plotFolder + annAccession + '/'
        figFile = plotFolder_add + 'promoterAUC.pdf'
        print(figFile)
        plt.savefig(figFile, bbox_inches='tight')
        plt.close('all')
        #plt.show()

        aucs[annAccession] = (geneAUC, promoterAUC)
        
outputFile = dataFolder + dataSubFolder + 'geneAndPromoterAUCs.pkl'
with open(outputFile, 'wb') as f:
    pickle.dump(aucs, f)


# get the AUC for Chromhmm files
########################################

aucs = {} # gene, promoter
clusterCount = 18
for annAccession in annAccessionList:
    
    if (annMeta[annAccession]['expressionFile'] != 'none'):

        expFileName = annMeta[annAccession]['expressionFile']
        expAccession = re.split('_|\.', expFileName)[2]
        inputFile = dataFolder + dataSubFolder + annAccession + '/geneExp_dict_' + expAccession + '.pkl'
        print(inputFile)
        with open(inputFile, 'rb') as pickledFile:
            expression = pickle.load(pickledFile)

        inputFile = dataFolder + dataSubFolder + annAccession + '/' + 'exp_promoter_labelCover_chmm.pkl'
        try:
            with open(inputFile, 'rb') as f:
                transPromoMat = pickle.load(f)  # promoter, genes
        except FileNotFoundError:
            print('no file')
        

        labelExpMat = transPromoMat['genes']
        sib = labelExpMat.sum(axis = 1)
        filterGene = sib > 0

        expArray = np.asarray([expression[x] for x in geneIDList])
        filterExp = expArray[filterGene,]
        filterExpMat = labelExpMat[filterGene, :]
        sumExp = filterExpMat.sum(axis = 1)
        filterExpMat = filterExpMat/sumExp[:,None]
    
        notExp = filterExp == 0

        #fig = plt.figure(figsize =(17, 7))
        #ax = fig.add_subplot(111)

        fig, axs = plt.subplots(nrows = 3, ncols=6, figsize =(12,9))

        geneAUC = []

        for i in range(clusterCount):
            print(i)
            filteredExpArray = filterExpMat[:, i]
            sib = np.argsort(filteredExpArray)[::-1]

            AUCCurve = np.zeros((len(sib), 2))
            xwalk = 0
            ywalk = 0
            area = 0 
            for j in range(len(sib)):
                if notExp[sib[j]] == True:
                    xwalk +=1
                    area+= ywalk
                else:
                    ywalk +=1
                    
                AUCCurve[j, 0] = xwalk
                AUCCurve[j, 1] = ywalk

            auc = area / (xwalk*ywalk)
            geneAUC.append(auc)
            
            pi = int(i /6)
            pj = i%6
            axs[pi, pj].plot(AUCCurve[:, 0], AUCCurve[:,1])
            if pj > 0:
                axs[pi, pj].set_yticks([])

            if pi < 2:
                axs[pi, pj].set_xticks([])

            clusterLabel = '%d' %(i)
            axs[pi, pj].text(xwalk*.35,ywalk/20, 'AUC: %.3f' %(auc))
            axs[pi, pj].text(xwalk*.012,ywalk*.9, '%s' %(chromLabels[i]))
        #plt.scatter(AUCCurve[:, 0], AUCCurve[:,1], s=4, ax=axs[1])
        #plt.boxplot(data)
        
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        tissue = tissue_info[annAccession]
        fig.suptitle('prediction for gene expression \n %s - %s' %(annAccession, tissue))
        fig.text(0.5, 0.04, 'zero expression', ha='center')
        fig.text(0.04, 0.5, 'expressed', va='center', rotation='vertical')

        plotFolder_add = plotFolder + annAccession + '/'
        figFile = plotFolder_add + 'geneAUC_chrom.pdf'
        print(figFile)
        plt.savefig(figFile, bbox_inches='tight')
        plt.close('all')

        #plt.show()

        labelExpMat = transPromoMat['promoter']
        sib = labelExpMat.sum(axis = 1)
        filterGene = sib > 0

        expArray = np.asarray([expression[x] for x in geneIDList])
        filterExp = expArray[filterGene,]
        filterExpMat = labelExpMat[filterGene, :]
        sumExp = filterExpMat.sum(axis = 1)
        filterExpMat = filterExpMat/sumExp[:,None]
    
        notExp = filterExp == 0

        #fig = plt.figure(figsize =(17, 7))
        #ax = fig.add_subplot(111)

        fig, axs = plt.subplots(nrows = 3, ncols=6, figsize =(12,9))

        promoterAUC = []

        for i in range(clusterCount):
            print(i)
            filteredExpArray = filterExpMat[:, i]
            sib = np.argsort(filteredExpArray)[::-1]

            AUCCurve = np.zeros((len(sib), 2))
            xwalk = 0
            ywalk = 0
            area = 0 
            for j in range(len(sib)):
                if notExp[sib[j]] == True:
                    xwalk +=1
                    area+= ywalk
                else:
                    ywalk +=1
                    
                AUCCurve[j, 0] = xwalk
                AUCCurve[j, 1] = ywalk

            auc = area / (xwalk*ywalk)

            promoterAUC.append(auc)
            
            pi = int(i /6)
            pj = i%6
            axs[pi, pj].plot(AUCCurve[:, 0], AUCCurve[:,1])
            if pj > 0:
                axs[pi, pj].set_yticks([])

            if pi < 3:
                axs[pi, pj].set_xticks([])

            clusterLabel = '%d' %(i)
            axs[pi, pj].text(xwalk*.35,ywalk/20, 'AUC: %.3f' %(auc))
            axs[pi, pj].text(xwalk*.012,ywalk*.9, '%s' %(chromLabels[i]))
        #plt.scatter(AUCCurve[:, 0], AUCCurve[:,1], s=4, ax=axs[1])
        #plt.boxplot(data)
        
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        tissue = tissue_info[annAccession]
        fig.suptitle('prediction for active promoter \n %s - %s' %(annAccession, tissue))
        fig.text(0.5, 0.04, 'zero exprssion', ha='center')
        fig.text(0.04, 0.5, 'expressed genes', va='center', rotation='vertical')

        plotFolder_add = plotFolder + annAccession + '/'
        figFile = plotFolder_add + 'promoterAUC_chrom.pdf'
        print(figFile)
        plt.savefig(figFile, bbox_inches='tight')
        plt.close('all')
        #plt.show()

        aucs[annAccession] = (geneAUC, promoterAUC)
        
outputFile = dataFolder + dataSubFolder + 'geneAndPromoterAUCs_chrom.pkl'
with open(outputFile, 'wb') as f:
    pickle.dump(aucs, f)

########################################    
# 4. get the AUC plot
########################################

# get the highest AUC for each of them, for each of them document which label had the highest AUC
inputFile = dataFolder + dataSubFolder + 'geneAndPromoterAUCs_chrom.pkl'


# get the label and the max value for each of the segway promoter and genes
inputFile = dataFolder + dataSubFolder + 'geneAndPromoterAUCs.pkl'
with open(inputFile, 'rb') as f:
    aucs = pickle.load(f)

segwayGenesAUC = {} # keeps the highest auc label and value
segwayPromoAUC = {} # keeps the highest auc label and value
segway3genes = {}
segway3promo = {}
for annAccession in annAccessionList:
    if annAccession in list(aucs.keys()):
        annauc = aucs[annAccession]

        label_term_mapping = {}
        mnemonics_file = dataFolder + dataSubFolder + annAccession + '/' + 'mnemonics_v04.txt'
        with open(mnemonics_file, 'r') as mnemonics:
            for line in mnemonics:
                print(line)
                label = line.strip().split()[0]
                term = line.strip().split()[1]
                label_term_mapping[label] = term

        genesAUC = annauc[0]
        promoterAUC = annauc[1]

        genesind = np.argmax(genesAUC)
        promoind = np.argmax(promoterAUC)

        genesval = genesAUC[genesind]
        promoval = promoterAUC[promoind]

        gl = label_term_mapping[str(genesind)]
        pl = label_term_mapping[str(promoind)]

        segwayGenesAUC[annAccession] = (genesind, gl, genesval)
        segwayPromoAUC[annAccession] = (promoind, pl, promoval)

        # get the three highest values above 70%

        sortedGenes = np.argsort(genesAUC)[::-1]
        sortedPromo = np.argsort(promoterAUC)[::-1]

        geneLabels = [label_term_mapping[str(i)] for i in sortedGenes[0:3]]
        promoLabels = [label_term_mapping[str(i)] for i in sortedPromo[0:3]]

        topGeneAUC = [genesAUC[i] for i in sortedGenes[0:3]]
        topPromoAUC = [promoterAUC[i] for i in sortedPromo[0:3]]

        segway3genes[annAccession] = (sortedGenes[0:3], geneLabels, topGeneAUC)
        segway3promo[annAccession] = (sortedPromo[0:3], promoLabels, topPromoAUC)


# get the highest AUC for each of them, for each of them document which label had the highest AUC
inputFile = dataFolder + dataSubFolder + 'geneAndPromoterAUCs_chrom.pkl'
with open(inputFile, 'rb') as f:
    aucs = pickle.load(f)

chromGenesAUC = {} # keeps the highest auc label and value
chromPromoAUC = {} # keeps the highest auc label and value
chrom3genes = {}
chrom3promo = {}
for annAccession in annAccessionList:
    if annAccession in list(aucs.keys()):
        annauc = aucs[annAccession]

        genesAUC = annauc[0]
        promoterAUC = annauc[1]

        genesind = np.argmax(genesAUC)
        promoind = np.argmax(promoterAUC)

        genesval = genesAUC[genesind]
        promoval = promoterAUC[promoind]

        chromGenesAUC[annAccession] = (chromLabels[genesind], genesval)
        chromPromoAUC[annAccession] = (chromLabels[promoind], promoval)

        # get the three highest values above 70%

        sortedGenes = np.argsort(genesAUC)[::-1]
        sortedPromo = np.argsort(promoterAUC)[::-1]

        geneLabels = [chromLabels[i] for i in sortedGenes[0:3]]
        promoLabels = [chromLabels[i] for i in sortedPromo[0:3]]

        topGeneAUC = [genesAUC[i] for i in sortedGenes[0:3]]
        topPromoAUC = [promoterAUC[i] for i in sortedGenes[0:3]]

        chrom3genes[annAccession] = (sortedGenes[0:3], geneLabels, topGeneAUC)
        chrom3promo[annAccession] = (sortedPromo[0:3], promoLabels, topPromoAUC)


# aaaaand we are done with it.

# colors:
colorDict_chrom = {'TssA':[255,0,0], 'TssFlnk': [255,69,0], 'TssFlnkU'	: [255,69,0], 
	           'TssFlnkD':[255,69,0],'Tx':[0,128,0],'TxWk':[0,100,0],'EnhG1':[194,225,5],                   
                   'EnhG2':[194,225,5],'EnhA1':[255,195,77],'EnhA2':[255,195,77],
                   'EnhWk':[255,255,0],'ZNF/Rpts':[102,205,170],'Het':[138,145,208],
                   'TssBiv':[205,92,92],'EnhBiv' :[189,183,107],'ReprPC' :[128,128,128],
                   'ReprPCWk':[192,192,192],'Quies':[255,255,255]}

colorDict_segway = {'Transcribed': [0,128,0], 'Enhancer':[194,225,5], 'EnhancerLow': [255,255,0],
                    'Promoter':[255,0,0], 'PromoterFlanking':[255,69,0], 'K9K36':[102,205,170],
                    'Bivalent':[189,183,107], 'CTCF':[0,0,120], 'FacultativeHet':[138,145,208],
                    'ConstitutiveHet':[138,145,208], 'Quiescent':[20,20,20]}


#segway3promo
allkeys = list(segway3promo.keys())
allhighVal = [segway3promo[allkeys[i]][2][0] for i in range(len(allkeys))]
asHighVal = np.argsort(allhighVal)[::-1]
segpromoHV = np.sort(allhighVal)[::-1]
plt.plot(segpromoHV)
plt.show()

# walk on keys, in the order of asHighVal. for each key, plot the 3 scatters with the relevant color
fig = plt.figure(figsize =(17, 7))
ax = fig.add_subplot(111)

for i,ind in enumerate(asHighVal):
    annAccession = allkeys[ind]
    values = segway3promo[annAccession][2]
    labels = segway3promo[annAccession][1]
    colors = [colorDict_segway[l] for l in labels]

    c = np.asarray(colors[0])/255
    plt.scatter(i, values[0], color=c)
    if values[1] > .7:
        c = np.asarray(colors[1])/255
        plt.scatter(i, values[1], color=c)
    if values[2] > .7:
        c = np.asarray(colors[2])/255
        plt.scatter(i, values[2], color=c)
plt.title('segway promoter predictions')
plt.show()


#segway3genes
allkeys = list(segway3genes.keys())
allhighVal = [segway3genes[allkeys[i]][2][0] for i in range(len(allkeys))]
asHighVal = np.argsort(allhighVal)[::-1]

# walk on keys, in the order of asHighVal. for each key, plot the 3 scatters with the relevant color
fig = plt.figure(figsize =(17, 7))
ax = fig.add_subplot(111)

for i,ind in enumerate(asHighVal):
    annAccession = allkeys[ind]
    values = segway3genes[annAccession][2]
    labels = segway3genes[annAccession][1]
    colors = [colorDict_segway[l] for l in labels]

    c = np.asarray(colors[0])/255
    plt.scatter(i, values[0], color=c)
    if values[1] > .7:
        c = np.asarray(colors[1])/255
        plt.scatter(i, values[1], color=c)
    if values[2] > .7:
        c = np.asarray(colors[2])/255
        plt.scatter(i, values[2], color=c)

plt.title('segway gene prediction')
plt.show()

    
        #fig = plt.figure(figsize =(17, 7))
        #ax = fig.add_subplot(111)
    # plot the plot
    

# sort by the max value
# plot by that
#chrom3promo
allkeys = list(chrom3promo.keys())
allhighVal = [chrom3promo[allkeys[i]][2][0] for i in range(len(allkeys))]
asHighVal = np.argsort(allhighVal)[::-1]
chrompromoHV = np.sort(allhighVal)[::-1]
plt.plot(chrompromoHV)
plt.show()

# walk on keys, in the order of asHighVal. for each key, plot the 3 scatters with the relevant color
fig = plt.figure(figsize =(17, 7))
ax = fig.add_subplot(111)

for i,ind in enumerate(asHighVal):
    annAccession = allkeys[ind]
    values = chrom3promo[annAccession][2]
    labels = chrom3promo[annAccession][1]
    colors = [colorDict_chrom[l] for l in labels]

    c = np.asarray(colors[0])/255
    plt.scatter(i, values[0], color=c)
    if values[1] > .7:
        c = np.asarray(colors[1])/255
        plt.scatter(i, values[1], color=c)
    if values[2] > .7:
        c = np.asarray(colors[2])/255
        plt.scatter(i, values[2], color=c)
plt.title('chrom promo prediction')
plt.show()

#chrom3genes
allkeys = list(chrom3genes.keys())
allhighVal = [chrom3genes[allkeys[i]][2][0] for i in range(len(allkeys))]
asHighVal = np.argsort(allhighVal)[::-1]

# walk on keys, in the order of asHighVal. for each key, plot the 3 scatters with the relevant color
fig = plt.figure(figsize =(17, 7))
ax = fig.add_subplot(111)

for i,ind in enumerate(asHighVal):
    annAccession = allkeys[ind]
    values = chrom3genes[annAccession][2]
    labels = chrom3genes[annAccession][1]
    colors = [colorDict_chrom[l] for l in labels]

    c = np.asarray(colors[0])/255
    plt.scatter(i, values[0], color=c)
    if values[1] > .7:
        c = np.asarray(colors[1])/255
        plt.scatter(i, values[1], color=c)
    if values[2] > .7:
        c = np.asarray(colors[2])/255
        plt.scatter(i, values[2], color=c)
plt.title('chrom gene prediction')
plt.show()


# for active promoter how t
allkeys = list(segway3promo.keys())
allhighVal = [segway3promo[allkeys[i]][2][0] for i in range(len(allkeys))]
asHighVal = np.argsort(allhighVal)[::-1]
segpromoHV = np.sort(allhighVal)[::-1]
plt.plot(segpromoHV)

allkeys = list(chrom3promo.keys())
allhighVal = [chrom3promo[allkeys[i]][2][0] for i in range(len(allkeys))]
asHighVal = np.argsort(allhighVal)[::-1]
chrompromoHV = np.sort(allhighVal)[::-1]
plt.plot(chrompromoHV)

plt.show()

# for genes how they do
allkeys = list(segway3genes.keys())
allhighVal = [segway3genes[allkeys[i]][2][0] for i in range(len(allkeys))]
asHighVal = np.argsort(allhighVal)[::-1]
seggenesHV = np.sort(allhighVal)[::-1]
plt.plot(seggenesHV)

allkeys = list(chrom3genes.keys())
allhighVal = [chrom3genes[allkeys[i]][2][0] for i in range(len(allkeys))]
asHighVal = np.argsort(allhighVal)[::-1]
chromgenesHV = np.sort(allhighVal)[::-1]
plt.plot(chromgenesHV)
plt.show()

data = []
data.append(seggenesHV)
data.append(chromgenesHV)

plt.boxplot(data)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)

plt.show()
from scipy import stats
stats.ttest_ind(chromgenesHV, seggenesHV)

 
### Draft 
########################################
# get the distribution of each of the samples ratio for each of the labels for exp vs non-exp genes

        book = pd.DataFrame(clusterMats[2])
        sns.heatmap(book)
        figFile = dataFolder + dataSubFolder + annAccession + '/exp2cluster_heatmap.pdf'
        plt.savefig(figFile)
        plt.close('all')
        book = pd.DataFrame(classMats[2])
        sns.heatmap(book)
        figFile = dataFolder + dataSubFolder + annAccession + '/exp2class_heatmap.pdf'
        plt.savefig(figFile)
        plt.close('all')
        book = pd.DataFrame(clusterMats[0])
        sns.heatmap(book)
        figFile = dataFolder + dataSubFolder + annAccession + '/exp0cluster_heatmap.pdf'
        plt.savefig(figFile)
        plt.close('all')
        book = pd.DataFrame(classMats[0])
        sns.heatmap(book)
        figFile = dataFolder + dataSubFolder + annAccession + '/exp0class_heatmap.pdf'
        plt.savefig(figFile)
        plt.close('all')


        print('annotation summary saved in %s' %(outputFile))


# For the samples that have transcriptomic data:

# Get the prediction for the expressed genes

# Get the prediciton for the expressed genes. from chromhmm 

# How well do they do regarding the transcription data? active prmoter


# How well do they do regarding the transcription data? transcribed region etc.
