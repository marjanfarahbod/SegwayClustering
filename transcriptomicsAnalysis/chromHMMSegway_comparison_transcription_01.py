# how do they compare regarding the transcription data.
# How well do they do regarding the transcription data? active prmoter
# How well do they do regarding the transcription data? transcribed region etc.
#
# Note: expressed genes act very much the same between the low and high expression.
# Note: For each gene, get the ratio of its promoter region and its region covered with each of the labels
# Note: Examine the distribution of these labels between the expressed and non expressed genes.
# 
# 0. Initials
# 0.0 The function for extracting data - Segway
# 0.1 The function for extracting data - chromhmm
# * Segway data extraction - function for 0.0 in transcription_overlap.py
# * CHromHMM data extraction -  function for 0.1 in transcription_overlap.py
# 1. The AUC curve for genes which were included
# 2. get the AUC plot
# 3. AUCs for the transcribed and active promoter regions
# 4. The AUCs box plots 
# 5. Regression: see code chromHMMSegway_comparison_transcription_regression.py
# 6. Regression plots: see code chromHMMSegway_comparison_transcription_regression.py

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
from QC_transcriptionComparison_util import Gene, Exon, Annotation, AnnotationClass

# Segway states:
segwayStates = ['Enhancer', 'EnhancerLow', 'Promoter', 'PromoterFlanking', 'Transcribed', 'CTCF', 'K9K36', 'Bivalent', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent']
print(len(segwayStates))

book =  'EnhA1 EnhA2 EnhBiv EnhG1 EnhG2 EnhWk Het Quies ReprPC ReprPCWk TssA TssBiv TssFlnk TssFlnkD TssFlnkU Tx TxWk ZNF/Rpts'
chromLabels = book.split()
print(len(chromLabels))

# General data folder
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'

# GTF data structure
fileName = dataFolder + '/geneLists.pkl'
with open(fileName, 'rb') as pickleFile:
    geneListsAndIDs = pickle.load(pickleFile)
    
geneList = geneListsAndIDs[0]
geneIDList = geneListsAndIDs[1]
del geneListsAndIDs

# load the meta file 
inputFileName = 'all235Annot_meta_corrected.pkl'
#inputFileName = 'all235Annot_meta.pkl'
inputFile = dataFolder +  inputFileName
with open(inputFile, 'rb') as f:
    allMeta = pickle.load(f)

accessionList = list(allMeta.keys())

# How well do they do regarding the transcription data? active prmoter
# How well do they do regarding the transcription data? transcribed region etc.

# expressed genes act very much the same between the low and high expression.

# for each gene, get the ratio of its promoter region and its region covered with each of the labels
# examine the distribution of these labels between the expressed and non expressed genes.
########################################
# 000 why are the genes with no labels? 

sib = labelExpMat.sum(axis = 1) # is this the same for all Segway files? 19339
filterGene = sib > 0 # genes included in the annotation
count = 0
for cgi in range(26017):

    'we are in the gene territory, we are walking on the genes'
        
    geneID = geneIDList[cgi]
    gene_chr = geneList[geneIDList[cgi]].chrom
        
    strand = geneList[geneID].strand
    gene_start = geneList[geneID].start
    gene_end = geneList[geneID].end

    gene_length = gene_end - gene_start - 300
    if gene_length <= 0:
        count+=1
        

########################################
# 0.0 The function - Segway
########################################
from transcription_overlap import EPLabelCover_Segway

'''
Seway bed files are already sorted and filtered. Please see the QC_transcriptionComparison_03.py section 5
output is saved in : annotationFolder + 'exp_promoter_labelCover_promLength%d.pkl' %(promLength)
'''

# load the meta file 
inputFileName = 'all235Annot_meta_corrected.pkl'
#inputFileName = 'all235Annot_meta.pkl'
inputFile = dataFolder +  inputFileName
with open(inputFile, 'rb') as f:
    allMeta = pickle.load(f)

accessionList = list(allMeta.keys())
accessionList.remove(accessionList[205])

count = 0
for accession in accessionList:
    annotation = allMeta[accession]
    count+=1

    #if ((('38batch' in annotation['folder']) or ('May11' in annotation['folder'])) and not(annotation['RNAseqFile'] == 'none')):
    if ('testBatch105' in annotation['folder']) and not(annotation['RNAseqFile'] == 'none') and not(annotation['chromFile'] == 'none'):
        #RNAFile = annotation['RNAseqFile'][0]

        if len(annotation['RNAseqFile'][0]) > 1:
            RNAFile = annotation['RNAseqFile'][0]
        else:
            RNAFile = annotation['RNAseqFile']
        print(RNAFile)

        print(count)
        print(accession)
        print(RNAFile)

        annotationFolder = annotation['folder']
        print(annotationFolder)

        expAccession = RNAFile[-15:-4]
        expFile = annotationFolder +  'geneExp_dict_' + expAccession + '.pkl'
        annFile = annotation['bedFile']
        mnemFile = annotationFolder + 'mnemonics_v04.txt'
        #geneList
        #geneIDList =
        promLength = 3000
        tssDown = 0
        EPLabelCover_Segway(annotationFolder, expFile, annFile, mnemFile, geneList, geneIDList, promLength, tssDown)

########################################
# 0.1 The function - chromhmm
########################################

'''
output is saved in annotationFolder + 'exp_promoter_labelCover_chmm.pkl'
'''

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

count = 0
for accession in accessionList:
    annotation = allMeta[accession]

    #if (('38batch' in annotation['folder']) or ('May11' in annotation['folder'])) and not(annotation['RNAseqFile'] == 'none') and not(annotation['chromFile'] == 'none'):
    if not(annotation['RNAseqFile'] == 'none') and not(annotation['chromFile'] == 'none'):
        count+=1
        print(count)
        if len(annotation['RNAseqFile'][0]) > 1:
            RNAFile = annotation['RNAseqFile'][0]
        else:
            RNAFile = annotation['RNAseqFile']
        print(RNAFile)
        print(accession)

        annotationFolder = annotation['folder']
        print(annotationFolder)

        expAccession = RNAFile[-15:-4]
        expFile = annotationFolder +  'geneExp_dict_' + expAccession + '.pkl'
        chmmFileName = annotation['chromFile'].split('/')[-1]

        if (('38batch' in annotation['folder']) or ('May11' in annotation['folder'])):
            annFile = annotationFolder + 'sorted_' + chmmFileName
        else:
            annFile = annotationFolder + chmmFileName
            
        #geneList
        #geneIDList =
        promLength = 3000

        print(expAccession)
        print(expFile)
        print(annFile)

        EPLabelCover_chromHMM(annotationFolder, expFile, annFile, chromLabels, geneList, geneIDList, promLength)
        # chromhmm has missing regions in the genome. Its coverage is ...

########################################        
# 1. The AUC curve for genes which were included
########################################

'''
update all the plots 
'''

inputFileName = 'all235Annot_meta_corrected.pkl'
#inputFileName = 'all235Annot_meta.pkl'
inputFile = dataFolder +  inputFileName
with open(inputFile, 'rb') as f:
    allMeta = pickle.load(f)
    
accessionList = list(allMeta.keys())

aucs = {} # gene, promoter
count = 0
for accession in accessionList:
    annotation = allMeta[accession]
    print(annotation)
    print(count)
    count +=1
    
    #if (annMeta[annAccession]['expressionFile'] != 'none'): # for samples with expression data
    if not(annotation['RNAseqFile'] == 'none') and not(annotation['chromFile'] == 'none'):

        print(accession)
        annotationFolder = annotation['folder']
        print(annotationFolder)

        if len(annotation['RNAseqFile'][0]) > 1:
            RNAFile = annotation['RNAseqFile'][0]
        else:
            RNAFile = annotation['RNAseqFile']
        print(RNAFile)

        expAccession = RNAFile[-15:-4]
        expFile = annotationFolder +  'geneExp_dict_' + expAccession + '.pkl'

        with open(expFile, 'rb') as pickledFile:
            expression = pickle.load(pickledFile) # load the expression file

        # load the mnemonics
        label_term_mapping = {}
        mnemonics_file = annotationFolder + 'mnemonics_v04.txt'
        with open(mnemonics_file, 'r') as mnemonics:
            for line in mnemonics:
                print(line)
                label = line.strip().split()[0]
                term = line.strip().split()[1]
                label_term_mapping[label] = term

        clusterCount = len(label_term_mapping)

        # load the recorded table
        inputFile = annotationFolder + 'exp_promoter_labelCover_promLength3000.pkl'
        with open(inputFile, 'rb') as f:
            transPromoMat = pickle.load(f)  # promoter, genes

        labelExpMat = transPromoMat['genes'] # matrix of genes by labels. How much of the gene body was covered by the label
        sib = labelExpMat.sum(axis = 1)
        filterGene = sib > 0 # genes included in the annotation

        #c = 0
        #for id in geneIDList:
        #    if not(id in book):
        #        c+=1

        expArray = np.asarray([expression[x] for x in geneIDList]) # dict to array
        filterExp = expArray[filterGene,] # getting the expression for genes included in the annotation 
        filterExpMat = labelExpMat[filterGene, :] # getting the label coverage for the same group of genes
        sumExp = filterExpMat.sum(axis = 1) 
        filterExpMat = filterExpMat/sumExp[:,None] # normalizing the gene length by the label coverage
    
        notExp = filterExp == 0 # selecting the genes with zero expression

        #fig = plt.figure(figsize =(17, 7))
        #ax = fig.add_subplot(111)

        # AUC plot for each of the labels
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
        tissue =  '%s  %s  %s' %(annotation['tissueInfo'][0], annotation['tissueInfo'][1], annotation['tissueInfo'][2])
        fig.suptitle('prediction for gene expression \n %s - %s' %(accession, tissue))
        fig.text(0.5, 0.04, 'zero expression', ha='center')
        fig.text(0.04, 0.5, 'expressed', va='center', rotation='vertical')

        plotFolder_add = annotation['plotFolder']
        figFile = plotFolder_add + 'geneAUC_v04.pdf'
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
        tissue =  '%s  %s  %s' %(annotation['tissueInfo'][0], annotation['tissueInfo'][1], annotation['tissueInfo'][2])
        fig.suptitle('prediction for gene expression \n %s - %s' %(accession, tissue))
        fig.text(0.5, 0.04, 'zero expression', ha='center')
        fig.text(0.04, 0.5, 'expressed', va='center', rotation='vertical')

        plotFolder_add = annotation['plotFolder']
        figFile = plotFolder_add + 'promoterAUC_v04_3000bp.pdf'
        print(figFile)
        plt.savefig(figFile, bbox_inches='tight')
        plt.close('all')
        #plt.show()

        aucs[accession] = (geneAUC, promoterAUC)
        
outputFile = dataFolder + 'geneAndPromoterAUCs_3000_segway_v04.pkl'
with open(outputFile, 'wb') as f:
    pickle.dump(aucs, f)


# get the AUC for Chromhmm files
########################################

aucs = {} # gene, promoter
clusterCount = 18
for accession in accessionList:
    annotation = allMeta[accession]
    
    #if (annMeta[annAccession]['expressionFile'] != 'none'): # for samples with expression data
    if not(annotation['RNAseqFile'] == 'none') and not(annotation['chromFile'] == 'none'):

        print(accession)
        annotationFolder = annotation['folder']
        print(annotationFolder)

        if len(annotation['RNAseqFile'][0]) > 1:
            RNAFile = annotation['RNAseqFile'][0]
        else:
            RNAFile = annotation['RNAseqFile']
        print(RNAFile)

        expAccession = RNAFile[-15:-4]
        expFile = annotationFolder +  'geneExp_dict_' + expAccession + '.pkl'

        with open(expFile, 'rb') as pickledFile:
            expression = pickle.load(pickledFile) # load the expression file

        inputFile = annotationFolder + 'exp_promoter_labelCover_chmm_promLength3000.pkl'
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
        tissue =  '%s  %s  %s' %(annotation['tissueInfo'][0], annotation['tissueInfo'][1], annotation['tissueInfo'][2])
        fig.suptitle('prediction for gene expression \n %s - %s' %(accession, tissue))
        fig.text(0.5, 0.04, 'zero expression', ha='center')
        fig.text(0.04, 0.5, 'expressed', va='center', rotation='vertical')

        plotFolder_add = annotation['plotFolder']
        figFile = plotFolder_add + 'geneAUC_chrom_v04.pdf'
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
        tissue =  '%s  %s  %s' %(annotation['tissueInfo'][0], annotation['tissueInfo'][1], annotation['tissueInfo'][2])
        fig.suptitle('prediction for gene expression \n %s - %s' %(accession, tissue))
        fig.text(0.5, 0.04, 'zero expression', ha='center')
        fig.text(0.04, 0.5, 'expressed', va='center', rotation='vertical')

        plotFolder_add = annotation['plotFolder']
        figFile = plotFolder_add + 'promoterAUC_chrom_v04_3000bp.pdf'
        print(figFile)
        plt.savefig(figFile, bbox_inches='tight')
        plt.close('all')
        #plt.show()

        aucs[accession] = (geneAUC, promoterAUC)

outputFile = dataFolder + 'geneAndPromoterAUCs_3000_chrom.pkl'
with open(outputFile, 'wb') as f:
    pickle.dump(aucs, f)

########################################    
# 2. get the AUC plot
########################################

# 2.0 Get the highest label and values - for the scatter plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# get the highest AUC for each of them, for each of them document which label had the highest AUC
inputFile = dataFolder + 'geneAndPromoterAUCs_3000_chrom.pkl'

# get the label and the max value for each of the segway promoter and genes
inputFile = dataFolder + 'geneAndPromoterAUCs_3000_segway_v04.pkl'
with open(inputFile, 'rb') as f:
    aucs = pickle.load(f)

segwayGenesAUC = {} # keeps the highest auc label and value
segwayPromoAUC = {} # keeps the highest auc label and value
segway3genes = {}
segway3promo = {}
for accession in accessionList:
    if accession in list(aucs.keys()):
        annotation = allMeta[accession]
        
        annotationFolder = annotation['folder']
        annauc = aucs[accession]

        label_term_mapping = {}
        mnemonics_file = annotationFolder + 'mnemonics_v04.txt'
        with open(mnemonics_file, 'r') as mnemonics:
            for line in mnemonics:
                print(line)
                label = line.strip().split()[0]
                term = line.strip().split()[1]
                label_term_mapping[label] = term

        genesAUC = annauc[0] # auc for all labels for the sample - genes
        promoterAUC = annauc[1] # auc for all labels for the sample - promoter

        genesind = np.argmax(genesAUC) # label with the highest auc for the genes
        promoind = np.argmax(promoterAUC) # label with the highest auc for the promoters

        genesval = genesAUC[genesind] # that auc value for genes
        promoval = promoterAUC[promoind] # that auc value for the genes

        gl = label_term_mapping[str(genesind)]
        pl = label_term_mapping[str(promoind)]

        segwayGenesAUC[accession] = (genesind, gl, genesval)
        segwayPromoAUC[accession] = (promoind, pl, promoval)

        # get the three highest values above 70% - everything above, but for the first three

        sortedGenes = np.argsort(genesAUC)[::-1]
        sortedPromo = np.argsort(promoterAUC)[::-1]

        geneLabels = [label_term_mapping[str(i)] for i in sortedGenes[0:3]]
        promoLabels = [label_term_mapping[str(i)] for i in sortedPromo[0:3]]

        topGeneAUC = [genesAUC[i] for i in sortedGenes[0:3]]
        topPromoAUC = [promoterAUC[i] for i in sortedPromo[0:3]]

        segway3genes[accession] = (sortedGenes[0:3], geneLabels, topGeneAUC)
        segway3promo[accession] = (sortedPromo[0:3], promoLabels, topPromoAUC)

# same thing above, but for chromhmm
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# get the highest AUC for each of them, for each of them document which label had the highest AUC
inputFile = dataFolder +  'geneAndPromoterAUCs_3000_chrom.pkl'
with open(inputFile, 'rb') as f:
    aucs = pickle.load(f)

chromGenesAUC = {} # keeps the highest auc label and value
chromPromoAUC = {} # keeps the highest auc label and value
chrom3genes = {}
chrom3promo = {}
for accession in accessionList:
    if accession in list(aucs.keys()): 
        annauc = aucs[accession]

        genesAUC = annauc[0]
        promoterAUC = annauc[1]

        genesind = np.argmax(genesAUC)
        promoind = np.argmax(promoterAUC)

        genesval = genesAUC[genesind]
        promoval = promoterAUC[promoind]

        chromGenesAUC[accession] = (chromLabels[genesind], genesval)
        chromPromoAUC[accession] = (chromLabels[promoind], promoval)

        # get the three highest values above 70%

        sortedGenes = np.argsort(genesAUC)[::-1]
        sortedPromo = np.argsort(promoterAUC)[::-1]

        geneLabels = [chromLabels[i] for i in sortedGenes[0:3]]
        promoLabels = [chromLabels[i] for i in sortedPromo[0:3]]

        topGeneAUC = [genesAUC[i] for i in sortedGenes[0:3]]
        topPromoAUC = [promoterAUC[i] for i in sortedGenes[0:3]]

        chrom3genes[accession] = (sortedGenes[0:3], geneLabels, topGeneAUC)
        chrom3promo[accession] = (sortedPromo[0:3], promoLabels, topPromoAUC)


# aaaaand we are done with it.

# 2.1 The scatter plot, there are four of them. 
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# >>>>>>>>>>>>>>>
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


# >>>>>>>>>>>>>>> plot 1/4: segway3promo
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


# >>>>>>>>>>>>>>> plot 2/4: #segway3genes
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
# >>>>>>>>>>>>>>> plot 3/4: chrom3promo
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

# >>>>>>>>>>>>>>> plot 4/4: chrom3genes
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


########################################
# 3. AUCs for the transcribed and active promoter regions
########################################

'''
This is where I sum them up. The values are the answer to question: how well do annotations predict transcribed and active promoter regions? For this, I sum up all Segway transcription labels, for Chrom I sum up the txweak and tx. For the promoter regions, from Segway I sum up all promoter and from ChromHMM I sum up TSSA and a couple of flanks. 
Some Segway samples might not have transcribed regions and I will skip them. 
'''

inputFileName = 'all235Annot_meta_corrected.pkl'
#inputFileName = 'all235Annot_meta.pkl'
inputFile = dataFolder +  inputFileName
with open(inputFile, 'rb') as f:
    allMeta = pickle.load(f)
    
accessionList = list(allMeta.keys())

aucs = {} # gene, promoter
count = 0
for accession in accessionList:
    annotation = allMeta[accession]
    print(annotation)
    print(count)
    count +=1
    
    #if (annMeta[annAccession]['expressionFile'] != 'none'): # for samples with expression data
    if not(annotation['RNAseqFile'] == 'none') and not(annotation['chromFile'] == 'none'):

        print(accession)
        annotationFolder = annotation['folder']
        print(annotationFolder)

        if len(annotation['RNAseqFile'][0]) > 1:
            RNAFile = annotation['RNAseqFile'][0]
        else:
            RNAFile = annotation['RNAseqFile']
        print(RNAFile)

        expAccession = RNAFile[-15:-4]
        expFile = annotationFolder +  'geneExp_dict_' + expAccession + '.pkl'

        with open(expFile, 'rb') as pickledFile:
            expression = pickle.load(pickledFile) # load the expression file

        # load the mnemonics
        label_term_mapping = {}
        mnemonics_file = annotationFolder + 'mnemonics_v04.txt'
        with open(mnemonics_file, 'r') as mnemonics:
            for line in mnemonics:
                print(line)
                label = line.strip().split()[0]
                term = line.strip().split()[1]
                label_term_mapping[label] = term

        clusterCount = len(label_term_mapping)

        # load the recorded table
        inputFile = annotationFolder + 'exp_promoter_labelCover_promLength3000.pkl'
        with open(inputFile, 'rb') as f:
            transPromoMat = pickle.load(f)  # promoter, genes

        labelExpMat = transPromoMat['genes'] # matrix of genes by labels. How much of the gene body was covered by the label
        sib = labelExpMat.sum(axis = 1)
        filterGene = sib > 0 # genes included in the annotation

        #c = 0
        #for id in geneIDList:
        #    if not(id in book):
        #        c+=1

        expArray = np.asarray([expression[x] for x in geneIDList]) # dict to array
        filterExp = expArray[filterGene,] # getting the expression for genes included in the annotation 
        filterExpMat = labelExpMat[filterGene, :] # getting the label coverage for the same group of genes
        sumExp = filterExpMat.sum(axis = 1) 
        filterExpMat = filterExpMat/sumExp[:,None] # normalizing the gene length by the label coverage
    
        notExp = filterExp == 0 # selecting the genes with zero expression

        #fig = plt.figure(figsize =(17, 7))
        #ax = fig.add_subplot(111)

        # AUC plot for each of the labels
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # fig, axs = plt.subplots(nrows = 4, ncols=4, figsize =(10,9))

        # get the clusters for which labels are transcribed
        # the gene
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        tLabels = []
        for i in range(clusterCount):
            if label_term_mapping[str(i)] == 'Transcribed':
                tLabels.append(i)

        if len(tLabels) == 0: 
            continue

        geneAUC = []

        filteredExpArray = np.sum(filterExpMat[:, tLabels], axis = 1)
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

        # the promoter
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        labelExpMat = transPromoMat['promoter']
        sib = labelExpMat.sum(axis = 1)
        filterGene = sib > 0

        expArray = np.asarray([expression[x] for x in geneIDList])
        filterExp = expArray[filterGene,]
        filterExpMat = labelExpMat[filterGene, :]
        sumExp = filterExpMat.sum(axis = 1)
        filterExpMat = filterExpMat/sumExp[:,None]
    
        notExp = filterExp == 0

        pLabels = []
        for i in range(clusterCount):
            if label_term_mapping[str(i)] == 'Promoter':
                pLabels.append(i)

        if len(tLabels) == 0: 
            continue

        filteredExpArray = np.sum(filterExpMat[:, pLabels], axis = 1)
        sib = np.argsort(filteredExpArray)[::-1]
        
        promoterAUC = []

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
            
        aucs[accession] = (geneAUC, promoterAUC)
        
outputFile = dataFolder + 'geneAndPromoterAUCs_sumAllTransPromo_3000_segway_v04.pkl'
with open(outputFile, 'wb') as f:
    pickle.dump(aucs, f)


# get the AUC for Chromhmm files
########################################
print(chromLabels)

aucs = {} # gene, promoter
clusterCount = 18
for accession in accessionList:
    annotation = allMeta[accession]
    
    #if (annMeta[annAccession]['expressionFile'] != 'none'): # for samples with expression data
    if not(annotation['RNAseqFile'] == 'none') and not(annotation['chromFile'] == 'none'):

        print(accession)
        annotationFolder = annotation['folder']
        print(annotationFolder)

        if len(annotation['RNAseqFile'][0]) > 1:
            RNAFile = annotation['RNAseqFile'][0]
        else:
            RNAFile = annotation['RNAseqFile']
        print(RNAFile)

        expAccession = RNAFile[-15:-4]
        expFile = annotationFolder +  'geneExp_dict_' + expAccession + '.pkl'

        with open(expFile, 'rb') as pickledFile:
            expression = pickle.load(pickledFile) # load the expression file

        inputFile = annotationFolder + 'exp_promoter_labelCover_chmm_promLength3000.pkl'
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

        tLabels = [15, 16]
        filteredExpArray = np.sum(filterExpMat[:, tLabels], axis = 1)
        sib = np.argsort(filteredExpArray)[::-1]
        
        geneAUC = []

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
            
        # promoter
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
        labelExpMat = transPromoMat['promoter']
        sib = labelExpMat.sum(axis = 1)
        filterGene = sib > 0

        expArray = np.asarray([expression[x] for x in geneIDList])
        filterExp = expArray[filterGene,]
        filterExpMat = labelExpMat[filterGene, :]
        sumExp = filterExpMat.sum(axis = 1)
        filterExpMat = filterExpMat/sumExp[:,None]
    
        notExp = filterExp == 0

        pLabels = [10, 12]
        filteredExpArray = np.sum(filterExpMat[:, pLabels], axis = 1)
        sib = np.argsort(filteredExpArray)[::-1]

        promoterAUC = []

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
            
        aucs[accession] = (geneAUC, promoterAUC)

outputFile = dataFolder + 'geneAndPromoterAUCs_sumTxTxweak_TssATssFlank_3000_chrom.pkl'
with open(outputFile, 'wb') as f:
    pickle.dump(aucs, f)

########################################
# 4. plots (the two box plots)
########################################
inputFile = dataFolder + 'geneAndPromoterAUCs_sumTxTxweak_TssATssFlank_3000_chrom.pkl'
with open(inputFile, 'rb') as f:
    chromAucs = pickle.load(f)

inputFile = dataFolder + 'geneAndPromoterAUCs_sumAllTransPromo_3000_segway_v04.pkl'
with open(inputFile, 'rb') as f:
    segAucs = pickle.load(f)

accessionList = list(segAucs.keys())

promoMat = np.zeros((88, 2))
geneMat = np.zeros((88, 2))
#chromAucs = aucs
for i, accession in enumerate(accessionList):
    # get promo matrix

    promoMat[i, 0] = segAucs[accession][0][0]
    promoMat[i, 1] = chromAucs[accession][0][0]

    
    geneMat[i, 0] = segAucs[accession][1][0]
    geneMat[i, 1] = chromAucs[accession][1][0]


plt.boxplot(promoMat)
plt.show()


 
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
