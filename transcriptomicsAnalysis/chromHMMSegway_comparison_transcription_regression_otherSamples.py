
########################################
# 7.0 That regressiong thing - initial
########################################

'''
Segway needs summing up of the states. Chromhmm has already the same count of states
'''
from numpy import random
import pickle
import re
import numpy as np
import pandas as pd
import subprocess as sp
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
from sklearn.metrics import classification_report
from sklearn.metrics import accuracy_score

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'

# load the meta file 
inputFileName = 'all235Annot_meta_corrected.pkl'
#inputFileName = 'all235Annot_meta.pkl'
inputFile = dataFolder +  inputFileName
with open(inputFile, 'rb') as f:
    allMeta = pickle.load(f)

# GTF data structure
fileName = dataFolder + '/geneLists.pkl'
with open(fileName, 'rb') as pickleFile:
    geneListsAndIDs = pickle.load(pickleFile)
    
geneList = geneListsAndIDs[0]
geneIDList = geneListsAndIDs[1][0:24771]
del geneListsAndIDs

# Segway states:
segwayStates = ['Enhancer', 'EnhancerLow', 'Promoter', 'PromoterFlanking', 'Transcribed', 'CTCF', 'K9K36', 'Bivalent', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent']
print(len(segwayStates))
segwayStateCount = len(segwayStates)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#define a filter for zero genes. Chrom has much more zero genes cause it doesnt have coverage for some part of the genome.

# for samples with both chrom and RNAseq, load the file and get the overlapping genes. We can get a list of overlapping genes for each of the samples. 

count = 0
geneFilterAll = {}
for accession in accessionList:
    annotation = allMeta[accession]

    annotationFolder = annotation['folder']

    if not(annotation['RNAseqFile'] == 'none') and not(annotation['chromFile'] == 'none'):
        count+=1

        #if count == 10:
        #    break

        if len(annotation['RNAseqFile'][0]) > 1:
            RNAFile = annotation['RNAseqFile'][0]
            print(RNAFile)
        else:
            RNAFile = annotation['RNAseqFile']
            print(RNAFile)

        expAccession = RNAFile[-15:-4]
        expFile = annotationFolder +  'geneExp_dict_' + expAccession + '.pkl'
        
        with open(expFile, 'rb') as pickledFile:
            expression = pickle.load(pickledFile) # load the expression file

        # 0.3 Get the promoter/gene body coverage of the labels
        # >>>>>>>>>>>>>>>>>>>>>>>>>>
        inputFile = annotationFolder + 'exp_promoter_labelCover_promLength3000_v02.pkl'
        with open(inputFile, 'rb') as f:
            segMats = pickle.load(f)  # promoter, genes

        inputFile = annotationFolder + 'exp_promoter_labelCover_chmm_promLength3000_v02.pkl'
        with open(inputFile, 'rb') as f:
            chromMats = pickle.load(f)  # promoter, genes

        chromGenes = chromMats['genes']
        segGenes = segMats['genes']

        chromIns = (np.sum(chromGenes, axis=1) >0).astype(int)
        segIns = (np.sum(segGenes, axis=1) > 0).astype(int)

        chromPromo = chromMats['promoter']
        segPromo = segMats['promoter']
        
        chromInsP = (np.sum(chromPromo, axis=1) >0).astype(int)
        segInsP = (np.sum(segPromo, axis=1) > 0).astype(int)

        sib = chromIns + segIns + chromInsP + segInsP
        print(sum(sib==4))
    
        geneFilterAll[accession] = sib == 4

########################################
# 7.1 Get the most similar and different samples regarding gene expression
########################################

rnaAccessionList = list(aucs.keys())
expSim = np.zeros((len(rnaAccessionList), len(rnaAccessionList)))
expCorr = np.zeros((len(rnaAccessionList), len(rnaAccessionList)))
for a, accessionA in enumerate(rnaAccessionList):

    print('a', a)
    
    annotation = allMeta[accessionA]

    geneFilterA = geneFilterAll[accessionA]
        
    annotationFolder = annotation['folder']

    # 0.2 get the expression data
    # >>>>>>>>>>>>>>>>>>>>>>>>>>
    if len(annotation['RNAseqFile'][0]) > 1:
        RNAFile = annotation['RNAseqFile'][0]
    else:
        RNAFile = annotation['RNAseqFile']
    #print(RNAFile)

    expAccession = RNAFile[-15:-4]
    expFile = annotationFolder +  'geneExp_dict_' + expAccession + '.pkl'

    with open(expFile, 'rb') as pickledFile:
        expressionA = pickle.load(pickledFile) # load the expression file

    expArrayA = np.asarray([expressionA[x] for x in geneIDList]) # genes in the transcriptomic data that are in the gene map file

    for b in range(a, len(rnaAccessionList)):
        print(b)
        accessionB = rnaAccessionList[b]
        geneFilterB = geneFilterAll[accessionB]
        annotation = allMeta[accessionB]
        annotationFolder = annotation['folder']
        
        if len(annotation['RNAseqFile'][0]) > 1:
            RNAFile = annotation['RNAseqFile'][0]
        else:
            RNAFile = annotation['RNAseqFile']
            
        #print(RNAFile)

        expAccession = RNAFile[-15:-4]
        expFile = annotationFolder +  'geneExp_dict_' + expAccession + '.pkl'

        with open(expFile, 'rb') as pickledFile:
            expressionB = pickle.load(pickledFile) # load the expression file

        expArrayB = np.asarray([expressionB[x] for x in geneIDList]) # genes in the transcriptomic data that are in the gene map file

        both = np.logical_and(geneFilterA, geneFilterB)

        ea = expArrayA[both]
        eb = expArrayB[both]

        bothExp = np.logical_and(ea>0, eb>0)
        expSim[a, b] = np.sum(bothExp)
        expSim[b, a] = np.sum(bothExp)

        #expCorr[a,b] = np.corrcoef(ea, eb)[0,1]
        #expCorr[b,a] = expCorr[a,b]

sns.heatmap(expCorr)
plt.show()

myData = {'matrix': expCorr, 'accessions':rnaAccessionList}
file = dataFolder + 'expressionCorrelationForRNAsamples.pkl'
with open(file, 'wb') as f:
    pickle.dump(myData, f)

myData = {'matrix': expSim, 'accessions':rnaAccessionList}
file = dataFolder + 'expressionOverlapForRNAsamples.pkl'
with open(file, 'wb') as f:
    pickle.dump(myData, f)

    
sns.heatmap(expSim)
plt.show()


########################################
# 7.2 Get the whole matrix
########################################
geneCount = 24700
rnaAccessionList = list(aucs.keys())
expMat = np.zeros((geneCount, len(rnaAccessionList)))
#expCorr = np.zeros((len(rnaAccessionList), len(rnaAccessionList)))
for a, accessionA in enumerate(rnaAccessionList):

    print('a', a)
    
    annotation = allMeta[accessionA]

    geneFilterA = geneFilterAll[accessionA]
        
    annotationFolder = annotation['folder']

    # 0.2 get the expression data
    # >>>>>>>>>>>>>>>>>>>>>>>>>>
    if len(annotation['RNAseqFile'][0]) > 1:
        RNAFile = annotation['RNAseqFile'][0]
    else:
        RNAFile = annotation['RNAseqFile']
    #print(RNAFile)

    expAccession = RNAFile[-15:-4]
    expFile = annotationFolder +  'geneExp_dict_' + expAccession + '.pkl'

    with open(expFile, 'rb') as pickledFile:
        expressionA = pickle.load(pickledFile) # load the expression file

    expArray = np.asarray([expressionA[x] for x in geneIDList[0:24700]]) # genes in the transcriptomic data that are in the gene map file

    expMat[:,a] = expArray

book = np.copy(expMat)
book[book>0] = 1
sns.heatmap(book)
sns.heatmap(expMat)
plt.show()

expSampleCount = np.sum(book, axis=1)

plt.plot(np.sort(kado))
plt.show()

genesExp = np.where(kado == 10)[0] # in the whole list of genes, those expressed in 10 cells

myMeds = np.zeros(len(genesExp))
for i,s in enumerate(genesExp):
    myExps = expMat[s, :]
    myMeds[i] = np.median(myExps[myExps>0])

print(myMeds[np.where(myMeds > 0.4)])
# 19 is the gene in the list of genes that are expressed in 10 samples 
geneInd = (genesExp[19])
geneExpVector = expMat[geneInd, :]
spp2AccessionList = []
for i in range(len(geneExpVector)):
    if geneExpVector[i] > 1:
        spp2AccessionList.append(rnaAccessionList[i])
        print(geneExpVector[i])
        print(rnaAccessionList[i])
        print(allMeta[rnaAccessionList[i]])

        # get the genes that are expressed in liver,
        # get your liver specific disease
        # get your SNP stuff
['ENCSR789ICN', 'ENCSR900FCU', 'ENCSR214AMQ', 'ENCSR721USS', 'ENCSR891ENU']
4163
print(geneIDList[geneInd])
print(geneList[geneIDList[geneInd]])
print(expMat[geneInd,:])
        
########################################
# 5.1 Segway logistic regression
########################################

# I think I want to use a switch structure.
# For all samples, distribution plots for all the analyses. 1. Unfiltered, Filter gene (lowest 15%, lowest 10%), Expression cut (lowest 10%, lowest 15%), balanced, unbalanced, promoter region, gene body. 5 * 2 * 2 =  20. Box plots.
# The gene body/promoter region I will do in one after another, no switch. Each of them will have five switch for the expression level things. The balanced/unbalanced will also have 2 switch cases. In the loop below, the first few sections are initials. 

# generate list of tupples for the keys. (['gene','promoter'], [1,2,3,4,5], ['balanced', 'unbalanced'])
# regionMode, expMode, regMode
switchLists = [['genes','promoter'], [1,2,3,4,5], ['balanced', 'unbalanced']]

modeKeyList = []
for i in range(2):
    regionMode = switchLists[0][i]
    for j in range(5):
        expMode = switchLists[1][j]
        for k in range(2):
            regMode = switchLists[2][k]
            myKey = (regionMode, expMode, regMode)
            modeKeyList.append(myKey)

# I need to have a list for each of the keys above, that is 20 lists. Each list holds the tupple of the regression outputs. The matrix, the string and the auc.

# TODO: make a loop for the keys

regressionOutput = {}
for i in range(len(modeKeyList)):
    regressionOutput[modeKeyList[i]] = []

regressionModel = {}
for i in range(len(modeKeyList)):
    regressionModel[modeKeyList[i]] = []

counter = 0
for accession in accessionList:
    
    if accession in list(aucs.keys()):
        annotation = allMeta[accession]

        geneFilter = geneFilterAll[accession]
        
        annotationFolder = annotation['folder']

        # 0.1 get the mnemonics
        # >>>>>>>>>>>>>>>>>>>>>>>>>>
        label_term_mapping = {}
        mnemonics_file = annotationFolder + 'mnemonics_v04.txt'
        with open(mnemonics_file, 'r') as mnemonics:
            for line in mnemonics:
                print(line)
                label = line.strip().split()[0]
                term = line.strip().split()[1]
                label_term_mapping[label] = term

        # 0.2 get the expression data
        # >>>>>>>>>>>>>>>>>>>>>>>>>>
        if len(annotation['RNAseqFile'][0]) > 1:
            RNAFile = annotation['RNAseqFile'][0]
        else:
            RNAFile = annotation['RNAseqFile']
        print(RNAFile)

        expAccession = RNAFile[-15:-4]
        expFile = annotationFolder +  'geneExp_dict_' + expAccession + '.pkl'

        with open(expFile, 'rb') as pickledFile:
            expression = pickle.load(pickledFile) # load the expression file

        # 0.3 Get the promoter/gene body coverage of the labels
        # >>>>>>>>>>>>>>>>>>>>>>>>>>
        inputFile = annotationFolder + 'exp_promoter_labelCover_promLength3000_v02.pkl'
        with open(inputFile, 'rb') as f:
            transPromoMat = pickle.load(f)  # promoter, genes

        # TODO: add the KEY
        for key in modeKeyList:
            regionMode = key[0]
            expMode = key[1]
            regMode = key[2]

            # 1. For Genes OR promoters
            # >>>>>>>>>>>>>>>>>>>>>>>>>>

            if regionMode == 'genes':
                labelExpMat = transPromoMat['genes']
            else:
                labelExpMat = transPromoMat['promoter']
        
            #sib = labelExpMat.sum(axis = 1)
            #filterGene = sib > 0 # genes that are present in the annotation

            expArray = np.asarray([expression[x] for x in geneIDList]) # genes in the transcriptomic data that are in the gene map file
            filterExp = expArray[geneFilter,] # from these, we filter for those that have annotation

            # 2. The expression filter switch mechanism, it has 5 switches 
            # >>>>>>>>>>>>>>>>>>>>>>>>>>

            # the outcome of each switch is the "notExp" and "filterExpMat"
            # the outcome of this section is "notExp", and the filtered (based on switches) and processed "expMat"
        
            # >>>>> 2.1 Switch 1
            if expMode == 1: # default
                notExp = filterExp == 0
                
                filterExpMat = labelExpMat[geneFilter, :] # now filtering the labelExpMat for genes that have annotation
            # >>>>> 2.2 Switch 2
            expValues = filterExp[filterExp > 0]
            if expMode == 2: # 10% low expressed filtered
                # get the 10% quantile of the filterExp
                q10 = np.quantile(expValues, .1)
                filterForLowExp = ~(((filterExp > 0).astype(int) + (filterExp < q10).astype(int)) == 2)
                filFilExp = filterExp[filterForLowExp,] ## >>>>> removing those with low exp
                filterExp = filFilExp
            
                notExp = filterExp == 0
            
                filterExpMat = labelExpMat[geneFilter, :] # now filtering the labelExpMat for genes that have annotation
                filFilExpMat = filterExpMat[filterForLowExp, :]## >>>>> filter for low Exp
                filterExpMat = filFilExpMat

            # >>>>> 2.3 Switch 3
            if expMode == 3: # 15% low expressed filtered
                q15 = np.quantile(expValues, .15)
                filterForLowExp = ~(((filterExp > 0).astype(int) + (filterExp < q15).astype(int)) == 2)
                filFilExp = filterExp[filterForLowExp,] ## >>>>> removing those with low exp
                filterExp = filFilExp
            
                notExp = filterExp == 0
            
                filterExpMat = labelExpMat[geneFilter, :] # now filtering the labelExpMat for genes that have annotation
                filFilExpMat = filterExpMat[filterForLowExp, :]## >>>>> filter for low Exp
                filterExpMat = filFilExpMat

            # >>>>> 2.4 Switch 4
            if expMode == 4: # 10% low expressed considered zero expressed
                q10 = np.quantile(expValues, .1)
                notExp = filterExp < q10
            
                filterExpMat = labelExpMat[geneFilter, :] # now filtering the labelExpMat for genes that have annotation

            # >>>>> 2.5 Switch 5
            if expMode == 5: # 15% low expressed considered zero expressed
                q15 = np.quantile(expValues, .15)
                notExp = filterExp < q15
            
                filterExpMat = labelExpMat[geneFilter, :] # now filtering the labelExpMat for genes that have annotation

            sumExp = filterExpMat.sum(axis = 1)
            filterExpMat = filterExpMat/sumExp[:,None] # normalizing the coverage to get ratio of the coverage by the labe
        
            expMat = np.zeros((filterExpMat.shape[0], segwayStateCount))

            # get the feature mat
            for i,label in enumerate(label_term_mapping):
                index = segwayStates.index(label_term_mapping[label])
                expMat[:, index] += filterExpMat[:, i]

            # split the feature mat
            X_train, X_test, y_train, y_test = train_test_split(expMat, notExp, test_size=0.25, random_state=11)

            # 3. The regression balance/unbalanced swithc mechanism, it has 2
            # >>>>>>>>>>>>>>>>>>>>>>>>>>

            if regMode == 'balanced': # if we want balanced regression
                logreg = LogisticRegression(random_state=11, class_weight='balanced')

            else: # if we want none balanced
                logreg = LogisticRegression(random_state=11)


            # 4. the outcome of the regression
            # >>>>>>>>>>>>>>>>>>>>>>>>>>

            logreg.fit(X_train, y_train)

            y_pred = logreg.predict(X_test)
            y_prob = logreg.predict_proba(X_test)
            cnf_matrix = metrics.confusion_matrix(y_test, y_pred) # a matrix

            print(cnf_matrix)
        
            target_names = ['expressed', 'not expressed']
            #print(classification_report(y_test, y_pred, target_names=target_names)) # an str
            report =(classification_report(y_test, y_pred, target_names=target_names)) 

            auc = accuracy_score(y_test, y_pred) # a float

            regressionOutput[key].append((accession, auc, cnf_matrix, report))
            regressionModel[key].append((accession, logreg))
            print(counter)
            print(key)
            print(accession)
            counter +=1

### save the data

file = dataFolder + 'Segway_logisticRegression_output_v04_expFile_v02.pkl'
with open(file, 'wb') as f:
    pickle.dump(regressionOutput, f)

with open(file, 'rb') as f:
    regOutput = pickle.load(f)

### save the model
file = dataFolder + 'Segway_logisticRegression_model_v04_expFile_v02.pkl'
with open(file, 'wb') as f:
    pickle.dump(regressionModel, f)
    
file = dataFolder + 'Segway_logisticRegression_model_v04.pkl'
with open(file, 'rb') as f:
    regModel = pickle.load(f)

# 5.1. * the ROC plot
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

y_pred_proba = logreg.predict_proba(X_test)[::,1]
fpr, tpr, _ = metrics.roc_curve(y_test,  y_pred_proba)

#create ROC curve
plt.plot(fpr,tpr)
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.show()


########################################
# 5.2 ChromHMM logistic regression (almost exactly the same)
########################################

# regionMode, expMode, regMode
switchLists = [['genes','promoter'], [1,2,3,4,5], ['balanced', 'unbalanced']]

modeKeyList = []
for i in range(2):
    regionMode = switchLists[0][i]
    for j in range(5):
        expMode = switchLists[1][j]
        for k in range(2):
            regMode = switchLists[2][k]
            myKey = (regionMode, expMode, regMode)
            modeKeyList.append(myKey)

regressionOutput = {}
for i in range(len(modeKeyList)):
    regressionOutput[modeKeyList[i]] = []

regressionModel = {}
for i in range(len(modeKeyList)):
    regressionModel[modeKeyList[i]] = []

counter = 0
for accession in accessionList:
    
    if accession in list(aucs.keys()):
        annotation = allMeta[accession]

        geneFilter = geneFilterAll[accession]
        
        annotationFolder = annotation['folder']

        # 0.2 get the expression data
        # >>>>>>>>>>>>>>>>>>>>>>>>>>
        if len(annotation['RNAseqFile'][0]) > 1:
            RNAFile = annotation['RNAseqFile'][0]
        else:
            RNAFile = annotation['RNAseqFile']
        print(RNAFile)

        expAccession = RNAFile[-15:-4]
        expFile = annotationFolder +  'geneExp_dict_' + expAccession + '.pkl'

        with open(expFile, 'rb') as pickledFile:
            expression = pickle.load(pickledFile) # load the expression file

        # 0.3 Get the promoter/gene body coverage of the labels
        # >>>>>>>>>>>>>>>>>>>>>>>>>>
        inputFile = annotationFolder + 'exp_promoter_labelCover_chmm_promLength3000_v02.pkl'
        with open(inputFile, 'rb') as f:
            transPromoMat = pickle.load(f)  # promoter, genes

        for key in modeKeyList:
            regionMode = key[0]
            expMode = key[1]
            regMode = key[2]

            # 1. For Genes OR promoters
            # >>>>>>>>>>>>>>>>>>>>>>>>>>

            if regionMode == 'genes':
                labelExpMat = transPromoMat['genes']
            else:
                labelExpMat = transPromoMat['promoter']
            #sib = labelExpMat.sum(axis = 1)
            #filterGene = sib > 0 # genes that are present in the annotation

            expArray = np.asarray([expression[x] for x in geneIDList]) # genes in the transcriptomic data that are in the gene map file
            filterExp = expArray[geneFilter,] # from these, we filter for those that have annotation

            # 2. The expression filter switch mechanism, it has 5 switches 
            # >>>>>>>>>>>>>>>>>>>>>>>>>>

            # the outcome of each switch is the "notExp" and "filterExpMat"
            # the outcome of this section is "notExp", and the filtered (based on switches) and processed "expMat"
        
            # >>>>> 2.1 Switch 1
            if expMode == 1: # default
                notExp = filterExp == 0
                
                filterExpMat = labelExpMat[geneFilter, :] # now filtering the labelExpMat for genes that have annotation
            # >>>>> 2.2 Switch 2
            expValues = filterExp[filterExp > 0]
            if expMode == 2: # 10% low expressed filtered
                # get the 10% quantile of the filterExp
                q10 = np.quantile(expValues, .1)
                filterForLowExp = ~(((filterExp > 0).astype(int) + (filterExp < q10).astype(int)) == 2)
                filFilExp = filterExp[filterForLowExp,] ## >>>>> removing those with low exp
                filterExp = filFilExp
            
                notExp = filterExp == 0
            
                filterExpMat = labelExpMat[geneFilter, :] # now filtering the labelExpMat for genes that have annotation
                filFilExpMat = filterExpMat[filterForLowExp, :]## >>>>> filter for low Exp
                filterExpMat = filFilExpMat

            # >>>>> 2.3 Switch 3
            if expMode == 3: # 15% low expressed filtered
                q15 = np.quantile(expValues, .15)
                filterForLowExp = ~(((filterExp > 0).astype(int) + (filterExp < q15).astype(int)) == 2)
                filFilExp = filterExp[filterForLowExp,] ## >>>>> removing those with low exp
                filterExp = filFilExp
            
                notExp = filterExp == 0
            
                filterExpMat = labelExpMat[geneFilter, :] # now filtering the labelExpMat for genes that have annotation
                filFilExpMat = filterExpMat[filterForLowExp, :]## >>>>> filter for low Exp
                filterExpMat = filFilExpMat

            # >>>>> 2.4 Switch 4
            if expMode == 4: # 10% low expressed considered zero expressed
                q10 = np.quantile(expValues, .1)
                notExp = filterExp < q10
            
                filterExpMat = labelExpMat[geneFilter, :] # now filtering the labelExpMat for genes that have annotation

            # >>>>> 2.5 Switch 5
            if expMode == 5: # 15% low expressed considered zero expressed
                q15 = np.quantile(expValues, .15)
                notExp = filterExp < q15
            
                filterExpMat = labelExpMat[geneFilter, :] # now filtering the labelExpMat for genes that have annotation

            sumExp = filterExpMat.sum(axis = 1)
            filterExpMat = filterExpMat/sumExp[:,None] # normalizing the coverage to get ratio of the coverage by the labe
        
            expMat = filterExpMat

            # split the feature mat
            X_train, X_test, y_train, y_test = train_test_split(expMat, notExp, test_size=0.25, random_state=11)

            # 3. The regression balance/unbalanced swithc mechanism, it has 2
            # >>>>>>>>>>>>>>>>>>>>>>>>>>

            if regMode == 'balanced': # if we want balanced regression
                logreg = LogisticRegression(random_state=11, class_weight='balanced')

            else: # if we want none balanced
                logreg = LogisticRegression(random_state=11)


            # 4. the outcome of the regression
            # >>>>>>>>>>>>>>>>>>>>>>>>>>

            logreg.fit(X_train, y_train)

            y_pred = logreg.predict(X_test)
            y_prob = logreg.predict_proba(X_test)
            cnf_matrix = metrics.confusion_matrix(y_test, y_pred) # a matrix

            print(cnf_matrix)
        
            target_names = ['expressed', 'not expressed']
            #print(classification_report(y_test, y_pred, target_names=target_names)) # an str
            report = (classification_report(y_test, y_pred, target_names=target_names)) 

            auc = accuracy_score(y_test, y_pred) # a float

            regressionOutput[key].append((accession, auc, cnf_matrix, report))
            regressionModel[key].append((accession, logreg))
            print(counter)
            print(key)
            print(accession)
            counter +=1

### save the data

file = dataFolder + 'Chrom_logisticRegression_output_v04_expFile_v02.pkl'
with open(file, 'wb') as f:
    pickle.dump(regressionOutput, f)

with open(file, 'rb') as f:
    regOutput = pickle.load(f)

### save the model

file = dataFolder + 'Chrom_logisticRegression_model_v04_expFile_v02.pkl'
with open(file, 'wb') as f:
    pickle.dump(regressionModel, f)

with open(file, 'rb') as f:
    regModel = pickle.load(f)


########################################
# 5.3 Box plots for the LogR  
########################################
plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/'
# Segway: I am goint to do 4 panel plots: promoter, gene, balanced, unbalanced. Each panel will have the 5 box plots. There are list of 20 keys.
switchLists = [['genes','promoter'], [1,2,3,4,5], ['balanced', 'unbalanced']]

file = dataFolder + 'Segway_logisticRegression_output_v04_expFile_v02.pkl'
file = dataFolder + 'Chrom_logisticRegression_output_v04_expFile_v02.pkl'
with open(file, 'rb') as f:
    regOutput = pickle.load(f)

modeKeyList = []
aucMats = []
titles = []
for i in range(2):
    regionMode = switchLists[0][i]
    for k in range(2):
        regMode = switchLists[2][k]
        aucMat = np.zeros((88, 5))
        for j in range(5):
            expMode = switchLists[1][j]
            myKey = (regionMode, expMode, regMode)
            output = regOutput[myKey]
            aucMat[:,j] = [kado[1] for kado in output]
            print(myKey)

        aucMats.append(aucMat)
        titles.append(myKey)
        
fig, axs = plt.subplots(1,4, figsize=(14,4))
bplist = []
for i in range(4):
    print(i)
    bp=axs[i].boxplot(aucMats[i], patch_artist=True)
    axs[i].set_ylim([.5,.95])
    title = titles[i][0] + ' ' + titles[i][2]
    axs[i].set_title(title)
    bplist.append(bp)
    
colors = ['pink', 'lightblue', 'lightgreen', 'goldenrod', 'royalblue']
for bplot in bplist:
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)

legends = ['default: exp == 0, exp > 0', 'exp == 0, exp > 10 centile exp', 'exp == 0, exp > 15 centile exp',
           'exp =< 10 centile exp, exp > 10 centile exp', 'exp <= 15 centile exp, exp > 15 centile exp']

from matplotlib.patches import Patch

legend_elements = [Patch(facecolor='pink', label = legends[0]),
                   Patch(facecolor='lightblue', label = legends[1]),
                   Patch(facecolor='lightgreen', label = legends[2]),
                   Patch(facecolor='goldenrod', label = legends[3]),
                   Patch(facecolor='royalblue', label = legends[4])]

axs[3].legend(handles=legend_elements)
#leg = axs[3].get_legend()
#leg.legendHandles[0].set_color('pink')
#leg.legendHandles[1].set_color('lightblue')

plt.show()
figFile = plotFolder + 'Segway_logisticRegression_self_v02.pdf'
figFile = plotFolder + 'Chrom_logisticRegression_self_v02.pdf'
print(figFile)
plt.savefig(figFile, bbox_inches='tight')
plt.close('all')


