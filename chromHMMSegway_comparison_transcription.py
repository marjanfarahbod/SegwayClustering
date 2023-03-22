# how do they compare regarding the transcription data.
# How well do they do regarding the transcription data? active prmoter
# How well do they do regarding the transcription data? transcribed region etc.

# expressed genes act very much the same between the low and high expression.

# for each gene, get the ratio of its promoter region and its region covered with each of the labels
# examine the distribution of these labels between the expressed and non expressed genes.

# do this for both chromhmm and segway. hahaha.

# 0. Initials 

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

segwayLabels = ['Quiescent', 'ConstitutiveHet', 'FacultativeHet', 'Transcribed', 'Promoter', 'Enhancer', 'RegPermissive', 'Bivalent', 'LowConfidence']

# How well do they do regarding the transcription data? active prmoter
# How well do they do regarding the transcription data? transcribed region etc.

# expressed genes act very much the same between the low and high expression.

# for each gene, get the ratio of its promoter region and its region covered with each of the labels
# examine the distribution of these labels between the expressed and non expressed genes.

# do this for both chromhmm and segway. hahaha.
########################################

extension = 3000 # count of basepairs monitored before and after the gene coordinates

annAccessionCount = 0
for annAccession in annAccessionList:
    
    # for each annotation
    if (annMeta[annAccession]['expressionFile'] != 'none'):
        
        print(annAccession)
        annAccessionCount += 1
        print(annAccessionCount)
        
        expFileName = annMeta[annAccession]['expressionFile']
        expAccession = re.split('_|\.', expFileName)[2]
        inputFile = dataFolder + dataSubFolder + annAccession + '/geneExp_dict_' + expAccession + '.pkl'
        print(inputFile)
        with open(inputFile, 'rb') as pickledFile:
            expression = pickle.load(pickledFile)

        # using the sorted and filtered .bed file instead of the initial one
        annFile = dataFolder + dataSubFolder + annAccession + '/' + annMeta[annAccession]['bedFile'].split('.')[0] + '_filteredSorted.bed.gz'

        # TODO: unzip the annFile
        os.system('gunzip %s' %(annFile))

        annFile = annFile[0:-3]

        # loading the annotation summary file
        bedFileName = annMeta[annAccession]['bedFile'].split('.')[0]
        inputFile = dataFolder + dataSubFolder + annAccession + '/' + bedFileName + '_annotationSummary.pkl'
        with open(inputFile, 'rb') as pickledFile:
            summaryAnnotation = pickle.load(pickledFile)

        # count of lines for annotation file
        annLineCount = int(sp.getoutput('wc -l %s' %(annFile)).split()[0])

        clusterCount = len(summaryAnnotation['clusters'])
        sortedClusterList = list(range(0, clusterCount))

        print(sortedClusterList)
            
        cgi = 0 # walks on the geneIDList
        ann_start = 0 # the start of the current annotation
        ann_end = 0 # the end of the current annotation
        ann_chr = 'chr'
        ann_line_count = 0 # this is just to check the progress through the annotation file
        previous_class = ''

        labelExpMat = np.zeros((len(geneIDList), clusterCount))
        labelPromoterMat = np.zeros((len(geneIDList), clusterCount))
        
        annLineInd = 1 # this keeps the annotation line index, the first line is empty, thus the index starts at 1
        genomic_region_start_annLineInd = 0

        previous_gene_chr = 'chr1'
        previous_extension_end = 0 # doesn't matter since chr condition is never true until the first gene is processed
        previous_ann_chr = 'chr1'

        #while cgi < len(geneIDList) and annLineInd < annLineCount: # >>>> test: modify the condition for the test runs
        while cgi < 26017: # >>>>>>>>>> MAIN
        #sib = 1
        #while sib > 0 and cgi <10:
        #while cgi < 18431:# >>>>
            
            
            print('cgi')
            print(cgi)

            'we are in the gene territory, we are walking on the genes'
   
            geneID = geneIDList[cgi]
            gene_chr = geneList[geneIDList[cgi]].chrom
            #print(geneID) # >>>> test
            
            strand = geneList[geneID].strand
            gene_start = geneList[geneID].start
            gene_end = geneList[geneID].end

            #print('geneStuff')
            #print(gene_start)
            #print(gene_end)

            gene_length = gene_end - gene_start 
            promoter_length = 1400
            gene_coverage = 0
            promoter_coverage = 0
            #gene_length = gene_end - gene_start
            #gene_length_unit = int(gene_length/100)
            #if gene_length_unit == 0:
            #    #print('gene smaller than 100bp, skipping it') # >>>> test
            #    cgi = cgi + 1
            #    continue

            #gene_length_last_unit = gene_length - (99* gene_length_unit)
            # TODO: ideally: something to fix: the gene_length_last_unit for negative strand versus positive strand

            if strand == '+':
                gene_start = gene_start - promoter_length

            if strand == '-':
                gene_end = gene_end + promoter_length
            
            #extension_start = gene_start - extension
            #extension_end = gene_end + extension

            #print('%ss, %s, %s' %(gene_chr, extension_start, extension_end)) # >>>> test

            
            #''' picking the label/class matrix based on the gene expression level'''
            #TODO catch exception for when geneID is not in expression
            #gene_exp = expression[geneID]
            #if gene_exp == 0:
            #    expMatInd = 0
            #elif gene_exp > 1:
            #    expMatInd = 2
            #else:
            #    expMatInd = 1

            #print('expMatInd')
            #print(expMatInd) # >>>>> test

            # reading the next annotation
            line = linecache.getline(annFile, annLineInd)
            annLineInd +=1
            ann_line_count += 1
            fields = line.strip().split()

            previous_ann_chr = ann_chr
            ann_chr = fields[0]
            ann_start = int(fields[1])
            ann_end = int(fields[2])
            #print('first ann')
            #print(fields)

            if (gene_chr != previous_gene_chr): # in case of chromosome change because gene moved to the next chromosome
                print('gene chr not equal to previous gene chr')
                print('gene_chr %s' %(gene_chr))
                print('p_gene_chr %s' %(previous_gene_chr))
                print('ann_chr %s' %(ann_chr))
                
                while (ann_chr != gene_chr): # move on in the annotation until we reach to the next chromosome
                    #print('reading annotations until it is equal') # >>>> test
                    line = linecache.getline(annFile, annLineInd)
                    fields = line.strip().split()
                    annLineInd +=1
                    ann_line_count += 1

                    previous_ann_chr = ann_chr                
                    ann_chr = fields[0]
                    ann_start = int(fields[1])
                    ann_end = int(fields[2])

                print(ann_chr)
                print(previous_ann_chr)

                
            if (ann_chr != gene_chr): # if annotation moved to the next chromosome, but gene has not yet moved to the next chromosome
                annLineInd = annLineInd - 2
                line = linecache.getline(annFile, annLineInd)
                annLineInd +=1
                ann_line_count += 1
                fields = line.strip().split()

                ann_chr = fields[0]
                ann_start = int(fields[1])
                ann_end = int(fields[2])
                print('marker 01')


            ''' 
            if this gene starts somewhere in the preivous genomic region, 
            I will go back in the annotation file to the beginning of annotation for the previous gene
            changed this th the next while loop - basically we will go back until the annotations start before the gene territory

            '''

            while (ann_start > gene_start) and (gene_chr == ann_chr): # in case of overlapping genes
                #print('ann start greater than extension start, getting back in annotation until it is not') # >>>> test
                #print(annLineInd) # >>>> test
                #print('%ss, %s, %s' %(ann_chr, ann_start, ann_end)) # >>>> test
                
                annLineInd = max(annLineInd - 5,1)
                line = linecache.getline(annFile, annLineInd)
                annLineInd +=1
                ann_line_count += 1
                fields = line.strip().split()

                ann_chr = fields[0]
                ann_start = int(fields[1])
                ann_end = int(fields[2])
                #print('overlapping genes')
                
                #print('overlapping genes here') # >>>> test
                #print(annLineInd) # >>>> test

            #while ((ann_start < extension_end) or not(gene_chr == ann_chr)) and geneMatWalkIndex < 160: 
            while ((ann_start < gene_end) and (gene_chr == ann_chr)) and gene_coverage+promoter_coverage < gene_length+promoter_length: #and geneMatWalkIndex < 160: # while we still have annotation before the gene

                #print('main while')

                # in the next sections we are processing the annotation
                if ann_chr == gene_chr: # if we are in the same choromosome # I think we are, that is the condition
                    if not((ann_end < gene_start) or ann_start > gene_end):
                        #print('genomic region')
                    #if ((ann_start < gene_end and ann_start > extension_start) or 
                    #   (ann_end < extension_end and ann_end > extension_start) or
                    #   (ann_start < extension_start and ann_end > extension_end)):

                        ''' We are in the genonimc region (with extension)'''
                        #print('annotation') # >>>> test
                        #print('%ss, %s, %s' %(ann_chr, ann_start, ann_end)) # >>>> test

                        ''' Taking off for filling the matrices... '''
                        
                        adjusted_ann_start = max(0, ann_start - gene_start)
                        adjusted_ann_end = min(ann_end - gene_start, gene_end - gene_start)
                        adjusted_ann_length = adjusted_ann_end - adjusted_ann_start

                        ann_cluster = fields[3].split('_')[0]
                        clusterInd = int(ann_cluster)

                        if gene_strand == '+':  # munching 100bp s from annotation, filling the geneMat
                            #print('positive strand')

                            if promoter_coverage < promoter_length:
                                new_coverage = min(adjusted_ann_length, promoter_length - promoter_coverage)
                                labelPromoterMat[cgi, clusterInd] += new_coverage
                                adjusted_ann_length = adjusted_ann_length - new_coverage
                                promoter_coverage += new_coverage
                            if gene_coverage < gene_length and adjusted_ann_length > 0:
                                new_coverage = min(adjusted_ann_length, gene_length - gene_coverage)
                                labelExpMat[cgi, clusterInd] += new_coverage
                                adjusted_ann_length = adjusted_ann_length - new_coverage
                                gene_coverage += new_coverage
                            #coverPromoter
 
                            #coverGene
                        if gene_strand == '-':
                            #print('negative strand')
                            if gene_coverage < gene_length:
                                new_coverage = min(adjusted_ann_length, gene_length - gene_coverage)
                                labelExpMat[cgi, clusterInd] += new_coverage
                                adjusted_ann_length = adjusted_ann_length - new_coverage
                            if promoter_coverage < promoter_length and adjusted_ann_length >0:
                                new_coverage = min(adjusted_ann_length, promoter_length - promoter_coverage)
                                labelPromoterMat[cgi, clusterInd] += new_coverage
                                adjusted_ann_length = adjusted_ann_length - new_coverage


                
                if gene_coverage+promoter_coverage < gene_length+promoter_length:
                    #print('not finished with gene coverage yet')
                    #print('geneMatWalkInd') # >>>> test
                    #print(geneMatWalkIndex) # >>>> test
                    #print('annLineInd') # >>>> test
                    #print(annLineInd) # >>>> test
                    line = linecache.getline(annFile, annLineInd)
                    annLineInd +=1
                    #print(annLineInd)
                    ann_line_count += 1
                    fields = line.strip().split()

                    ann_chr = fields[0]
                    ann_start = int(fields[1])
                    ann_end = int(fields[2])
                    #print('just going fw')
                    #print(line)
                    #print('%ss, %s, %s' %(ann_chr, ann_start, ann_end)) # >>>> test

            sib = labelExpMat[cgi,:].sum()
            print(sib)
            cgi += 1 # next gene
            #print(cgi)
            previous_gene_end = gene_end
            previous_gene_chr = gene_chr
            

        linecache.clearcache()

        transPromoMat = {"genes": labelExpMat, "promoter": labelPromoterMat}
        outputFile = dataFolder + dataSubFolder + annAccession + '/' + 'exp_promoter_labelCover.pkl'
        with open(outputFile, 'wb') as f:
            pickle.dump(transPromoMat, f)

        labelExpMat = labelExpMat/labelExpMat.sum(axis=1)[:,None]

        
inputFile = dataFolder + dataSubFolder + annAccession + '/' + 'exp_promoter_labelCover.pkl'
with open(inputFile, 'rb') as f:
    transPromoMat = pickle.load( f)

labelExpMat = transPromoMat['promoter']
sib = labelExpMat.sum(axis = 1)
sib[sib == 0] = 1

labelExpMat = labelExpMat/sib[:,None]

####

    # load the mnomics, and extract the enhancer labels for mnemonics
    label_term_mapping = {}
    mnemonics_file = dataFolder + dataSubFolder + annAccession + '/' + 'mnemonics_v04.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term

    expArray = np.asarray([expression[x] for x in geneIDList])
    
    plt.hist(expArray)
    plt.show()

    notExp = expArray == 0
    exp = expArray > 0

    fig = plt.figure(figsize =(17, 7))
    ax = fig.add_subplot(111)
 
    data = []
    for i in range(labelCount):
        expVals = labelExpMat[exp, i]
        notExpVals = labelExpMat[notExp, i]
        data.append(expVals)
        data.append(notExpVals)
        
        plt.boxplot(data)
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)

        plt.show()
        
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)


    book = lableExpMat[]


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



