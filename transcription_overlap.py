# The two functions for chromhmm and Segway for transcription overlap.
# 1. Segway overlap of the labels in genomic regions
# 2. ChromHMM overlap of the labels in genomic regions
# 3. Segway enrichment of the labels in genomic regions (from QC_transcriptionComparison_03.py)
# 4. Segway enrichment for labels in genomic regions (no transcript file required)
import pickle
import os
import numpy as np
import pandas as pd
import subprocess as sp
import linecache

########################################
# 1. Segway overlap of the labels in genomic regions
########################################

def EPLabelCover_Segway(annotationFolder, expFile, annFile, mnemFile, geneList, geneIDList, promLength):

    # load expression data
    with open(expFile, 'rb') as f:
        expression = pickle.load(f)

    # prepare ann
    if annFile.endswith('.gz'):
        os.system('gunzip %s' %(annFile))
    else:
        os.system('gunzip %s.gz' %(annFile))

    print(annFile)

    annLineCount = int(sp.getoutput('wc -l %s' %(annFile)).split()[0])

    # get the mnemonics
    label_term_mapping = {}
    with open(mnemFile, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term

    clusterCount = len(label_term_mapping)
    sortedClusterList = list(range(0,clusterCount))

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

    while cgi < 26017: # >>>>>>>>>> MAIN

        'we are in the gene territory, we are walking on the genes'
        
        geneID = geneIDList[cgi]
        gene_chr = geneList[geneIDList[cgi]].chrom
            
        strand = geneList[geneID].strand
        gene_start = geneList[geneID].start
        gene_end = geneList[geneID].end

        gene_length = gene_end - gene_start - 300
        promoter_length = 1500
        gene_coverage = 0
        promoter_coverage = 0

        if strand == '+':
            gene_start = gene_start - promoter_length + 300

        if strand == '-':
            gene_end = gene_end + promoter_length - 300
                

        # reading the next annotation
        line = linecache.getline(annFile, annLineInd)
        annLineInd +=1
        ann_line_count += 1
        fields = line.strip().split()
            
        previous_ann_chr = ann_chr
        ann_chr = fields[0]
        ann_start = int(fields[1])
        ann_end = int(fields[2])

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
                
            annLineInd = max(annLineInd - 5,1)
            line = linecache.getline(annFile, annLineInd)
            annLineInd +=1
            ann_line_count += 1
            fields = line.strip().split()

            ann_chr = fields[0]
            ann_start = int(fields[1])
            ann_end = int(fields[2])

        while ((ann_start < gene_end) and (gene_chr == ann_chr)) and gene_coverage+promoter_coverage < gene_length+promoter_length: #and geneMatWalkIndex < 160: # while we still have annotation before the gene

            # in the next sections we are processing the annotation
            if not((ann_end < gene_start) or ann_start > gene_end):

                ''' We are in the genonimc region (with extension)'''

                ''' Taking off for filling the matrices... '''
                        
                adjusted_ann_start = max(0, ann_start - gene_start)
                adjusted_ann_end = min(ann_end - gene_start, gene_end - gene_start)
                adjusted_ann_length = adjusted_ann_end - adjusted_ann_start

                ann_cluster = fields[3].split('_')[0]
                clusterInd = int(ann_cluster)

                if strand == '+':  # munching 100bp s from annotation, filling the geneMat

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
 

                if strand == '-':

                    if gene_coverage < gene_length:
                        new_coverage = min(adjusted_ann_length, gene_length - gene_coverage)
                        labelExpMat[cgi, clusterInd] += new_coverage
                        adjusted_ann_length = adjusted_ann_length - new_coverage
                        gene_coverage += new_coverage
                    if promoter_coverage < promoter_length and adjusted_ann_length >0:
                        new_coverage = min(adjusted_ann_length, promoter_length - promoter_coverage)
                        labelPromoterMat[cgi, clusterInd] += new_coverage
                        adjusted_ann_length = adjusted_ann_length - new_coverage
                        promoter_coverage += new_coverage


                
            if gene_coverage+promoter_coverage < gene_length+promoter_length:
                line = linecache.getline(annFile, annLineInd)
                annLineInd +=1

                ann_line_count += 1
                fields = line.strip().split()

                ann_chr = fields[0]
                ann_start = int(fields[1])
                ann_end = int(fields[2])

        cgi += 1 # next gene
        previous_gene_end = gene_end
        previous_gene_chr = gene_chr
            

    linecache.clearcache()
    os.system('gzip %s' %(annFile))

    transPromoMat = {"genes": labelExpMat, "promoter": labelPromoterMat}
    outputFile = annotationFolder  + 'exp_promoter_labelCover_promLength%d.pkl' %(promLength)
    with open(outputFile, 'wb') as f:
        pickle.dump(transPromoMat, f)

    print(outputFile)

########################################
# 2. ChromHMM Overlap of the labels in genomic regions
########################################

def EPLabelCover_chromHMM(annotationFolder, expFile, annFile, chromLabels, geneList, geneIDList, promLength):

    # load expression data
    with open(expFile, 'rb') as f:
        expression = pickle.load(f)

    # prepare ann
    if annFile.endswith('.gz'):
        os.system('gunzip %s' %(annFile))
        annFile = annFile[0:-3]
    else:
        os.system('gunzip %s.gz' %(annFile))


    clusterCount = len(chromLabels)
    annLineCount = int(sp.getoutput('wc -l %s' %(annFile)).split()[0])

    sortedClusterList = chromLabels

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
    while cgi < 25800: # >>>>>>>>>> MAIN
        #sib = 1
        #while sib > 0 and cgi <10:
        #while cgi < 18431:# >>>>
        
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
        promoter_length = 1500
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
            
        line = linecache.getline(annFile, annLineInd)
        annLineInd +=1
        ann_line_count += 1
        fields = line.strip().split()

        previous_ann_chr = ann_chr
        ann_chr = fields[0]
        ann_start = int(fields[1])
        ann_end = int(fields[2])

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
            #print('marker 01')


        ''' 
        if this gene starts somewhere in the preivous genomic region, 
        I will go back in the annotation file to the beginning of annotation for the previous gene
        changed this th the next while loop - basically we will go back until the annotations start before the gene territory

        '''
            
        while (ann_start > gene_start) and (gene_chr == ann_chr): # in case of overlapping genes
                
            annLineInd = max(annLineInd - 5,1)
            line = linecache.getline(annFile, annLineInd)
            annLineInd +=1
            ann_line_count += 1
            fields = line.strip().split()

            ann_chr = fields[0]
            ann_start = int(fields[1])
            ann_end = int(fields[2])
                

        while ((ann_start < gene_end) and (gene_chr == ann_chr)) and gene_coverage+promoter_coverage < gene_length+promoter_length: #and geneMatWalkIndex < 160: # while we still have annotation before the gene

            if not((ann_end < gene_start) or ann_start > gene_end):
                
                ''' We are in the genonimc region (with extension)'''

                ''' Taking off for filling the matrices... '''
                        
                adjusted_ann_start = max(0, ann_start - gene_start)
                adjusted_ann_end = min(ann_end - gene_start, gene_end - gene_start)
                adjusted_ann_length = adjusted_ann_end - adjusted_ann_start

                ann_cluster = fields[3]
                clusterInd = chromLabels.index(ann_cluster)

                if strand == '+':  # munching 100bp s from annotation, filling the geneMat
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
 
                if strand == '-':
                    if gene_coverage < gene_length:
                        new_coverage = min(adjusted_ann_length, gene_length - gene_coverage)
                        labelExpMat[cgi, clusterInd] += new_coverage
                        adjusted_ann_length = adjusted_ann_length - new_coverage
                        gene_coverage += new_coverage
                    if promoter_coverage < promoter_length and adjusted_ann_length >0:
                        new_coverage = min(adjusted_ann_length, promoter_length - promoter_coverage)
                        labelPromoterMat[cgi, clusterInd] += new_coverage
                        adjusted_ann_length = adjusted_ann_length - new_coverage
                        promoter_coverage += new_coverage


                
            if gene_coverage+promoter_coverage < gene_length+promoter_length:

                line = linecache.getline(annFile, annLineInd)
                annLineInd +=1

                ann_line_count += 1
                fields = line.strip().split()

                ann_chr = fields[0]
                ann_start = int(fields[1])
                ann_end = int(fields[2])
                

        cgi += 1 # next gene
        #print(cgi)
        previous_gene_end = gene_end
        previous_gene_chr = gene_chr
            

    linecache.clearcache()

    os.system('gzip %s' %(annFile))

    transPromoMat = {"genes": labelExpMat, "promoter": labelPromoterMat}
    outputFile = annotationFolder + 'exp_promoter_labelCover_chmm_promLength%d.pkl' %(promLength)
    with open(outputFile, 'wb') as f:
        pickle.dump(transPromoMat, f)


########################################
# 3. Segway enrichment of the labels in genomic regions (from QC_transcriptionComparison_03.py)
########################################
# this is only for the raw .bed files, where labels are only cluster numbers
def SegwayTranscriptionEnrichment(annotationFolder, annFile, RNAseqFile, extension, geneList, geneIDList, mnemFile):

        
    with open(RNAseqFile, 'rb') as pickledFile:
        expression = pickle.load(pickledFile)

    # unzip the bed file
    
    # prepare ann
    if annFile.endswith('.gz'):
        os.system('gunzip %s' %(annFile))
        annFile = annFile[0:-3]
    else:
        os.system('gunzip %s.gz' %(annFile))

    annLineCount = int(sp.getoutput('wc -l %s' %(annFile)).split()[0])

    # get the mnemonics
    label_term_mapping = {}
    with open(mnemFile, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term

    clusterCount = len(label_term_mapping)
    sortedClusterList = list(range(0,clusterCount))

    cgi = 0 # walks on the geneIDList
    ann_start = 0 # the start of the current annotation
    ann_end = 0 # the end of the current annotation
    ann_chr = 'chr'
    ann_line_count = 0 # this is just to check the progress through the annotation file
        
    clusterMats = [np.zeros((clusterCount, 160)),np.zeros((clusterCount, 160)),np.zeros((clusterCount, 160))]

    annLineInd = 1 # this keeps the annotation line index, the first line is empty, thus the index starts at 1
    genomic_region_start_annLineInd = 0

    previous_gene_chr = 'chr1'
    previous_extension_end = 0 # doesn't matter since chr condition is never true until the first gene is processed
    previous_ann_chr = 'chr1'

    while cgi < 26017: # >>>>>>>>>> MAIN
            
        #print('cgi')
        print(cgi)

        'we are in the gene territory, we are walking on the genes'
   
        geneID = geneIDList[cgi]
        gene_chr = geneList[geneIDList[cgi]].chrom
            
        gene_strand = geneList[geneID].strand
        gene_start = geneList[geneID].start
        gene_end = geneList[geneID].end
        gene_length = gene_end - gene_start
        gene_length_unit = int(gene_length/100)
        if gene_length_unit == 0:
            cgi = cgi + 1
            continue

        gene_length_last_unit = gene_length - (99* gene_length_unit)
        # TODO: ideally: something to fix: the gene_length_last_unit for negative strand versus positive strand

        extension_start = gene_start - extension
        extension_end = gene_end + extension

        ''' picking the label/class matrix based on the gene expression level'''
        #TODO catch exception for when geneID is not in expression
        gene_exp = expression[geneID]
        if gene_exp == 0:
            expMatInd = 0
        elif gene_exp > 1:
            expMatInd = 2
        else:
            expMatInd = 1

        geneMatWalkIndex = 0
        previous_fill = 0 # at the begining of each gene, annotations are full
            
        # reading the next annotation
        line = linecache.getline(annFile, annLineInd)
        annLineInd +=1
        ann_line_count += 1
        fields = line.strip().split()

        previous_ann_chr = ann_chr
        ann_chr = fields[0]
        ann_start = int(fields[1])
        ann_end = int(fields[2])

        if (gene_chr != previous_gene_chr): # in case of chromosome change because gene moved to the next chromosome
            #print('gene chr not equal to previous gene chr')
            #print('gene_chr %s' %(gene_chr))
            #print('p_gene_chr %s' %(previous_gene_chr))
            #print('ann_chr %s' %(ann_chr))
            #print(cgi)
            #print(annLineInd)
                
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

                #print(ann_chr)
                #print(previous_ann_chr)
                #print(gene_chr)
                #print(annLineInd)
                #print(cgi)
                #print('in the IF')

                
        if (ann_chr != gene_chr): # if annotation moved to the next chromosome, but gene has not yet moved to the next chromosome
            annLineInd = annLineInd - 2
            line = linecache.getline(annFile, annLineInd)
            annLineInd +=1
            ann_line_count += 1
            fields = line.strip().split()

            ann_chr = fields[0]
            ann_start = int(fields[1])
            ann_end = int(fields[2])


            ''' 
            if this gene starts somewhere in the preivous genomic region, 
            I will go back in the annotation file to the beginning of annotation for the previous gene
            changed this th the next while loop - basically we will go back until the annotations start before the gene territory

            '''

        while (ann_start > extension_start) and (gene_chr == ann_chr): # in case of overlapping genes
                
            annLineInd = max(annLineInd - 5,1)
            line = linecache.getline(annFile, annLineInd)
            annLineInd +=1
            ann_line_count += 1
            fields = line.strip().split()

            ann_chr = fields[0]
            ann_start = int(fields[1])
            ann_end = int(fields[2])
                
            #while ((ann_start < extension_end) or not(gene_chr == ann_chr)) and geneMatWalkIndex < 160: 
        while ((ann_start < extension_end) and (gene_chr == ann_chr)) and geneMatWalkIndex < 160:
                
            '''
            NOTE: the second condition is for when we are at the end of a gene's territory, and then at the end of a chromosome, so when we go to the next gene. At some point I changed the second condition from "or" to "and" 
            ann_start is not smaller than extension_end, and we need to read the annotations until we are on the same chromosome
            The condition to be in one gene's territory (and not the next one), and it is for reading the annotations
            while in THIS gene territory, read annotations until we reach to the next genomic region (plus extension)
      
            '''

            # in the next sections we are processing the annotation
            if ann_chr == gene_chr: # if we are in the same choromosome
                if ((ann_start < extension_end and ann_start > extension_start) or 
                    (ann_end < extension_end and ann_end > extension_start) or
                    (ann_start < extension_start and ann_end > extension_end)):

                    ''' We are in the genonimc region (with extension)'''
                    
                    ''' Taking off for filling the matrices... '''
                        
                    adjusted_ann_start = max(0, ann_start - extension_start)
                    adjusted_ann_end = min(ann_end - extension_start, extension_end - extension_start)
                    adjusted_ann_length = adjusted_ann_end - adjusted_ann_start

                    ann_cluster = fields[3].split('_')[0]
                    clusterInd = int(ann_cluster)

                    # expMatInd

                    if gene_strand == '+':  # munching 100bp s from annotation, filling the geneMat
                        while(adjusted_ann_length > 0) and geneMatWalkIndex < 30:
                            if (adjusted_ann_length >= 100*(1- previous_fill)):
                                clusterMats[expMatInd][clusterInd][geneMatWalkIndex] += 1 - previous_fill
                                adjusted_ann_length -= 100*(1- previous_fill)
                                previous_fill = 0
                                geneMatWalkIndex +=1
                            else:
                                clusterMats[expMatInd][clusterInd][geneMatWalkIndex] += adjusted_ann_length/100
                                previous_fill += adjusted_ann_length/100
                                adjusted_ann_length = 0
                                if previous_fill > 1:
                                    previous_fill = 1

                        while(adjusted_ann_length > 0) and geneMatWalkIndex < 129:
                            if (adjusted_ann_length >= gene_length_unit*(1- previous_fill)):
                                clusterMats[expMatInd][clusterInd][geneMatWalkIndex] += 1 - previous_fill
                                adjusted_ann_length -= gene_length_unit*(1 - previous_fill)
                                previous_fill = 0
                                geneMatWalkIndex +=1
                            else:
                                clusterMats[expMatInd][clusterInd][geneMatWalkIndex] += adjusted_ann_length/gene_length_unit
                                previous_fill += adjusted_ann_length/gene_length_unit
                                adjusted_ann_length = 0
                                if previous_fill > 1:
                                    previous_fill = 1

                        if (adjusted_ann_length > 0) and geneMatWalkIndex == 129:
                            if (adjusted_ann_length >= gene_length_last_unit*(1- previous_fill)):
                                clusterMats[expMatInd][clusterInd][geneMatWalkIndex] += 1 - previous_fill
                                adjusted_ann_length -= gene_length_last_unit*(1 - previous_fill)
                                previous_fill = 0
                                geneMatWalkIndex +=1
                            else:
                                clusterMats[expMatInd][clusterInd][geneMatWalkIndex] += adjusted_ann_length/gene_length_last_unit
                                previous_fill += adjusted_ann_length/gene_length_last_unit
                                adjusted_ann_length = 0
                                if previous_fill > 1:
                                    previous_fill = 1

                        while(adjusted_ann_length > 0) and (geneMatWalkIndex < 160):
                            if (adjusted_ann_length >= 100*(1- previous_fill)):
                                clusterMats[expMatInd][clusterInd][geneMatWalkIndex] += 1 - previous_fill
                                adjusted_ann_length -= 100*(1- previous_fill)
                                previous_fill = 0
                                geneMatWalkIndex +=1
                            else:
                                clusterMats[expMatInd][clusterInd][geneMatWalkIndex] += adjusted_ann_length/100
                                previous_fill += adjusted_ann_length/100
                                adjusted_ann_length = 0
                                if previous_fill > 1:
                                    previous_fill = 1


                    if gene_strand == '-':

                        while(adjusted_ann_length > 0) and geneMatWalkIndex < 30:
                        
                            if (adjusted_ann_length >= 100*(1- previous_fill)):
                                clusterMats[expMatInd][clusterInd][159 - geneMatWalkIndex] += 1 - previous_fill
                                adjusted_ann_length -= 100*(1- previous_fill)
                                previous_fill = 0
                                geneMatWalkIndex +=1
                            else:
                                clusterMats[expMatInd][clusterInd][159 - geneMatWalkIndex] += adjusted_ann_length/100
                                previous_fill += adjusted_ann_length/100
                                adjusted_ann_length = 0
                                if previous_fill > 1:
                                    previous_fill = 1

                        if (adjusted_ann_length > 0) and geneMatWalkIndex == 30:
                            if (adjusted_ann_length >= gene_length_last_unit*(1- previous_fill)):
                                clusterMats[expMatInd][clusterInd][159 - geneMatWalkIndex] += 1 - previous_fill
                                adjusted_ann_length -= gene_length_last_unit*(1 - previous_fill)
                                previous_fill = 0
                                geneMatWalkIndex +=1
                            else:
                                clusterMats[expMatInd][clusterInd][159 - geneMatWalkIndex] += adjusted_ann_length/gene_length_last_unit
                                previous_fill += adjusted_ann_length/gene_length_last_unit
                                adjusted_ann_length = 0
                                if previous_fill > 1:
                                    previous_fill = 1


                        while(adjusted_ann_length > 0) and geneMatWalkIndex < 129:
                            if (adjusted_ann_length >= gene_length_unit*(1- previous_fill)):
                                clusterMats[expMatInd][clusterInd][159 - geneMatWalkIndex] += 1 - previous_fill
                                adjusted_ann_length -= gene_length_unit*(1 - previous_fill)
                                previous_fill = 0
                                geneMatWalkIndex +=1
                            else:
                                clusterMats[expMatInd][clusterInd][159 - geneMatWalkIndex] += adjusted_ann_length/gene_length_unit
                                previous_fill += adjusted_ann_length/gene_length_unit
                                adjusted_ann_length = 0
                                if previous_fill > 1:
                                    previous_fill = 1


                        while(adjusted_ann_length > 0) and (geneMatWalkIndex < 160):
                            if (adjusted_ann_length >= 100*(1- previous_fill)):
                                clusterMats[expMatInd][clusterInd][159 - geneMatWalkIndex] += 1 - previous_fill
                                adjusted_ann_length -= 100*(1- previous_fill)
                                previous_fill = 0
                                geneMatWalkIndex +=1
                            else:
                                clusterMats[expMatInd][clusterInd][159 - geneMatWalkIndex] += adjusted_ann_length/100
                                previous_fill += adjusted_ann_length/100
                                adjusted_ann_length = 0
                                if previous_fill > 1:
                                    previous_fill = 1


            if geneMatWalkIndex < 160: # we read annotations until we cover the genee

                line = linecache.getline(annFile, annLineInd)
                annLineInd +=1
                ann_line_count += 1
                fields = line.strip().split()

                ann_chr = fields[0]
                ann_start = int(fields[1])
                ann_end = int(fields[2])

            #if geneMatWalkIndex >= 159:
                #print('gene Done')
                #print(geneMatWalkIndex)


        cgi += 1 # next gene
        previous_extension_end = extension_end
        previous_gene_chr = gene_chr
            

    linecache.clearcache()

    os.system('gzip %s' %(annFile))

    expOverlapMats = {"clusterMats": clusterMats}
    outputFile = annotationFolder + 'defaultExp_5kg_expSummary_newSort.pkl'
    with open(outputFile, 'wb') as f:
        pickle.dump(expOverlapMats, f)


########################################
# 4. Segway enrichment for labels in genomic regions (no transcript file required)
########################################
# this is only for the raw .bed files, where labels are only cluster numbers
def SegwayGeneBodyEnrichment(annotationFolder, annFile, extension, geneList, geneIDList, mnemFile):

    # unzip the bed file
    # prepare ann
    if annFile.endswith('.gz'):
        os.system('gunzip %s' %(annFile))
        annFile = annFile[0:-3]
    else:
        os.system('gunzip %s.gz' %(annFile))

    annLineCount = int(sp.getoutput('wc -l %s' %(annFile)).split()[0])

    # get the mnemonics
    label_term_mapping = {}
    with open(mnemFile, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term

    clusterCount = len(label_term_mapping)
    sortedClusterList = list(range(0,clusterCount))

    cgi = 0 # walks on the geneIDList
    ann_start = 0 # the start of the current annotation
    ann_end = 0 # the end of the current annotation
    ann_chr = 'chr'
    ann_line_count = 0 # this is just to check the progress through the annotation file

    geneMat = np.zeros((clusterCount, 160))

    annLineInd = 1 # this keeps the annotation line index, the first line is empty, thus the index starts at 1
    genomic_region_start_annLineInd = 0

    previous_gene_chr = 'chr1'
    previous_extension_end = 0 # doesn't matter since chr condition is never true until the first gene is processed
    previous_ann_chr = 'chr1'

    while cgi < 26017: # >>>>>>>>>> MAIN
            
        #print('cgi')

        'we are in the gene territory, we are walking on the genes'
   
        geneID = geneIDList[cgi]
        gene_chr = geneList[geneIDList[cgi]].chrom
            
        gene_strand = geneList[geneID].strand
        gene_start = geneList[geneID].start
        gene_end = geneList[geneID].end
        gene_length = gene_end - gene_start
        gene_length_unit = int(gene_length/100)
        if gene_length_unit == 0:
            cgi = cgi + 1
            continue

        gene_length_last_unit = gene_length - (99* gene_length_unit)
        # TODO: ideally: something to fix: the gene_length_last_unit for negative strand versus positive strand

        extension_start = gene_start - extension
        extension_end = gene_end + extension

        ''' picking the label/class matrix based on the gene expression level'''

        geneMatWalkIndex = 0
        previous_fill = 0 # at the begining of each gene, annotations are full
            
        # reading the next annotation
        line = linecache.getline(annFile, annLineInd)
        annLineInd +=1
        ann_line_count += 1
        fields = line.strip().split()

        previous_ann_chr = ann_chr
        ann_chr = fields[0]
        ann_start = int(fields[1])
        ann_end = int(fields[2])

        if (gene_chr != previous_gene_chr): # in case of chromosome change because gene moved to the next chromosome
            #print('gene chr not equal to previous gene chr')
            #print('gene_chr %s' %(gene_chr))
            #print('p_gene_chr %s' %(previous_gene_chr))
            #print('ann_chr %s' %(ann_chr))
            #print(cgi)
            #print(annLineInd)
                
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

                #print(ann_chr)
                #print(previous_ann_chr)
                #print(gene_chr)
                #print(annLineInd)
                #print(cgi)
                #print('in the IF')

                
        if (ann_chr != gene_chr): # if annotation moved to the next chromosome, but gene has not yet moved to the next chromosome
            annLineInd = annLineInd - 2
            line = linecache.getline(annFile, annLineInd)
            annLineInd +=1
            ann_line_count += 1
            fields = line.strip().split()

            ann_chr = fields[0]
            ann_start = int(fields[1])
            ann_end = int(fields[2])


            ''' 
            if this gene starts somewhere in the preivous genomic region, 
            I will go back in the annotation file to the beginning of annotation for the previous gene
            changed this th the next while loop - basically we will go back until the annotations start before the gene territory

            '''

        while (ann_start > extension_start) and (gene_chr == ann_chr): # in case of overlapping genes
                
            annLineInd = max(annLineInd - 5,1)
            line = linecache.getline(annFile, annLineInd)
            annLineInd +=1
            ann_line_count += 1
            fields = line.strip().split()

            ann_chr = fields[0]
            ann_start = int(fields[1])
            ann_end = int(fields[2])
                
            #while ((ann_start < extension_end) or not(gene_chr == ann_chr)) and geneMatWalkIndex < 160: 
        while ((ann_start < extension_end) and (gene_chr == ann_chr)) and geneMatWalkIndex < 160:
                
            '''
            NOTE: the second condition is for when we are at the end of a gene's territory, and then at the end of a chromosome, so when we go to the next gene. At some point I changed the second condition from "or" to "and" 
            ann_start is not smaller than extension_end, and we need to read the annotations until we are on the same chromosome
            The condition to be in one gene's territory (and not the next one), and it is for reading the annotations
            while in THIS gene territory, read annotations until we reach to the next genomic region (plus extension)
      
            '''

            # in the next sections we are processing the annotation
            if ann_chr == gene_chr: # if we are in the same choromosome
                if ((ann_start < extension_end and ann_start > extension_start) or 
                    (ann_end < extension_end and ann_end > extension_start) or
                    (ann_start < extension_start and ann_end > extension_end)):

                    ''' We are in the genonimc region (with extension)'''
                    
                    ''' Taking off for filling the matrices... '''
                        
                    adjusted_ann_start = max(0, ann_start - extension_start)
                    adjusted_ann_end = min(ann_end - extension_start, extension_end - extension_start)
                    adjusted_ann_length = adjusted_ann_end - adjusted_ann_start

                    ann_cluster = fields[3].split('_')[0]
                    clusterInd = int(ann_cluster)

                    # expMatInd

                    if gene_strand == '+':  # munching 100bp s from annotation, filling the geneMat
                        while(adjusted_ann_length > 0) and geneMatWalkIndex < 30:
                            if (adjusted_ann_length >= 100*(1- previous_fill)):
                                geneMat[clusterInd][geneMatWalkIndex] += 1 - previous_fill
                                adjusted_ann_length -= 100*(1- previous_fill)
                                previous_fill = 0
                                geneMatWalkIndex +=1
                            else:
                                geneMat[clusterInd][geneMatWalkIndex] += adjusted_ann_length/100
                                previous_fill += adjusted_ann_length/100
                                adjusted_ann_length = 0
                                if previous_fill > 1:
                                    previous_fill = 1

                        while(adjusted_ann_length > 0) and geneMatWalkIndex < 129:
                            if (adjusted_ann_length >= gene_length_unit*(1- previous_fill)):
                                geneMat[clusterInd][geneMatWalkIndex] += 1 - previous_fill
                                adjusted_ann_length -= gene_length_unit*(1 - previous_fill)
                                previous_fill = 0
                                geneMatWalkIndex +=1
                            else:
                                geneMat[clusterInd][geneMatWalkIndex] += adjusted_ann_length/gene_length_unit
                                previous_fill += adjusted_ann_length/gene_length_unit
                                adjusted_ann_length = 0
                                if previous_fill > 1:
                                    previous_fill = 1

                        if (adjusted_ann_length > 0) and geneMatWalkIndex == 129:
                            if (adjusted_ann_length >= gene_length_last_unit*(1- previous_fill)):
                                geneMat[clusterInd][geneMatWalkIndex] += 1 - previous_fill
                                adjusted_ann_length -= gene_length_last_unit*(1 - previous_fill)
                                previous_fill = 0
                                geneMatWalkIndex +=1
                            else:
                                geneMat[clusterInd][geneMatWalkIndex] += adjusted_ann_length/gene_length_last_unit
                                previous_fill += adjusted_ann_length/gene_length_last_unit
                                adjusted_ann_length = 0
                                if previous_fill > 1:
                                    previous_fill = 1

                        while(adjusted_ann_length > 0) and (geneMatWalkIndex < 160):
                            if (adjusted_ann_length >= 100*(1- previous_fill)):
                                geneMat[clusterInd][geneMatWalkIndex] += 1 - previous_fill
                                adjusted_ann_length -= 100*(1- previous_fill)
                                previous_fill = 0
                                geneMatWalkIndex +=1
                            else:
                                geneMat[clusterInd][geneMatWalkIndex] += adjusted_ann_length/100
                                previous_fill += adjusted_ann_length/100
                                adjusted_ann_length = 0
                                if previous_fill > 1:
                                    previous_fill = 1


                    if gene_strand == '-':

                        while(adjusted_ann_length > 0) and geneMatWalkIndex < 30:
                        
                            if (adjusted_ann_length >= 100*(1- previous_fill)):
                                geneMat[clusterInd][159 - geneMatWalkIndex] += 1 - previous_fill
                                adjusted_ann_length -= 100*(1- previous_fill)
                                previous_fill = 0
                                geneMatWalkIndex +=1
                            else:
                                geneMat[clusterInd][159 - geneMatWalkIndex] += adjusted_ann_length/100
                                previous_fill += adjusted_ann_length/100
                                adjusted_ann_length = 0
                                if previous_fill > 1:
                                    previous_fill = 1

                        if (adjusted_ann_length > 0) and geneMatWalkIndex == 30:
                            if (adjusted_ann_length >= gene_length_last_unit*(1- previous_fill)):
                                geneMat[clusterInd][159 - geneMatWalkIndex] += 1 - previous_fill
                                adjusted_ann_length -= gene_length_last_unit*(1 - previous_fill)
                                previous_fill = 0
                                geneMatWalkIndex +=1
                            else:
                                geneMat[clusterInd][159 - geneMatWalkIndex] += adjusted_ann_length/gene_length_last_unit
                                previous_fill += adjusted_ann_length/gene_length_last_unit
                                adjusted_ann_length = 0
                                if previous_fill > 1:
                                    previous_fill = 1


                        while(adjusted_ann_length > 0) and geneMatWalkIndex < 129:
                            if (adjusted_ann_length >= gene_length_unit*(1- previous_fill)):
                                geneMat[clusterInd][159 - geneMatWalkIndex] += 1 - previous_fill
                                adjusted_ann_length -= gene_length_unit*(1 - previous_fill)
                                previous_fill = 0
                                geneMatWalkIndex +=1
                            else:
                                geneMat[clusterInd][159 - geneMatWalkIndex] += adjusted_ann_length/gene_length_unit
                                previous_fill += adjusted_ann_length/gene_length_unit
                                adjusted_ann_length = 0
                                if previous_fill > 1:
                                    previous_fill = 1


                        while(adjusted_ann_length > 0) and (geneMatWalkIndex < 160):
                            if (adjusted_ann_length >= 100*(1- previous_fill)):
                                geneMat[clusterInd][159 - geneMatWalkIndex] += 1 - previous_fill
                                adjusted_ann_length -= 100*(1- previous_fill)
                                previous_fill = 0
                                geneMatWalkIndex +=1
                            else:
                                geneMat[clusterInd][159 - geneMatWalkIndex] += adjusted_ann_length/100
                                previous_fill += adjusted_ann_length/100
                                adjusted_ann_length = 0
                                if previous_fill > 1:
                                    previous_fill = 1


            if geneMatWalkIndex < 160: # we read annotations until we cover the genee

                line = linecache.getline(annFile, annLineInd)
                annLineInd +=1
                ann_line_count += 1
                fields = line.strip().split()

                ann_chr = fields[0]
                ann_start = int(fields[1])
                ann_end = int(fields[2])

            #if geneMatWalkIndex >= 159:
                #print('gene Done')
                #print(geneMatWalkIndex)


        cgi += 1 # next gene
        previous_extension_end = extension_end
        previous_gene_chr = gene_chr
            

    linecache.clearcache()

    os.system('gzip %s' %(annFile))

    outputFile = annotationFolder + 'default_5kg_geneBodyCoverage_newSort.pkl'
    with open(outputFile, 'wb') as f:
        pickle.dump(geneMat, f)




