# This I am getting from QC_transcriptionComparison.py, but just the transcription file running
# This code is different than the original in that it works for a limited region of the genome, and for cluster labels only
#
# Segments:
# 0. Initials
# 1. Specifics
# 2. Main


########################################
# 0. Initials 
########################################

import linecache

# General data folder
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'

# GTF data structure
fileName = dataFolder + '/geneLists.pkl'
with open(fileName, 'rb') as pickleFile:
    geneListsAndIDs = pickle.load(pickleFile)
    
geneList = geneListsAndIDs[0]
geneIDList = geneListsAndIDs[1]
del geneListsAndIDs

segwayLabels = ['Quiescent', 'ConstitutiveHet', 'FacultativeHet', 'Transcribed', 'Promoter', 'Enhancer', 'RegPermissive', 'Bivalent', 'LowConfidence']

########################################
# 1. Specifics
########################################

# .bed files - replicates
bedFile1 = dataFolder + 'Mehdi_testRun/segway_1.bed'
bedFile2 = dataFolder + 'Mehdi_testRun/segway_2.bed' 


# expression file - just two expresion for GM12878
expFile1 = dataFolder + 'Mehdi_testRun/ENCFF345SHY.tsv'
expFile1 = dataFolder + 'Mehdi_testRun/ENCFF240WBI.tsv'


########################################
# 2. Pre processings
########################################

# preping expression data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# file1        
expression1 = {}
with open(expFile1, 'r') as expFile:
    line = expFile.readline() # there is one header line

    for line in expFile:
        fields = line.strip().split()
        geneID = fields[0]
        transcriptID = fields[1]

        if geneID in expression1:
            expression1[geneID] += np.log10(float(fields[5]) + 1)
        else:
            expression1[geneID] = np.log10(float(fields[5]) + 1)


fileName = dataFolder + 'Mehdi_testRun/geneExp_dict_ENCFF345SHY.pkl'
with open(fileName, 'wb') as f:
    pickle.dump(expression1, f)

# loading
fileName = dataFolder + 'Mehdi_testRun/geneExp_dict_ENCFF240WBI.pkl'
with open(fileName, 'rb') as pickledFile:
    expression2 = pickle.load(pickledFile)


# file2
expression2 = {}
with open(expFile1, 'r') as expFile:
    line = expFile.readline() # there is one header line

    for line in expFile:
        fields = line.strip().split()
        geneID = fields[0]
        transcriptID = fields[1]

        if geneID in expression2:
            expression2[geneID] += np.log10(float(fields[5]) + 1)
        else:
            expression2[geneID] = np.log10(float(fields[5]) + 1)


fileName = dataFolder + 'Mehdi_testRun/geneExp_dict_ENCFF240WBI.pkl'
with open(fileName, 'wb') as f:
    pickle.dump(expression2, f)

    
# loading classes and labels
fileName = dataFolder + 'Mehdi_testRun/geneExp_dict_ENCFF240WBI.pkl'
with open(fileName, 'rb') as pickledFile:
    expression2 = pickle.load(pickledFile)


# preping annotation data
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

annotation_generalInfo_clusters(bedFile2)

def annotation_generalInfo_clusters(bedFile):

    splitFileName = bedFile.split('/')
    fileName = splitFileName[-1].split('.')[0]
    index = bedFile.index(fileName)
    outputFolder = bedFile[0:index]

    clusters = {}
    with open(bedFile, 'r') as annotations:

        header = annotations.readline()
        
        for line in annotations:
            fields = line.strip().split()
            if fields[3] in clusters.keys():

                clusters[fields[3]].bp_count += int(fields[2]) - int(fields[1])
                clusters[fields[3]].region_count += 1
                
            else:
                
                called = fields[3]
                biolabel = ''
                cluster = ''
                color = ''
                bp_count = int(fields[2]) - int(fields[1])
                region_count = 1
                region_dist = 1
                clusters[fields[3]] = Annotation(called, biolabel, cluster, color, bp_count, region_count, region_dist)

    fileName = outputFolder + fileName + '_annotationSummary.pkl'
    with open(fileName, 'wb') as f:
        pickle.dump(clusters, f)

    print('annotation summary saved in %s' %(fileName))


########################################
# 3. Main
########################################

# 3.0 getting the genes which are in the pilot region
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

annLineInd = 1 # this keeps the annotation line index, the first line is empty, thus the index starts at 1
genomic_region_start_annLineInd = 0
firstGenomicAnn = True

previous_gene_chr = 'chr'
previous_extension_end = 0 # doesn't matter since chr condition is never true until the first gene is processe

line = linecache.getline(annFile, annLineInd)
fields = line.strip().split()

ann_chr = fields[0]
ann_start = int(fields[1])
ann_end = int(fields[2])

#while cgi < len(geneList):

chrlist = []
for i in range(22):
    chrlist.append('chr' + str(i+1))

while cgi < 10 and (ann_chr in chrlist) and (gene_chr in chrlist)

   firstGenomicAnn = True
   
   geneID = geneIDList[cgi]
   gene_chr = geneList[geneIDList[cgi]].chrom
   print(geneID)

   gene_strand = geneList[geneID].strand
   gene_start = geneList[geneID].start
   gene_end = geneList[geneID].end
   
   extension_start = gene_start - extension
   extension_end = gene_end + extension
   extension_length = extension_end - extension_start

   ''' 
   if this gene starts somewhere in the preivous genomic region, 
   I will go back in the annotation file to the beginning of annotation for the previous gene

   '''
   if (gene_chr == previous_gene_chr) and (previous_extension_end > extension_start):
      annLineInd = genomicRegionStartAnnLineInd
      # TODO here read the new annotation


   # if annotation is ahead, go to the next gene
   if (ann_start > extension_end) or (int(ann_chr[3:] > int(gene_chr[3:]))):
       cgi += 1

   else:
       genecover = 0
       while (genecover < gene_length) or not((ann_start > extension_end) or (int(ann_chr[3:] > int(gene_chr[3:])))):
           ''' while we haven't walked on the gene or annotations have not passed the gene mov on in annotation'''
           
           annLineInd +=1
           line = linecache.getline(annFile, annLineInd)
           fields = line.strip().split()
           
           ann_chr = fields[0]
           ann_start = int(fields[1])
           ann_end = int(fields[2])

           if ann_chr == gene_chr: # if we are in the genomic region (with extension)
               if ((ann_start < extension_end and ann_start > extension_start) or 
                   (ann_end < extension_end and ann_end > extension_start) or
                   (ann_start < extension_start and ann_end > extension_end)):
                   ''' We are in the genonimc region (with extension)'''
            
                   if firstGenomicAnn:
                       ''' Keeping the line index of the first annotation of the genomic region for overlapping genes'''
                       genomicRegionStartAnnLineInd = annLineInd -1
                       firstGenomicAnn = False


                       ''' Taking off for filling the matrices... '''
                       adjusted_ann_start = max(0, ann_start - extension_start)
                       adjusted_ann_end = min(ann_end - extension_start, extension_end - extension_start)
                       adjusted_ann_length = adjusted_ann_end - adjusted_ann_start

                       genecover += adjusted_ann_end - adjusted_ann_start 

        if genecover > gene_length:
            pilotGeneList.append(geneIDList[cgi])

        cgi += 1

# 3.1 doing it with the selected gene list      
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      
'''

This is different than the main pipeline, since the annotations cover only the pilot region

'''

annotationFileSummary = 'segway_1_annotationSummary.pkl'
fileName = dataFolder + 'Mehdi_testRun/' + annotationFileSummary
with open(fileName, 'rb') as pickledFile:
    clusters = pickle.load(pickledFile)
labels = clusters
labelList = list(labels.keys())
labelCount = len(labelList)
labelListIndMap = {}
for i in range(len(labelList)):
   labelListIndMap[labelList[i]] = i


expressionFile = 'geneExp_dict_ENCFF240WBI.pkl'
fileName = dataFolder + 'Mehdi_testRun/' + expressionFile
with open(fileName, 'rb') as pickledFile:
    expression = pickle.load(pickledFile)

    
extension = 3000 # count of basepairs monitored before and after the gene coordinates
cgi = 0 # walks on the geneIDList
ann_start = 0 # the start of the current annotation
ann_end = 0 # the end of the current annotation
ann_chr = 'chr'
ann_line_count = 0 # this is just to check the progress through the annotation file
previous_class = ''


labelMats = [np.zeros((labelCount, 160)),np.zeros((labelCount, 160)),np.zeros((labelCount, 160))]
classMats = [np.zeros((classCount, 160)),np.zeros((classCount, 160)),np.zeros((classCount, 160))]

# this is to use for the negative strand genes - not using now
#tempLabelMats = [np.zeros((labelCount, 160)),np.zeros((labelCount, 160)),np.zeros((labelCount, 160))]
#tempClassMats = [np.zeros((classCount, 160)),np.zeros((classCount, 160)),np.zeros((classCount, 160))]

annLineInd = 1 # this keeps the annotation line index, the first line is empty, thus the index starts at 1
genomic_region_start_annLineInd = 0
firstGenomicAnn = True

previous_gene_chr = 'chr'
previous_extension_end = 0 # doesn't matter since chr condition is never true until the first gene is processed

#while cgi < len(geneIDList):
while cgi < 1000:
    print(cgi)

   'we are in the gene territory'
   firstGenomicAnn = True
   
   geneID = geneIDList[cgi]
   gene_chr = geneList[geneIDList[cgi]].chrom
   print(geneID)

   gene_strand = geneList[geneID].strand
   gene_start = geneList[geneID].start
   gene_end = geneList[geneID].end
   gene_length = gene_end - gene_start
   gene_length_unit = int(gene_length/100)
   gene_length_last_unit = gene_length - (99* gene_length_unit)
   # TODO: something to fix: the gene_length_last_unit for negative strand versus positive strand

   extension_start = gene_start - extension
   extension_end = gene_end + extension

   ''' 
   if this gene starts somewhere in the preivous genomic region, 
   I will go back in the annotation file to the beginning of annotation for the previous gene

   '''
   if (gene_chr == previous_gene_chr) and (previous_extension_end > extension_start):
      annLineInd = genomicRegionStartAnnLineInd
   

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
            
   ann_chr = fields[0]
   ann_start = int(fields[1])
   ann_end = int(fields[2])
        
   while ((ann_start < extension_end) or not(gene_chr == ann_chr)) and geneMatWalkIndex < 160: 
           
      '''
      NOTE: the second condition is for when we are at the end of a gene's territory, and then at the end of a chromosome, so when we go to the next gene, 
      ann_start is not smaller than extension_end, and we need to read the annotations until we are on the same chromosome

      The condition to be in one gene's territory (and not the next one), and it is for reading the annotations
      while in THIS gene territory, read annotations until we reach to the genomic region (plus extension)
      
      '''

      #  print('hereIam') #>>>>>>>>>> test
      #  line = annotations.readline()
            

      if ann_chr == gene_chr: # if we are in the genomic region (with extension)
         if ((ann_start < extension_end and ann_start > extension_start) or 
             (ann_end < extension_end and ann_end > extension_start) or
             (ann_start < extension_start and ann_end > extension_end)):

                   
            ''' We are in the genonimc region (with extension)'''
            
            if firstGenomicAnn:
               ''' Keeping the line index of the first annotation of the genomic region for overlapping genes'''
               genomicRegionStartAnnLineInd = annLineInd -1
               firstGenomicAnn = False


            ''' Taking off for filling the matrices... '''
            adjusted_ann_start = max(0, ann_start - extension_start)
            adjusted_ann_end = min(ann_end - extension_start, extension_end - extension_start)
            adjusted_ann_length = adjusted_ann_end - adjusted_ann_start

            ann_label = fields[3]
            labelInd = labelListIndMap(ann_label)
            
            ann_class = ann_label.split('_')[1]
            classInd = classListIndMap[ann_class]
            
            # expMatInd
            if gene_strand == '+':
               while(adjusted_ann_length > 0) and geneMatWalkIndex < 30:
                  if (adjusted_ann_length >= 100*(1- previous_fill)):
                     classMats[expMatInd][classInd][geneMatWalkIndex] += 1 - previous_fill
                     labelMats[expMatInd][labelInd][geneMatWalkIndex] += 1 - previous_fill
                     adjusted_ann_length -= 100*(1- previous_fill)
                     previous_fill = 0
                     geneMatWalkIndex +=1
                  else:
                     classMats[expMatInd][classInd][geneMatWalkIndex] += adjusted_ann_length/100
                     labelMats[expMatInd][labelInd][geneMatWalkIndex] += adjusted_ann_length/100
                     previous_fill += adjusted_ann_length/100
                     adjusted_ann_length = 0
                     if previous_fill > 1:
                        previous_fill = 1

               while(adjusted_ann_length > 0) and geneMatWalkIndex < 129:
                  if (adjusted_ann_length >= gene_length_unit*(1- previous_fill)):
                     classMats[expMatInd][classInd][geneMatWalkIndex] += 1 - previous_fill
                     labelMats[expMatInd][labelInd][geneMatWalkIndex] += 1 - previous_fill
                     adjusted_ann_length -= gene_length_unit*(1 - previous_fill)
                     previous_fill = 0
                     geneMatWalkIndex +=1
                  else:
                     classMats[expMatInd][classInd][geneMatWalkIndex] += adjusted_ann_length/gene_length_unit
                     labelMats[expMatInd][labelInd][geneMatWalkIndex] += adjusted_ann_length/gene_length_unit
                     previous_fill += adjusted_ann_length/gene_length_unit
                     adjusted_ann_length = 0
                     if previous_fill > 1:
                        previous_fill = 1

               if (adjusted_ann_length > 0) and geneMatWalkIndex == 129:
                     if (adjusted_ann_length >= gene_length_last_unit*(1- previous_fill)):
                        classMats[expMatInd][classInd][geneMatWalkIndex] += 1 - previous_fill
                        labelMats[expMatInd][labelInd][geneMatWalkIndex] += 1 - previous_fill
                        adjusted_ann_length -= gene_length_last_unit*(1 - previous_fill)
                        previous_fill = 0
                        geneMatWalkIndex +=1
                     else:
                        classMats[expMatInd][classInd][geneMatWalkIndex] += adjusted_ann_length/gene_length_last_unit
                        labelMats[expMatInd][labelInd][geneMatWalkIndex] += adjusted_ann_length/gene_length_last_unit
                        previous_fill += adjusted_ann_length/gene_length_last_unit
                        adjusted_ann_length = 0
                        if previous_fill > 1:
                           previous_fill = 1
                         
               while(adjusted_ann_length > 0) and (geneMatWalkIndex < 160):
                  if (adjusted_ann_length >= 100*(1- previous_fill)):
                     classMats[expMatInd][classInd][geneMatWalkIndex] += 1 - previous_fill
                     labelMats[expMatInd][labelInd][geneMatWalkIndex] += 1 - previous_fill
                     adjusted_ann_length -= 100*(1- previous_fill)
                     previous_fill = 0
                     geneMatWalkIndex +=1
                  else:
                     classMats[expMatInd][classInd][geneMatWalkIndex] += adjusted_ann_length/100
                     labelMats[expMatInd][labelInd][geneMatWalkIndex] += adjusted_ann_length/100
                     previous_fill += adjusted_ann_length/100
                     adjusted_ann_length = 0
                     if previous_fill > 1:
                        previous_fill = 1


            if gene_strand == '-':
               while(adjusted_ann_length > 0) and geneMatWalkIndex < 30:
                  if (adjusted_ann_length >= 100*(1- previous_fill)):
                     classMats[expMatInd][classInd][159 - geneMatWalkIndex] += 1 - previous_fill
                     labelMats[expMatInd][labelInd][159 - geneMatWalkIndex] += 1 - previous_fill
                     adjusted_ann_length -= 100*(1- previous_fill)
                     previous_fill = 0
                     geneMatWalkIndex +=1
                  else:
                     classMats[expMatInd][classInd][159 - geneMatWalkIndex] += adjusted_ann_length/100
                     labelMats[expMatInd][labelInd][159 - geneMatWalkIndex] += adjusted_ann_length/100
                     previous_fill += adjusted_ann_length/100
                     adjusted_ann_length = 0
                     if previous_fill > 1:
                        previous_fill = 1
                        
               if (adjusted_ann_length > 0) and geneMatWalkIndex == 30:
                     if (adjusted_ann_length >= gene_length_last_unit*(1- previous_fill)):
                        classMats[expMatInd][classInd][159 - geneMatWalkIndex] += 1 - previous_fill
                        labelMats[expMatInd][labelInd][159 - geneMatWalkIndex] += 1 - previous_fill
                        adjusted_ann_length -= gene_length_last_unit*(1 - previous_fill)
                        previous_fill = 0
                        geneMatWalkIndex +=1
                     else:
                        classMats[expMatInd][classInd][159 - geneMatWalkIndex] += adjusted_ann_length/gene_length_last_unit
                        labelMats[expMatInd][labelInd][159 - geneMatWalkIndex] += adjusted_ann_length/gene_length_last_unit
                        previous_fill += adjusted_ann_length/gene_length_last_unit
                        adjusted_ann_length = 0
                        if previous_fill > 1:
                           previous_fill = 1


               while(adjusted_ann_length > 0) and geneMatWalkIndex < 129:
                  if (adjusted_ann_length >= gene_length_unit*(1- previous_fill)):
                     classMats[expMatInd][classInd][159 - geneMatWalkIndex] += 1 - previous_fill
                     labelMats[expMatInd][labelInd][159 - geneMatWalkIndex] += 1 - previous_fill
                     adjusted_ann_length -= gene_length_unit*(1 - previous_fill)
                     previous_fill = 0
                     geneMatWalkIndex +=1
                  else:
                     classMats[expMatInd][classInd][159 - geneMatWalkIndex] += adjusted_ann_length/gene_length_unit
                     labelMats[expMatInd][labelInd][159 - geneMatWalkIndex] += adjusted_ann_length/gene_length_unit
                     previous_fill += adjusted_ann_length/gene_length_unit
                     adjusted_ann_length = 0
                     if previous_fill > 1:
                        previous_fill = 1

                         
               while(adjusted_ann_length > 0) and (geneMatWalkIndex < 160):
                  if (adjusted_ann_length >= 100*(1- previous_fill)):
                     classMats[expMatInd][classInd][159 - geneMatWalkIndex] += 1 - previous_fill
                     labelMats[expMatInd][labelInd][159 - geneMatWalkIndex] += 1 - previous_fill
                     adjusted_ann_length -= 100*(1- previous_fill)
                     previous_fill = 0
                     geneMatWalkIndex +=1
                  else:
                     classMats[expMatInd][classInd][159 - geneMatWalkIndex] += adjusted_ann_length/100
                     labelMats[expMatInd][labelInd][159 - geneMatWalkIndex] += adjusted_ann_length/100
                     previous_fill += adjusted_ann_length/100
                     adjusted_ann_length = 0
                     if previous_fill > 1:
                        previous_fill = 1


                # if gene strand is positive, or negative, flag which label list we use
                # the only difference between the two strands is that I am going to reverse the index of
                # the columns in the matrix

      if geneMatWalkIndex < 160:
         #print(geneMatWalkIndex)
         #print(annLineInd)
         line = linecache.getline(annFile, annLineInd)
         annLineInd +=1
         ann_line_count += 1
         fields = line.strip().split()

         ann_chr = fields[0]
         ann_start = int(fields[1])
         ann_end = int(fields[2])
        

   cgi += 1 # next gene
   previous_extension_end = extension_end
   previous_gene_chr = gene_chr
   ###
   #TODO: we are probably missing something at the end of the last gene
   ###



