# To get the transcript evaluation for Mehdi's work
###################################################
# 0. Initials
###################################################

from urllib.request import urlopen
import urllib
import json
import http.cookiejar
import requests
import re
import glob
import pickle
from QC_transcriptionComparison_util import Gene, Exon, Annotation, AnnotationClass
# TODO: still have the problem that can not import functions from QC_transcriptionComparison_util

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/Mehdi_testRun/tomarjan/'

###################################################
# 0. Initials
###################################################

# from each folder, get the tissue name, get the transcript data from the portal for that tissue

# get the list of folders:

sampleList = os.listdir(dataFolder + 'segway')

# make the annotation summary for Segway files:
for item in sampleList:
    bedFile = dataFolder + 'segway/' + item + '/segway.bed'

    if 'concat' in bedFile:
        continue
    print(bedFile)
    annotation_generalInfo_clusters_noclass(bedFile)

# sample code for checking annotation summary data
file = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/Mehdi_testRun/tomarjan/segway/HeLa-S3_rep1_rs7/segway_annotationSummary.pkl'
with open(file, 'rb') as input:
    clusters = pickle.load(input)

# get a list of tissues:
tissue_list = []
for item in sampleList:
    
    # skipping the concat runs
    # where to get the transcript data
    tissue = item.split('_')[0]
    if tissue == 'CD14-positive': # we don't have transcript for this one at the moment
        tissue = 'CD14-positive-monocyte'
    print(tissue)
    if tissue in tissue_list:
        continue
    else:
        tissue_list.append(tissue)


# make folders for list of tissues for transcript data
os.mkdir(dataFolder + 'RNA_seq')

for tissue in tissue_list:
    os.mkdir(dataFolder + 'RNA_seq/' + tissue)

cookies = {'required_cookie': 'values'}
headers = {'accept' : 'application/json'}
dataSubFolder = 'segway/'

extension = 3000

for tissue in tissue_list:

    # make a folder and stuff
    
    # skipping the concat runs
    #if 'concat' in item or 'rep2' in item:
    #    continue

    # where to get the transcript data
    #tissue = item.split('_')[0]
    #if tissue == 'CD14-positive': # we don't have transcript for this one at the moment
    #    tissue = 'CD14-positive-monocyte'

    term_name = tissue
    term_name_mod = term_name.replace('-', '+')

    search_url = 'https://www.encodeproject.org/search/?type=Experiment&biosample_ontology.term_name=%s&status=released&assay_title=total+RNA-seq' %(term_name_mod)

    response = requests.get(search_url, cookies=cookies, headers=headers)
    search_results = response.json()

    tcount = 0


    print('count of RNA-seq files: %d' %(len(search_results['@graph'])))

    for result in search_results['@graph']:
        if tcount == 3:
            break


        if 'NOT_COMPLIANT' in result['audit']:
            continue
            
        RNA_accession = result['accession']
        print('annotation accession' + accession)
        print(RNA_accession)
        
        # go to the RNA page and download the gene quantity files in the folder
        #            rna_url = 'https://www.encodeproject.org/experiments/ENCSR052FJA/'
        rna_url = 'https://www.encodeproject.org/experiments/%s' %(RNA_accession)
        response = requests.get(rna_url, cookies=cookies, headers=headers)
        rna_json = response.json()

        rna_fileList = rna_json['files']

        for file in rna_fileList:
            if tcount < 3:
                print('here')
                if file['output_type'] == 'gene quantifications':
                    
                    file_accession = file['accession']
                    downloadURL = 'https://www.encodeproject.org/files/%s/@@download/%s.tsv' %(file_accession, file_accession)
                    if 'preferred_default' in file:
                        print('there')
                        if file['preferred_default'] == True:
                            downloadPath = dataFolder + 'RNA_seq/' + tissue + '/preferred_default_' + file_accession + '.tsv'
                            print('doing the download')
                            urllib.request.urlretrieve(downloadURL, downloadPath)
                            tcount += 1
            else:
                break
            

# make the gene expression file

tissue = tissue_list[-1]
expFile_list = glob.glob(dataFolder + 'RNA_seq/' + tissue + '/preferred_default_*.tsv')

for expFile in expFile_list:
    expression = {}
    with open(expFile, 'r') as file:
        line = file.readline() # there is one header line
            
        for line in file:
            fields = line.strip().split()
            geneID = fields[0]
            transcriptID = fields[1]
                
            if geneID in expression:
                expression[geneID] += np.log10(float(fields[5]) + 1)
            else:
                expression[geneID] = np.log10(float(fields[5]) + 1)
                    
    # saving expression dict            
    expAccession = re.split('_|\.', expFile)[4]
    print(expAccession)
    outputFile = dataFolder + 'RNA_seq/' + tissue + '/geneExp_dict_' + expAccession + '.pkl'
    print('printing this file %s' %(outputFile))
         
    with open(outputFile, 'wb') as f:
        pickle.dump(expression, f)

# prepping the data: we only need folders starting the three tissues that we have
inTissues = []
inTissues.append(tissue_list[0])
inTissues.append(tissue_list[2])
inTissues.append(tissue_list[4])

sampleSubList = []
for item in sampleList:
    for tissue in inTissues:
        if tissue in item and not('concat' in item):
            sampleSubList.append(item)

sampleSubList = ['GM12878_rep1', 'K562_rep1_rs5', 'MCF-7_rep1', 'MCF-7_rep1_rs5', 'GM12878_rep1_rs7', 'K562_rep2', 'GM12878_rep2', 'MCF-7_rep2', 'MCF-7_rep1_rs7', 'K562_rep1_rs7', 'K562_rep1', 'GM12878_rep1_rs5']

for item in sampleSubList:
    
    print(item)
    sampleFolderAdd = dataFolder + 'segway/' + item + '/'
    bedFileName = 'segway.bed'

    bedFileAdd = sampleFolderAdd + bedFileName
    
    # sort the bed file
    sortedBedFile = sampleFolderAdd + 'sortedBedFile.bed'
    os.system('sort -V -k1,1 -k2,2 %s > %s' %(bedFileAdd, sortedBedFile))

    # filter the bed file
    filteredSortedBedFile = sampleFolderAdd + 'segway_filteredSorted.bed'
    os.system("awk -F'\t' '$1 !~ \"_\"' %s > %s" %(sortedBedFile, filteredSortedBedFile))

    os.remove(sortedBedFile)

# the actual comparison with transcriptomic data

# GTF data structure
fileName = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/geneLists.pkl'
with open(fileName, 'rb') as pickleFile:
    geneListsAndIDs = pickle.load(pickleFile)
    
geneList = geneListsAndIDs[0]
geneIDList = geneListsAndIDs[1]
del geneListsAndIDs

extension = 3000
counter = 0

tempList = sampleSubList[4:]
for item in sampleSubList:
    
for item in tempList:
    
    print(item)
    print(counter)
    counter +=1

    # get the list of exp files
    expFolder = dataFolder + 'RNA_seq/' + item.split('_')[0] + '/'
    expFileList = glob.glob(expFolder + 'geneExp_dict_*.pkl')

    # get the add of the bed file
    sampleFolderAdd = dataFolder + 'segway/' + item + '/'  # where the .bed file is
    bedFileName = 'segway_filteredSorted.bed'
    annFile = sampleFolderAdd + bedFileName

    # get the annotation summary
    inputFile = sampleFolderAdd + 'segway_annotationSummary.pkl'
    with open(inputFile, 'rb') as pickledFile:
        summaryAnnotation = pickle.load(pickledFile)

    annLineCount = int(sp.getoutput('wc -l %s' %(annFile)).split()[0])

    for expFile in expFileList:
        with open(expFile, 'rb') as pickledFile:
            expression = pickle.load(pickledFile)

        expAccession = expFile.split('/')[-1].split('_')[-1].split('.')[-2]

        clusters = summaryAnnotation # we only have clusters for these samples, no need to sort the clusters
        clusterList = list(clusters.keys())
        clusterCount = len(clusterList)

        clusterListIndMap = {}
        for i in range(len(clusterList)):
            clusterListIndMap[clusterList[i]] = i


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
        #while cgi < 2:
        #while cgi < 18431:# >>>>
            
            #print('cgi')
            #print(cgi)

            'we are in the gene territory, we are walking on the genes'
   
            geneID = geneIDList[cgi]
            gene_chr = geneList[geneIDList[cgi]].chrom
            #print(geneID) # >>>> test
            
            gene_strand = geneList[geneID].strand
            gene_start = geneList[geneID].start
            gene_end = geneList[geneID].end
            gene_length = gene_end - gene_start
            gene_length_unit = int(gene_length/100)
            if gene_length_unit == 0:
                #print('gene smaller than 100bp, skipping it') # >>>> test
                cgi = cgi + 1
                continue

            gene_length_last_unit = gene_length - (99* gene_length_unit)
            # TODO: ideally: something to fix: the gene_length_last_unit for negative strand versus positive strand

            extension_start = gene_start - extension
            extension_end = gene_end + extension

            #print('%ss, %s, %s' %(gene_chr, extension_start, extension_end)) # >>>> test

            
            ''' picking the label/class matrix based on the gene expression level'''
            #TODO catch exception for when geneID is not in expression
            gene_exp = expression[geneID]
            if gene_exp == 0:
                expMatInd = 0
            elif gene_exp > 1:
                expMatInd = 2
            else:
                expMatInd = 1

            #print(gene_exp)
            #print('expMatInd')
            #print(expMatInd) # >>>>> test

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


            ''' 
            if this gene starts somewhere in the preivous genomic region, 
            I will go back in the annotation file to the beginning of annotation for the previous gene
            changed this th the next while loop - basically we will go back until the annotations start before the gene territory

            '''

            while (ann_start > extension_start) and (gene_chr == ann_chr): # in case of overlapping genes
                #print('ann start greater than extension start, getting back in annotation until it is not') # >>>> test
                #print(annLineInd) # >>>> test
                #print('%ss, %s, %s' %(ann_chr, ann_start, ann_end)) # >>>> test
                
                annLineInd = max(annLineInd - 5, 1)
                line = linecache.getline(annFile, annLineInd)
                annLineInd +=1
                ann_line_count += 1
                fields = line.strip().split()

                ann_chr = fields[0]
                ann_start = int(fields[1])
                ann_end = int(fields[2])
                
                #print('overlapping genes here') # >>>> test
                #print(annLineInd) # >>>> test

            #while ((ann_start < extension_end) or not(gene_chr == ann_chr)) and geneMatWalkIndex < 160: 
            while ((ann_start < extension_end) and (gene_chr == ann_chr)) and geneMatWalkIndex < 160:
                
#            if (ann_chr != previous_ann_chr): # in case of chromosome change because annotation moved to the next chromosome
 #           if (ann_chr != gene_chr): # if annotation moved to the next chromosome, but gene has not yet moved to the next chromosome

                '''
                NOTE: the second condition is for when we are at the end of a gene's territory, and then at the end of a chromosome, so when we go to the next gene. At some point I changed the second condition from "or" to "and" 
                ann_start is not smaller than extension_end, and we need to read the annotations until we are on the same chromosome
                The condition to be in one gene's territory (and not the next one), and it is for reading the annotations
                while in THIS gene territory, read annotations until we reach to the next genomic region (plus extension)
      
                '''

                #print('in the gene region, reading annotations until we are out') # >>>> test
                #print('%ss, %s, %s' %(gene_chr, extension_start, extension_end)) # >>>> test

                # in the next sections we are processing the annotation
                if ann_chr == gene_chr: # if we are in the same choromosome
                    if ((ann_start < extension_end and ann_start > extension_start) or 
                        (ann_end < extension_end and ann_end > extension_start) or
                        (ann_start < extension_start and ann_end > extension_end)):


                        ''' We are in the genonimc region (with extension)'''
                        #print('annotation') # >>>> test
                        #print('%ss, %s, %s' %(ann_chr, ann_start, ann_end)) # >>>> test


                        ''' Taking off for filling the matrices... '''
                        
                        adjusted_ann_start = max(0, ann_start - extension_start)
                        adjusted_ann_end = min(ann_end - extension_start, extension_end - extension_start)
                        adjusted_ann_length = adjusted_ann_end - adjusted_ann_start

                        ann_cluster = fields[3]
                        clusterInd = clusterListIndMap[ann_cluster]

                        #ann_class = ann_cluster.split('_')[1]
                        #classInd = classListIndMap[ann_class]

                        # expMatInd

                        if gene_strand == '+':  # munching 100bp s from annotation, filling the geneMat
                            while(adjusted_ann_length > 0) and geneMatWalkIndex < 30:
                                if (adjusted_ann_length >= 100*(1- previous_fill)):
                                    #classMats[expMatInd][classInd][geneMatWalkIndex] += 1 - previous_fill
                                    clusterMats[expMatInd][clusterInd][geneMatWalkIndex] += 1 - previous_fill
                                    adjusted_ann_length -= 100*(1- previous_fill)
                                    previous_fill = 0
                                    geneMatWalkIndex +=1
                                else:
                                    #classMats[expMatInd][classInd][geneMatWalkIndex] += adjusted_ann_length/100
                                    clusterMats[expMatInd][clusterInd][geneMatWalkIndex] += adjusted_ann_length/100
                                    previous_fill += adjusted_ann_length/100
                                    adjusted_ann_length = 0
                                    if previous_fill > 1:
                                        previous_fill = 1

                            while(adjusted_ann_length > 0) and geneMatWalkIndex < 129:
                                if (adjusted_ann_length >= gene_length_unit*(1- previous_fill)):
                                    #classMats[expMatInd][classInd][geneMatWalkIndex] += 1 - previous_fill
                                    clusterMats[expMatInd][clusterInd][geneMatWalkIndex] += 1 - previous_fill
                                    adjusted_ann_length -= gene_length_unit*(1 - previous_fill)
                                    previous_fill = 0
                                    geneMatWalkIndex +=1
                                else:
                                    #classMats[expMatInd][classInd][geneMatWalkIndex] += adjusted_ann_length/gene_length_unit
                                    clusterMats[expMatInd][clusterInd][geneMatWalkIndex] += adjusted_ann_length/gene_length_unit
                                    previous_fill += adjusted_ann_length/gene_length_unit
                                    adjusted_ann_length = 0
                                    if previous_fill > 1:
                                        previous_fill = 1

                            if (adjusted_ann_length > 0) and geneMatWalkIndex == 129:
                                if (adjusted_ann_length >= gene_length_last_unit*(1- previous_fill)):
                                    #classMats[expMatInd][classInd][geneMatWalkIndex] += 1 - previous_fill
                                    clusterMats[expMatInd][clusterInd][geneMatWalkIndex] += 1 - previous_fill
                                    adjusted_ann_length -= gene_length_last_unit*(1 - previous_fill)
                                    previous_fill = 0
                                    geneMatWalkIndex +=1
                                else:
                                    #classMats[expMatInd][classInd][geneMatWalkIndex] += adjusted_ann_length/gene_length_last_unit
                                    clusterMats[expMatInd][clusterInd][geneMatWalkIndex] += adjusted_ann_length/gene_length_last_unit
                                    previous_fill += adjusted_ann_length/gene_length_last_unit
                                    adjusted_ann_length = 0
                                    if previous_fill > 1:
                                        previous_fill = 1

                            while(adjusted_ann_length > 0) and (geneMatWalkIndex < 160):
                                if (adjusted_ann_length >= 100*(1- previous_fill)):
                                    #classMats[expMatInd][classInd][geneMatWalkIndex] += 1 - previous_fill
                                    clusterMats[expMatInd][clusterInd][geneMatWalkIndex] += 1 - previous_fill
                                    adjusted_ann_length -= 100*(1- previous_fill)
                                    previous_fill = 0
                                    geneMatWalkIndex +=1
                                else:
                                    #classMats[expMatInd][classInd][geneMatWalkIndex] += adjusted_ann_length/100
                                    clusterMats[expMatInd][clusterInd][geneMatWalkIndex] += adjusted_ann_length/100
                                    previous_fill += adjusted_ann_length/100
                                    adjusted_ann_length = 0
                                    if previous_fill > 1:
                                        previous_fill = 1


                        if gene_strand == '-':
                           # print('here in the negative strand') # >>>> test
                            while(adjusted_ann_length > 0) and geneMatWalkIndex < 30:
                            #    print('first while') # >>>> test
                                if (adjusted_ann_length >= 100*(1- previous_fill)):
                                    #classMats[expMatInd][classInd][159 - geneMatWalkIndex] += 1 - previous_fill
                                    clusterMats[expMatInd][clusterInd][159 - geneMatWalkIndex] += 1 - previous_fill
                                    adjusted_ann_length -= 100*(1- previous_fill)
                                    previous_fill = 0
                                    geneMatWalkIndex +=1
                                else:
                                    #classMats[expMatInd][classInd][159 - geneMatWalkIndex] += adjusted_ann_length/100
                                    clusterMats[expMatInd][clusterInd][159 - geneMatWalkIndex] += adjusted_ann_length/100
                                    previous_fill += adjusted_ann_length/100
                                    adjusted_ann_length = 0
                                    if previous_fill > 1:
                                        previous_fill = 1

                            if (adjusted_ann_length > 0) and geneMatWalkIndex == 30:
                                if (adjusted_ann_length >= gene_length_last_unit*(1- previous_fill)):
                                    #classMats[expMatInd][classInd][159 - geneMatWalkIndex] += 1 - previous_fill
                                    clusterMats[expMatInd][clusterInd][159 - geneMatWalkIndex] += 1 - previous_fill
                                    adjusted_ann_length -= gene_length_last_unit*(1 - previous_fill)
                                    previous_fill = 0
                                    geneMatWalkIndex +=1
                                else:
                                    #classMats[expMatInd][classInd][159 - geneMatWalkIndex] += adjusted_ann_length/gene_length_last_unit
                                    clusterMats[expMatInd][clusterInd][159 - geneMatWalkIndex] += adjusted_ann_length/gene_length_last_unit
                                    previous_fill += adjusted_ann_length/gene_length_last_unit
                                    adjusted_ann_length = 0
                                    if previous_fill > 1:
                                        previous_fill = 1


                            while(adjusted_ann_length > 0) and geneMatWalkIndex < 129:
                                if (adjusted_ann_length >= gene_length_unit*(1- previous_fill)):
                                    #classMats[expMatInd][classInd][159 - geneMatWalkIndex] += 1 - previous_fill
                                    clusterMats[expMatInd][clusterInd][159 - geneMatWalkIndex] += 1 - previous_fill
                                    adjusted_ann_length -= gene_length_unit*(1 - previous_fill)
                                    previous_fill = 0
                                    geneMatWalkIndex +=1
                                else:
                                    #classMats[expMatInd][classInd][159 - geneMatWalkIndex] += adjusted_ann_length/gene_length_unit
                                    clusterMats[expMatInd][clusterInd][159 - geneMatWalkIndex] += adjusted_ann_length/gene_length_unit
                                    previous_fill += adjusted_ann_length/gene_length_unit
                                    adjusted_ann_length = 0
                                    if previous_fill > 1:
                                        previous_fill = 1


                            while(adjusted_ann_length > 0) and (geneMatWalkIndex < 160):
                                if (adjusted_ann_length >= 100*(1- previous_fill)):
                                    #classMats[expMatInd][classInd][159 - geneMatWalkIndex] += 1 - previous_fill
                                    clusterMats[expMatInd][clusterInd][159 - geneMatWalkIndex] += 1 - previous_fill
                                    adjusted_ann_length -= 100*(1- previous_fill)
                                    previous_fill = 0
                                    geneMatWalkIndex +=1
                                else:
                                    #classMats[expMatInd][classInd][159 - geneMatWalkIndex] += adjusted_ann_length/100
                                    clusterMats[expMatInd][clusterInd][159 - geneMatWalkIndex] += adjusted_ann_length/100
                                    previous_fill += adjusted_ann_length/100
                                    adjusted_ann_length = 0
                                    if previous_fill > 1:
                                        previous_fill = 1


                                # if gene strand is positive, or negative, flag which cluster list we use
                                # the only difference between the two strands is that I am going to reverse the index of

                if geneMatWalkIndex < 160: # we read annotations until we cover the genee
                    #print('geneMatWalkInd') # >>>> test
                    #print(geneMatWalkIndex) # >>>> test
                    #print('annLineInd') # >>>> test
                    #print(annLineInd) # >>>> test
                    line = linecache.getline(annFile, annLineInd)
                    annLineInd +=1
                    ann_line_count += 1
                    fields = line.strip().split()

                    ann_chr = fields[0]
                    ann_start = int(fields[1])
                    ann_end = int(fields[2])
                    #print('%ss, %s, %s' %(ann_chr, ann_start, ann_end)) # >>>> test

                #if geneMatWalkIndex >= 159:
                  #  print('gene Done')
                  #  print(geneMatWalkIndex)


            cgi += 1 # next gene
            previous_extension_end = extension_end
            previous_gene_chr = gene_chr
            

        linecache.clearcache()

        clusterOutput = {'clusterMats': clusterMats, 'clusterList': clusterListIndMap}
        outputFile = sampleFolderAdd + expAccession + '_expSummary.pkl'
        with open(outputFile, 'wb') as f:
            pickle.dump(clusterMats, f)

        
        book = pd.DataFrame(clusterMats[0])
        sns.heatmap(book)
        plt.show()
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



# get transcript data from previous segway runs (mapping of the Segway accession to the tissue)
##################################################


#not sure if this helps
bioinfo = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/testBatch105/fromAPI/biosample_tissue_info.pkl'
with open(bioinfo, 'rb') as inputFile:
    bio_list = pickle.load(inputFile)

# I have the transcript data for Segway runs and that includes almost all Mehdi samples except for one. So I will fetch the transcript data from these samples - nevermind. since most these are cell lines they don't specify the donor and did not show up in my segway search (since for that I used the donor)
tissueList = []
segAccession_tissue_map = {}
for key in list(bio_list.keys()):
    tissue = bio_list[key][0]
    tissueList.append(tissue)
    segAccession_tissue_map[tissue] = key

print(segAccession_tissue_map['MCF-7'])

for item in sampleList:

    # skipping the concat runs
    if 'concat' in item: 
        continue

    # where to get the transcript data
    tissue = item.split('_')[0]
    if tissue == 'CD14-positive': # we don't have transcript for this one at the moment
        continue 
    print(tissue)
    
    if tissue in tissueList:
        # segway folder
        segAccession = segAccession_tissue_map[tissue]
        RNA_folder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/testBatch105/fromAPI/' + segAccession + '/'
        RNA_file = glob.glob(RNA_folder + 'geneExp_dict_*.pkl')

    # if there is no transcript data
    if len(RNA_file) == 0:
        continue

    print(RNA_file)
    # now to the transcript comparison:

    
