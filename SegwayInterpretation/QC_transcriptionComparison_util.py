# the util file made from QC_transcriptionComparison.py


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Annotation and classes (biolabels)

import pickle

# keeping the unique cluster_class labels and their info for annotations
# I defined my own classes since namedtuples are immutable - I preferred namedtuples to dict since they are defined with fields - 

class Annotation(object):
    def __init__(self, called, biolabel, cluster, color, bp_count, region_count, region_dist):
        self.called = called
        self.biolabel = biolabel
        self.cluster = cluster
        self.color = color
        self.bp_count = bp_count
        self.region_count = region_count
        self.region_dist = region_dist # not used currently

    def __str__(self):
        return 'called = %s, biolabel = %s, cluster = %s, color = %s, bp_count = %d, region_count = %d, region_dist = cluster.region_dist' %(self.called, self.biolabel, self.cluster, self.color, self.bp_count, self.region_count)

class AnnotationClass(object):
    def __init__(self, biolabel, clusters, color, bp_count, region_count, region_dist):
        self.biolabel = biolabel
        self.clusters = clusters
        self.color = color
        self.bp_count = bp_count
        self.region_count = region_count # the merged count
        self.region_dist = region_dist # not used currently

    def __str__(self):
        return 'biolabel = %s, clusters = %s, color = %s, bp_count = %d, region_count = %d, region_dist = none' %(self.biolabel, self.clusters, self.color, self.bp_count, self.region_count)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> EpigenomeAnnotationMeta

# this is the information about a epigenome

class EpigenomeAnnotation(object):
    def __init__(self, accession, tissue, Uberon, donor_age, donor_sex, transcriptID):
        self.accession = accession
        self.tissue = tissue
        self.Uberon = Uberon
        self.donor_age = donor_age
        self.donor_sex = donor_sex
        self.transcriptID = transcriptID

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> AnnotationMeta

class SegwayAnnotationMeta(object):
    def __init__(self, accession, tissueInfo, donorID, donorAge, donorSex, bedFile, transcriptFile, ChromFile, reside):
        self.accession = accession # accession of the annotation
        self.tissueInfo = tissueInfo # tissue info
        self.donorID = donorID # donor info
        self.donorAge = donorAge # donor info
        self.donorSex = donorSex # donor info
        self.bedFile = bedFile # name of the bed file in the folder
        self.transcriptFile = transcriptFile # 'none'/name of the transcript file in the folder
        self.ChromFile = ChromFile # 'none'/name of the chromhmm file in the folder
        self.reside = reside # local address of the folder
        self.batch = batch # which run was it

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Genes and Exons
class Gene(object):
    def __init__(self, name, ENS_ID, gtype, chrom, start, end, strand):
        self.name = name 
        self.ENS_ID = ENS_ID # ensemble ID
        self.gtype = gtype
        self.chrom = chrom
        self.start = start # always smaller than end, direction is defined by strand
        self.end = end
        self.length = end - start
        self.strand = strand
        self.exon_list = [] # list of exons
        self.exon_count = len(self.exon_list)
        self.intron_list = [] # list of introns
####        self._flag

    def __str__(self):
        return 'name = %s, ENS_ID = %s, gtype = %s, chrom = %s, start = %d, end = %d, strand = %s, exon_count = %d' %(self.name, self.ENS_ID, self.gtype, self.chrom, self.start, self.end, self.strand, self.exon_count)

    def printExonList(self):
        for i in range(len(self.exon_list)):
            print(self.exon_list[i])

    def getIntronList(self):
        if len(self.exon_list) > 1:
            intronList = []
            for i in range(len(self.exon_list)-1):
                in_start = self.exon_list[i].end + 1
                in_end = self.exon_list[i+1].start - 1
                intron = Intron(in_start, in_end)
                intronList.append(intron)

            self.intron_list = intronList
        else:
            print('only one exon recorded for the gene')

            
    def printIntronList(self):
        for i in range(len(self.intron_list)):
            print(self.intron_list[i])

            
    def sortExonList(self):
        
        '''sorts exon list'''
        
        exon_sorted_list = []
        exon_start_list = []
        for i in range(len(self.exon_list)):
            exon_start_list.append(self.exon_list[i].start)

        exon_sorted_list = []
        inds = sorted(range(len(exon_start_list)), key = lambda k: exon_start_list[k])
        for i in range(len(inds)):
            exon_sorted_list.append(self.exon_list[inds[i]])

        self.exon_list = exon_sorted_list

            
    def exonOverlaps(self):

        '''reports overlapping exons'''
        #TODO: define and check the exon sort flag
        # to report if there are overlapping exons
        exonOverlapIncidents = 0
        for i in range(len(self.exon_list)-1):
            if self.exon_list[i].end > self.exon_list[i+1].start:
                exonOverlapIncidents +=1

        return(exonOverlapIncidents)
    

class Exon(object):
    def __init__(self, ENS_ID, start, end):
        self.ENS_ID = ENS_ID
        self.start = start
        self.end = end

    def __str__(self):
        return 'ENS_ID = %s, start = %d, end = %d' %(self.ENS_ID, self.start, self.end)

class Intron(object): # introns don't have ID, but just identified by exon_end+1 and exon_start-1
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __str__(self):
        return 'start = %d, end = %d' %(self.start, self.end)

#########################################
# transcriptomic files preprocessing
#########################################

def transcriptFile_preprocessing(folder, fileName):

    expFile = folder + fileName
    print(expFile)
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
    expAccession = re.split('_|\.', fileName)[2]
    outputFile = folder + 'geneExp_dict_' + expAccession + '.pkl'
    print('printing this file %s' %(outputFile))
         
    with open(outputFile, 'wb') as f:
        pickle.dump(expression, f)



#########################################
# annotation preprocessing
#########################################



#########################################
# Meta info from the annotations
#########################################

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# this is for Segway

def annotation_generalInfo_clusters(bedFileAdd):

    splitFileName = bedFileAdd.split('/')
    bedFileName = splitFileName[-1].split('.')[0]
    index = bedFileAdd.index(bedFileName + '.bed')
    outputFolder = bedFileAdd[0:index]

    previous_class = ''

    clusters = {}
    classes = {}
    c = 0
    with open(bedFileAdd, 'r') as annotations:

        # annotations have no header
        # header = annotations.readline()

#        f = open(bedFileAdd, 'r')
#        line = f.readline()
        
        for line in annotations:
            c += 1
            fields = line.strip().split()

            # doing the clusters first
            if fields[3] in clusters.keys():

                clusters[fields[3]].bp_count += int(fields[2]) - int(fields[1])
                clusters[fields[3]].region_count += 1
                
            else:
                
                called = fields[3]
                biolabel = fields[3].split('_')[1]
                cluster = fields[3].split('_')[0]
                color = fields[8]
                bp_count = int(fields[2]) - int(fields[1])
                region_count = 1
                region_dist = 1
                clusters[fields[3]] = Annotation(called, biolabel, cluster, color, bp_count, region_count, region_dist)


            # doing the class
            if previous_class == fields[3].split('_')[1]:
                classes[previous_class].bp_count +=  int(fields[2]) - int(fields[1])
            else:
                current_class = fields[3].split('_')[1]
                if current_class in classes.keys():
                    classes[current_class].bp_count +=  int(fields[2]) - int(fields[1])
                    classes[current_class].region_count += 1
                else:
                    clusterList = [] # not filling it now, it can be filled later using annotations 
                    biolabel = current_class
                    color = fields[8]
                    bp_count = int(fields[2]) - int(fields[1])
                    region_count = 1
                    region_dist = 1
                    classes[biolabel] = AnnotationClass(biolabel, clusterList, color, bp_count, region_count, region_dist)

            previous_class = fields[3].split('_')[1]

            
    # filling up the cluster distribution from distribution file
    
    segmentSizesFile = outputFolder + 'segment_sizes.tab.txt'
    with open(segmentSizesFile, 'r') as inputFile:
        lines = inputFile.readlines()
    
    clusterList = list(clusters.keys())
    for cluster in clusterList:
        cluster_number = int(cluster.split('_')[0])
        fields = lines[cluster_number + 2].strip().split('\t')
        info = {'num.segs': int(fields[1]), 'mean.len': float(fields[2]), 'median.len': float(fields[3]), 'stdev.len': float(fields[4]), 'num.bp': float(fields[5]), 'frac.bp': float(fields[6])}
        clusters[cluster].region_dist = info

    annotationSummary = {"classes": classes, "clusters": clusters}
    outputFile = outputFolder + bedFileName + '_annotationSummary.pkl'
    with open(outputFile, 'wb') as f:
        pickle.dump(annotationSummary, f)

    print('annotation summary saved in %s' %(outputFile))

########################################################
# this is for Segway
def annotation_generalInfo_clusters_noclass(bedFileAdd):

    splitFileName = bedFileAdd.split('/')
    bedFileName = splitFileName[-1].split('.')[0]
    index = bedFileAdd.index(bedFileName + '.bed')
    outputFolder = bedFileAdd[0:index]


    previous_class = ''

    clusters = {}
    c = 0

    # get count of header lines:
    print('header fix ...')
    headerCount = 0
    with open(bedFileAdd, 'r') as annotations:
        line = annotations.readline()
        fields = line.strip().split()
        while fields[0] != 'chr1':
            headerCount += 1
            line = annotations.readline()
            fields = line.strip().split()

    with open(bedFileAdd, 'r') as annotations:

        print('getting the annotation summary...')

        # annotations have no header
        # header = annotations.readline()

#        f = open(bedFileAdd, 'r')
#        line = f.readline()

        # removing headers:

        if headerCount > 0:
            for i in range(headerCount):
                line = annotations.readline()

        print('starting the annotation ...')
        for line in annotations:
            
            c += 1
            fields = line.strip().split()

            # doing the clusters first
            if fields[3] in clusters.keys():

                clusters[fields[3]].bp_count += int(fields[2]) - int(fields[1])
                clusters[fields[3]].region_count += 1
                
            else:
                
                called = fields[3]
                biolabel = 'none'
                cluster = fields[3]
                color = fields[8]
                bp_count = int(fields[2]) - int(fields[1])
                region_count = 1
                region_dist = 1
                clusters[fields[3]] = Annotation(called, biolabel, cluster, color, bp_count, region_count, region_dist)

    annotationSummary = clusters
    outputFile = outputFolder + bedFileName + '_annotationSummary.pkl'
    with open(outputFile, 'wb') as f:
        pickle.dump(annotationSummary, f)

    print('annotation summary saved in %s' %(outputFile))


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# this is for chmm
# chmm has no cluster, we only have class data
 
import pickle

def annotation_generalInfo_classes_chmm(bedFileAdd):

    splitFileName = bedFileAdd.split('/')
    bedFileName = splitFileName[-1].split('.')[0]
    index = bedFileAdd.index(bedFileName)
    outputFolder = bedFileAdd[0:index]

    previous_class = ''

    classes = {}
    c = 0
    with open(bedFileAdd, 'r') as annotations:

        # annotations have no header
        # header = annotations.readline()

        f = open(bedFileAdd, 'r')
        line = f.readline()
        
        for line in annotations:
            c += 1
            fields = line.strip().split()

            # doing the class 
            if previous_class == fields[3]:
                classes[previous_class].bp_count +=  int(fields[2]) - int(fields[1])
            else:
                current_class = fields[3]
                if current_class in classes.keys():
                    classes[current_class].bp_count +=  int(fields[2]) - int(fields[1])
                    classes[current_class].region_count += 1
                else:
                    clusterList = [] # not filling it now, it can be filled later using annotations 
                    biolabel = current_class
                    color = fields[8]
                    bp_count = int(fields[2]) - int(fields[1])
                    region_count = 1
                    region_dist = 1
                    classes[biolabel] = AnnotationClass(biolabel, clusterList, color, bp_count, region_count, region_dist)

            previous_class = fields[3]

    annotationSummary = {"classes": classes}
    outputFile = outputFolder + bedFileName + '_annotationSummary.pkl'
    with open(outputFile, 'wb') as f:
        pickle.dump(annotationSummary, f)

    print('annotation summary saved in %s' %(outputFile))


# function for the ccre
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 
import pickle

def annotation_generalInfo_classes_ccre(bedFileAdd):

    splitFileName = bedFileAdd.split('/')
    bedFileName = splitFileName[-1].split('.')[0]
    index = bedFileAdd.index(bedFileName)
    outputFolder = bedFileAdd[0:index]

    previous_class = ''

    classes = {}
    c = 0
    with open(bedFileAdd, 'r') as annotations:

        # annotations have no header
        # header = annotations.readline()

        f = open(bedFileAdd, 'r') # >>>> test
        line = f.readline() # >>>> test
        
        for line in annotations:
            c += 1
            fields = line.strip().split()

            # doing the class 
            current_class = fields[9]
            if current_class in classes.keys():
                classes[current_class].bp_count +=  int(fields[2]) - int(fields[1])
                classes[current_class].region_count += 1
            else:
                clusterList = [] # doesn't apply
                biolabel = current_class
                color = fields[8]
                bp_count = int(fields[2]) - int(fields[1])
                region_count = 1
                region_dist = 1
                classes[biolabel] = AnnotationClass(biolabel, clusterList, color, bp_count, region_count, region_dist)

    annotationSummary = {"classes": classes}
    outputFile = outputFolder + bedFileName + '_annotationSummary.pkl'
    with open(outputFile, 'wb') as f:
        pickle.dump(annotationSummary, f)

    print('annotation summary saved in %s' %(outputFile))

