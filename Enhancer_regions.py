# the goal of this code is to record the enhancer regions within chromosome 19.
### sections ###
#
# 0. Initials
# 1. Get chr19 files
# 2. First region analyses 
# 2.1 Getting the regions
# 2.2 Exploratory analyses for the region

#########################################
# 0. Initials
#########################################

import os
import shutil
import subprocess

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
dataSubFolder = 'testBatch105/fromAPI/'

inputFile = dataFolder + dataSubFolder + 'metaInfo.pkl'
with open(inputFile, 'rb') as f:
    annMeta = pickle.load(f)

annAccessionList = list(annMeta.keys())
annAccession = annAccessionList[104]
print(annAccession)

annAccession = 'ENCSR566HBT'
index = accession_index_map[annAccession]
print(index)

#########################################
# 1. Get chr19 files
#########################################

#coverage = np.zeros((max_region_count, int(region_length/10)+40))
# create the en file
for annAccession in annAccessionList:

    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
    print(sampleFolderAdd)

    index = accession_index_map[annAccession]
    print(index)

    # get the .bed file chr19 extract. 
    originalBedFile = annMeta[annAccession]['bedFile']
    originalBed_accession = originalBedFile.split('.')[0]
    segwayFile = dataFolder + dataSubFolder + annAccession + '/' + originalBed_accession + '_filteredSorted.bed.gz'

    # load the mnomics, and extract the enhancer labels for mnemonics
    label_term_mapping = {}
    mnemonics_file = sampleFolderAdd + 'mnemonics_v02.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term

    enLabels = []
    for label in label_term_mapping.keys():
        if 'Enhancer' in label_term_mapping[label]:
            enLabels.append(label)

    one_digit_label = []
    two_digit_label = []
    for label in enLabels:
        if len(label) == 2:
            two_digit_label.append(label)
        else:
            one_digit_label.append(label)

    label_string = ''
    if len(one_digit_label)>0:
        label_string = label_string + '['
        for label in one_digit_label:
            label_string = label_string + label

        label_string = label_string + ']'
        if len(two_digit_label)>0:
            label_string = label_string + '|'


    if len(two_digit_label)>0:
        for label in two_digit_label:
            label_string = label_string + '(' + label + ')' + '|'
        label_string = label_string[0:-1]

    
    os.system('gunzip %s' %(segwayFile))

    command = "grep -E 'chr19.*\\t(%s)_\' %s" %(label_string, segwayFile[0:-3])
    print(command)
    print(counter)
    out = sampleFolderAdd + 'chr19_enh.bed'
    f = open(out, 'w')
    subprocess.run(command, shell=True, stdout=f)

    os.system('gzip %s' %(segwayFile[0:-3]))

    print('chr19 enhancer regions selected')

    # if it is the first file:


#########################################
# 2. First region analyses 
#########################################


# 2.1 Getting the regions
#########################################

'''
pseudo code:

region_list = []
For each sample:
  For each enhancer_label:
     min_dis = min(min(distance(enahncer_label, (for region in region_list))), 1000)
       if min_dis = 1000:
          create new enhancer_region, add it to region_list
       else:
          add the enhancer_label to the closest region from region_list
          (v0: adjust the center, v1: not adjust the center)

'''


region_length = 12 # x10 is the length of the enhancer region
max_region_count = 40000 # the maximum count of regions, not sure if is enough for files, the estimate is to be 40,000 (60mil/3bil * 2mil) [chr19/whole genome/enhancer count]
region_coverage = np.zeros((max_region_count, region_length+12)) # I keep the coverage from 1000 each side from original center
sample_belong = np.zeros((max_region_count, 104))
# keep a center for the enhancer region, if we are far away from it (more than half, create another region)
counter = 1
# now process the regions:
# we keep a region center, and we work based on that.
centers = np.zeros((max_region_count))
sums = np.zeros((max_region_count))
region_start = np.zeros((max_region_count))
segment_counts = np.zeros((max_region_count))
cc_dist = 1000

# do the analysis
for annAccession in annAccessionList:

    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
    print(sampleFolderAdd)

    index = accession_index_map[annAccession]
    print(index)

    out = sampleFolderAdd + 'chr19_enh.bed'

    if index == 0:

        with open(out, 'r') as enFile:
            line = enFile.readline()
            fields = line.split()
            center = (int(fields[2]) - int(fields[1]))/2 + int(fields[1])
            region_start[0] = center - 800
            

    # centers # array, list of centers
    # sums # sum of segments in regions
    # segments_count # count of the segments in region, array
    # coverage #  matrix of coverage
    # counter #  the counter/walker on regions
    # region_start #  an array keep start of regions

    sample_line_count = 0
    with open(out, 'r') as enFile:
        #line = enFile.readline()
        #print(line)
        

        #fields = line.split()
        #center = int(fie)
        for line in enFile:
            sample_line_count +=1
            fields = line.split()
            seg_start = int(fields[1])
            seg_end = int(fields[2])
            seg_unit_length = int((seg_end - seg_start)/100)
            center = (seg_end - seg_start)/2 + seg_start

            center_list = centers[0:counter]
            diff = np.abs(center_list - center)

            minInd = np.argmin(diff)

            min_dist = np.min([diff[minInd], cc_dist])

            if min_dist == cc_dist: # new segment needs to be added
                sums[counter] = center
                segment_counts[counter] = 1
                centers[counter] = center
                region_start[counter] = center - 1200 
                cov_start = int((seg_start - region_start[counter])/100)
                cov_end = int(np.min([cov_start+seg_unit_length, 23]))
                region_coverage[counter, cov_start:cov_end] += 1
                sample_belong[counter, index] += 1
                counter +=1

            else: # add it to the already existing region, update region info
                sums[minInd] += center
                segment_counts[minInd] += 1
                centers[minInd] = int(sums[minInd]/segment_counts[minInd])
                cov_start = int((seg_start - region_start[minInd])/100)
                cov_end = int(np.min([cov_start+seg_unit_length, 23]))
                region_coverage[minInd, cov_start:cov_end] += 1
                sample_belong[minInd, index] +=1


            #if (np.max(region_coverage))>(index +1):
                #print(sample_line_count)
                #break
        
        print('line_count %d' %(sample_line_count))
        print('counter %d' %(counter))
        print('______________')

region_coverage = region_coverage[1:counter,]
region_start = region_start[1:counter]
centers = centers[1:counter]
sample_belong = sample_belong[1:counter,]

# 2.2 Exploratory analyses for the region
#########################################

# >>>> heatmap of covreage
sns.heatmap(region_coverage[1:100, ])
plt.show()

# >>>> saving the v0 output ccdist = 1200, coverage = 160
select_ind = np.arange(5,156, 10)
sum_rg = region_coverage[:, select_ind]
result = {}
result['region_coverage'] = sum_rg
result['region_start'] = region_start
result['centers'] = centers
result['sample_belong'] = sample_belong

file = dataFolder + 'enhancer_regions_v0.pkl'
with open(file, 'wb') as f:
    pickle.dump(result, f)
# <<< saving the v0 output

# >>>> saving the v1 output ccdist = 1000, covreage = 24 (240)
result = {}
result['region_coverage'] = region_coverage
result['region_start'] = region_start
result['centers'] = centers
result['sample_belong'] = sample_belong

file = dataFolder + 'enhancer_regions_v1.pkl'
with open(file, 'wb') as f:
    pickle.dump(result, f)
# <<< saving the v0 output


# >>>> distance of centers from original region start (original center)
centers_dis = centers - (region_start+1200) # this is for v1
print(sum(abs(centers_dis) >1000))
print(sum(abs(centers_dis) <1000)/counter)
plt.hist(centers_dis[1:])
plt.show()
# regions can overalp, but centers don't  for 99% of centers
# get the difference of the cents
book = np.sort(centers)
diff = np.zeros(len(book))
for i in range(len(diff)-1):
    diff[i] = book[i+1] - book[i]

maxInd = (np.where(diff == max(diff))) # there is 2mil gap
diff[diff > 5000] = 5000
plt.hist(diff)
plt.show()

# for the v1, 99% of regions would still include the initial center

# pick regions where center is within the region.
len(centers_dis)
centers_in_regions = abs(centers_dis) < 600
select_ind = np.arange(5,156, 10)
selected_regions = region_coverage[centers_in_regions,]
selected_regions = selected_regions[:,select_ind]
sns.heatmap(selected_regions[0:40,])
plt.show()

# >>>> count of regions covered with count of cell types - this changes between v0 and v1
max_cov = np.max(region_coverage, axis = 1)
plt.hist(max_cov)
plt.show()
print(sum(max_cov > 11))

    
# >>>> how many of regions have "peaks"
mean_cov = np.sum(region_coverage, axis = 1)/24
peaks= np.divide(max_cov, mean_cov)
plt.hist(peaks)
plt.show()

print(sum(peaks<2))
temp = (peaks<2) # WIDE peaks
temp_cov = region_coverage[temp,]
sns.heatmap(temp_cov[0:100,])
plt.show()
# IMPORTANT observation: regions with high peaks (20+ can be VERY wide, up to 2k, so the high/mean ratio can be very low <2, but we still observe a peak of 30+, in fact it seems to be more common like this)

max_temp_cov = np.max(temp_cov, axis = 1)
sorted_max_temp_cov = np.argsort(max_temp_cov)

sorted_temp_cov = temp_cov[sorted_max_temp_cov,]
sns.heatmap(sorted_temp_cov)
plt.show()

# now with the high peaks to mean ratio
print(sum(peaks>5)) # are these the slim peaks?
temp = peaks >= 5  # SHARP peaks
temp_cov = region_coverage[temp,]
sns.heatmap(temp_cov[0:100,])
plt.show()

p1 = peaks >= 2
p2 = peaks < 5
p3 = np.logical_and(p1, p2)
print(sum(p3))

temp_cov = region_coverage[p3,] 
sns.heatmap(temp_cov[0:100,])
plt.show() # MEDUM PEAKS 

sorted_peaks = np.argsort(peaks)

cum_peaks = np.zeros(len(sorted_peaks)-1)
cum_peaks[0] = sorted_peaks[0]
for i in range(len(cum_peaks)-1):
    cum_peaks[i+1] = cum_peaks[i] + sorted_peaks[i+1]

plt.plot(sorted_peaks[0:100])
plt.show()


# >>>> count of regions with samples - this does not change between v0 and v1
sample_belong_binary = sample_belong > 0 
sample_count = np.sum(sample_belong_binary, axis = 1)
plt.hist(sample_count)
plt.show()

# sort the start of regions, and get their distance 
book = np.sort(region_start)
diff = np.zeros(len(book))
for i in range(len(diff)-1):
    diff[i] = book[i+1] - book[i]

print(sum(diff > 5000)) # we have a few gaps like this, we will replace them with 5000 (54 of them)

maxInd = (np.where(diff == max(diff))) # there is 2mil gap
diff[diff > 5000] = 5000
plt.hist(diff)
plt.show()

maxCov_sortedInd = argsort()

sns.heatmap(region_coverage[1:300, :])
plt.show()

sns.heatmap(sample_belong[1:400],)
plt.show()

sima = kado[sib,]

sns.heatmap(kado)
plt.show()

method = 'ward'
sib = sns.clustermap(sima[0:20000,], method = method)
plt.show()

dendro_ind = sib.dendrogram_row.ordered_ind

# >>>> clustering of region-samples with regions including more than 10 samples
sample_present = sample_belong > 0
sample_present_count = np.sum(sample_present, axis=1)

sample_select10 = sample_present_count > 12
method = 'ward'
sib = sns.clustermap(sample_present[sample_select10,:])
plt.show()


# >>>> get the p-value of the peaks
# in each throw, the prob of that region be visited is 8%. what are the the odds of that region be visited 10 times, in 100 throws. and then there is an f value.

# what are teh odds of long region coverage
#alternative : {‘two-sided’, ‘greater’, ‘less’},
from scipy.stats import binom_test
p = binom_test(x= 11, n = 104, p= .08, alternative= 'greater')
print(p*counter)

for i in range(1,95):
    print(i)
    observed = sum(max_cov > i)
    print(sum(max_cov > i))
    p = binom_test(x= i, n = 104, p= .1, alternative= 'greater')
    expected = p*counter
    print(p*counter)

    print('FDR: %f' %(expected/observed))
    print('_____')

# we have 15094k peaks which are significant with FDR < .05, for 15+ values
# select these regions, how many samples contribute, what is the distance from center, and cluster
max_cov = np.max(region_coverage, axis = 1)
select_regions_ind = max_cov >=15

select_regions_belong = sample_belong[select_regions_ind,]
select_regions_belong[select_regions_belong>1] = 1
method = 'ward'
sib = sns.clustermap(select_regions_belong)
plt.show()

select_region_coverage = region_coverage[select_regions_ind,]
sns.heatmap(select_region_coverage[1:200,])
plt.show()



# TODO: get the location of the peaks, total number of the peaks etc. 

# TODO: that sanity check: which labels do not exist in any of the peaks. By theory, the low enhancers should not exsit in the peaks. But we don't know.

# Bits and pieces from the previous analyses
############################################
# how much of chr19 is covered with any enhancer label
#

# get the count bp covered in chr19 here.
file = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/testBatch105/fromAPI/ENCSR121RCG/chr19_whole.bed'
diff = 0
with open(file, 'r') as input:
    for line in input:
        fields = line.split()
        diff += int(fields[2]) - int(fields[1])

# Enhancer around genes
#########################################


# Enhancer around known enhancers
#########################################




# get the label from label from mnemonics: a dictionary from labels to terms
label_term_mapping = {}
mnemonics_file = sampleFolderAdd + 'mnemonics_v02.txt'
with open(mnemonics_file, 'r') as mnemonics:
    for line in mnemonics:
        #print(line)
        label = line.strip().split()[0]
        term = line.strip().split()[1]
        label_term_mapping[label] = term



# Draft

        

    p = subprocess.run(command, shell=True, capture_output = True)
    book = str(p.stdout)[2:-1]
    lines = book.splitlines()
    
    # run the command to get the lines with chr19 and the label that is enhancer
        command = "grep -E 'chr19.*\\t([%s]|13)_\' %s" %('297', segwayFile[0:-3])
    command = "grep \'chr19.*[%s]_\' %s > chr19_en.bed" %(label_string, segwayFile[0:-3])    
    # read the chr19 file line by line,

    # make a dict of tuples, for which the region belongs in. 


