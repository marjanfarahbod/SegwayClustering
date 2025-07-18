# the goal of this code is to record the enhancer regions within chromosome 19.
### sections ###
#
# 0. Initials
# 1. Get chr19 files 
# 2. First region analyses (obsolete, the whole section 2)
# 2.1 Getting the regions
# 2.2 Exploratory analyses for the region
# 3. [second region analysis] Chr vector, adding all the samples and some analyses there
# 3.1 adding all the smaples, geting the count of reg for each sample
# 3.2 coverage for count of repeats throught the chromosome
# 3.3 coverage and coverage of the overlap for each label
# 3.4 Examining the odd labels: labels with small overlap with others
# 4. peak identification
# 5. Peak label examination


#########################################
# 0. Initials
#########################################
import linecache
import pickle
import re
import numpy as np
import pandas as pd
import os
import shutil
import subprocess

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
dataSubFolder = 'testBatch105/fromAPI/'

inputFile = dataFolder + dataSubFolder + 'metaInfo.pkl'
with open(inputFile, 'rb') as f:
    annMeta = pickle.load(f)

sample_count = len(annMeta)

annAccessionList = list(annMeta.keys())
annAccession = annAccessionList[104]
print(annAccession)

inputFileName = 'all_annInfo_list.pkl'
inputFile = dataFolder + dataSubFolder + inputFileName
with open(inputFile, 'rb') as f:
    ann_info_list = pickle.load(f)

# get the mapping from accession to index
accession_index_map = {}
for ann in ann_info_list:
   accession_index_map[ ann['accession'] ] = ann['index']


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
    mnemonics_file = sampleFolderAdd + 'mnemonics_v04.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            print(line)
            label = line.strip().split()[0]
            term = line.strip().split()[1]
            label_term_mapping[label] = term

    enLabels = []
    for label in label_term_mapping.keys():
        if 'Enhancer' == label_term_mapping[label]:
            enLabels.append(label)

    one_digit_label = []
    two_digit_label = []
    for label in enLabels:
        if len(label) == 2:
            two_digit_label.append(label)
        else:
            one_digit_label.append(label)

    # making the label string for the label
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
    
    out = sampleFolderAdd + 'chr19_enhOnly.bed'
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

file = dataFolder + 'enhancer_regions_v1.pkl'
with open(file, 'rb') as f:
    result = pickle.load(f)

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
temp_ind = np.where(peaks <2)
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

#for i,k in enumerate(reversed(sib)):
#    print(i)
#    print(k)

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

#########################################
# 3. Chr vector, adding all the samples and some analyses there
#########################################

# 3.1 adding all the smaples, getting the count of reg for each sample
########################################

# 58597100
# 59000000

# arraylength
chr_length = 590000

chr_coverage = np.zeros(chr_length)
max_ind = 0
seg_end = 0
# do the analysis
sampleCoverage = np.zeros(len(annAccessionList)) # what percent is covered by the label for each sample
sampleMax = np.zeros(len(annAccessionList)) # what is the max length for ecah sample
sumHighCount = np.zeros(len(annAccessionList)) # sum of all labels greater than 2000 (20)
sumHighCount2x = np.zeros(len(annAccessionList)) # sum of all labels greater than 2000 (20)
disGreaterThan100 = np.zeros(len(annAccessionList)) # sum of all labels greater than 2000 (20)
prev_end = 10
disGreaterThan100_coverage = np.zeros(len(annAccessionList)) # sum of all labels greater than 2000 (20)
maxDists = np.zeros(len(annAccessionList)) # max dists for each sample
for annAccession in annAccessionList:

    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
    print(sampleFolderAdd)

    index = accession_index_map[annAccession]
    print(index)

    if(seg_end > max_ind):
        max_ind = seg_end

    out = sampleFolderAdd + 'chr19_enhOnly.bed'
    with open(out, 'r') as enFile:
        #line = enFile.readline()
        #print(line)

        #fields = line.split()
        #center = int(fie)
        for line in enFile:
            #sample_line_count +=1
            fields = line.split()
            seg_start = int(fields[1][0:-2])
            seg_end = int(fields[2][0:-2])-1
            seg_length = seg_end - seg_start 
            chr_coverage[seg_start:seg_end] += 1
            sampleCoverage[index] += seg_length
            if sampleMax[index] < seg_length:
                sampleMax[index] = seg_length
            if seg_length > 20:
                sumHighCount[index] +=1
            if seg_length > 40:
                sumHighCount2x[index] +=1

            if seg_start - prev_end >= 200:
                disGreaterThan100[index] +=1
                disGreaterThan100_coverage[index] += seg_start-prev_end-200
                if (seg_start-prev_end) > maxDists[index]:
                    maxDists[index] = seg_start-prev_end
                
            prev_end = seg_end

if(seg_end > max_ind):
    max_ind = seg_end

chr_coverage = chr_coverage[0:max_ind]
enhRatio = sampleCoverage/len(chr_coverage)

# plot the sorted ratio 

# 3.2 coverage for count of repeats throught the chromosome
######################################### 
# what percentage is zero?
print(sum(chr_coverage == 0)/max_ind)

# what percentage is one, two or zero:
print(sum(chr_coverage < 3)/max_ind)

coverage_count = np.zeros(sample_count) 
coverage_ratio = np.zeros(sample_count)
s_coverage = chr_coverage
for i in range(sample_count):
    print(i)
    coverage_count[i] = sum(s_coverage == i)
    coverage_ratio[i] = sum(s_coverage == i)/max_ind
    s_coverage = s_coverage[s_coverage>i]


cc_cum_sum = np.cumsum(coverage_count)
cr_cum_sum = np.cumsum(coverage_ratio)
import matplotlib.pyplot as plt
plt.plot(coverage_ratio[1:])
plt.plot(cr_cum_sum)
plt.plot([14, 14], [cr_cum_sum[14], cr_cum_sum[0]], 'k-', linewidth = .5)
plt.plot([0, 14], [cr_cum_sum[14], cr_cum_sum[14]], 'k-', linewidth = .5)

plt.plot([17, 17], [cr_cum_sum[17], cr_cum_sum[0]], 'k--', linewidth = .5)
plt.plot([0, 17], [cr_cum_sum[17], cr_cum_sum[17]], 'k--', linewidth = .5)
plt.xlabel('overlap count')
plt.ylabel('chr19 coverage')

plt.show()

# 3.3 coverage and coverage of the overlap for each label
########################################
# this is useful to idnetify odd labels, those which dwel in areas with lower overlap
# ideally, from each label we expect a certain coverage of the overlaps. 

# first, getting count of labels for each sample, for allocation

en_label_count = np.zeros(sample_count)
en_label_dict = {}
for annAccession in annAccessionList:

    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
    print(sampleFolderAdd)

    index = accession_index_map[annAccession]
    print(index)
    # load the mnomics, and extract the enhancer labels for mnemonics
    label_term_mapping = {}
    mnemonics_file = sampleFolderAdd + 'mnemonics_v02.txt'
    with open(mnemonics_file, 'r') as mnemonics:
        for line in mnemonics:
            #print(line)
            label = line.split()[0]
            term = line.split()[1]
            label_term_mapping[label] = term

    enLabels = []
    for label in label_term_mapping.keys():
        if 'Enhancer' in label_term_mapping[label]:
            enLabels.append(label)

    en_label_count[index] = len(enLabels)
    en_label_dict[index] = enLabels

sib = np.sort(en_label_count)
plt.plot(sib)
plt.show()

# NOTE: I used to do it for up to 20 overlaps, then I just switched to 85 for the record. However, 20 is good enough since most labels reside mostly in areas with less than 20 overlaps. 
# getting the indices in the chromosome for each count of overlap, we do it for up to 20
overlap_inds = {}
for i in range(85):
    print(i)
    overlap_inds[i] = np.where(chr_coverage == i+1)


en_label_count # for each sample, it has count of enhancer labels
en_label_dict # for each sample, it has list of enhancer labels
chr_coverage # to the length of the chromosome
total_label = int(sum(en_label_count))
label_coverages = np.zeros((total_label, 85)) # for each label, what portion is covered with the bps with overlap 0 to 85
label_ind = 0 # walking on the matrix 
label_IDs = {} # for the ind (label_ind), the sample and the name of the label
for annAccession in annAccessionList:

    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
    print(sampleFolderAdd)

    index = accession_index_map[annAccession]
    print(index)

    if en_label_count[index] == 0:
        continue

    # get the coverage for each of the labels in the sample
    temp_labels_coverage = np.zeros((int(en_label_count[index]), len(chr_coverage)))
    # getting the list of the labels for the sample
    en_label_list = en_label_dict[index]

    # getting the chromosome coverage for each label
    out = sampleFolderAdd + 'chr19_enh.bed'
    with open(out, 'r') as enFile:
        #line = enFile.readline()
        #print(line)

        #fields = line.split()
        #center = int(fie)
        for line in enFile:
            fields = line.split()
            seg_start = int(fields[1][0:-2])
            seg_end = int(fields[2][0:-2])-1
            label = fields[3].split('_')[0]
            label_sample_ind = en_label_list.index(label) # according to this list we go
            temp_labels_coverage[label_sample_ind, seg_start:seg_end] +=1 


    # get the coverage of the label for each count of overlap
    temp_coverage = np.zeros((len(en_label_list), 85))
    for i in range(len(overlap_inds)): 
        ois = overlap_inds[i][0]
        for j in range(len(en_label_list)):
            temp_coverage[j, i] = np.sum(temp_labels_coverage[j, ois])

    for i in range(len(en_label_list)):
        label_IDs[label_ind] = (index, en_label_list[i]) # the sample and the name of the label
        label_coverages[label_ind,] = temp_coverage[i,]
        label_ind += 1 # walking on the matrix 


# get the percentage of each overlap coverage for each of teh labels [1 to 20]
norm_lc = np.zeros(label_coverages.shape)
for i in range(total_label):
    norm_lc[i,] = label_coverages[i,]/sum(label_coverages[i,])

sns.heatmap(norm_lc)
plt.show()

countratio_lc = np.zeros(label_coverages.shape)
for i in range(20):
    countratio_lc[:,i] = label_coverages[:, i]/coverage_count[i+1]
    
import seaborn as sns
sns.heatmap(countratio_lc[:,0:20])
plt.show()

# get the ratio of label residing in non-significant regions
kado = np.sum(norm_lc[:, 0:15], axis = 1)
plt.plot(np.sort(kado))
plt.show()

print(np.where(kado >.60))

x = np.random.uniform(-.5, .5, len(kado))
plt.scatter(x, kado, alpha=.2)
plt.xlim((-5,5))
plt.show()


# 3.4 Examining the odd labels: labels with small overlap with others
#########################################################

for i in range(325):
    #rcum_sum = np.cumsum(norm_lc[i,])
    #plt.plot(rcum_sum, 'k', alpha = .07)
    plt.plot(norm_lc[i,], 'k', alpha = .07)
    maxy = max(norm_lc[i,])
    maxx = np.where(norm_lc[i,]==maxy)
    plt.plot(maxx, maxy, 'b.')

plt.show()

label_IDs_list = []
for i in range(len(label_IDs)):
    id_text = str(label_IDs[i][0]) + '___' + label_IDs[i][1]
    label_IDs_list.append(id_text)

plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/testBatch105/'
# method = 'ward'n
norm_lc[norm_lc>.11] = .11
df = pd.DataFrame(norm_lc, index = label_IDs_list)
sib = sns.clustermap(df, figsize=(6,40), col_cluster=False, cmap='magma', cbar_pos=(.02, .985,.03,.01), dendrogram_ratio=[.25,.01])
sns.set(font_scale=.4)
figfile = plotFolder + 'enhancerlabel_overlap_coverage.pdf'
plt.savefig(figfile)
plt.show()

dendro_ind = sib.dendrogram_row.reordered_ind


# TODO: here find the odd labels: I will do a pairwise overlap test for the odd labels

# get the index of less overlap labels. For each of them, do the overlap with the files.
book = dendro_ind[263:]
res_list = [label_IDs_list[i] for i in book]
print(res_list)

overlap_mat = np.zeros((len(res_list), len(res_list)))
for i,label1 in enumerate(res_list):
    # get the index
    print(label1)
    print(i)
    print('_________________')
    index1 = label1.split('___')[0]
    labelID1 = label1.split('___')[1]

    annAccession_1 = annAccessionList[int(index1)]
    sampleFolderAdd_1 = dataFolder + dataSubFolder + annAccession_1 + '/' 
    bed_1 = sampleFolderAdd_1 + 'chr19_enh.bed'

    # get the name of the .bed file for chr19
    for j in range(i+1, len(res_list)):
        label2 = res_list[j]
        print(label2)
        # get the index
        index2 = label2.split('___')[0]
        labelID2 = label2.split('___')[1]
        
        annAccession_2 = annAccessionList[int(index2)]
        sampleFolderAdd_2 = dataFolder + dataSubFolder + annAccession_2 + '/' 
        bed_2 = sampleFolderAdd_2 + 'chr19_enh.bed'

        with open(bed_1, 'r') as f_1, open(bed_2, 'r') as f_2:
            
            lines_1 = f_1.readlines()
            lines_2 = f_2.readlines()

        w1 = 0
        w2 = 0
        overlap = 0
        while w1 < len(lines_1) and w2 < len(lines_2):
            fields_1 = lines_1[w1].split()
            ll_1 = fields_1[3].split('_')[0]
            while ll_1 != labelID1 and w1<len(lines_1)-1:
                w1 += 1
                fields_1 = lines_1[w1].split()
                ll_1 = fields_1[3].split('_')[0]

            if ll_1 == labelID1:
                s1 = int(fields_1[1])
                e1 = int(fields_1[2])-100
            else:
                print('end of file 1')
                break
            
            fields_2 = lines_2[w2].split()
            ll_2 = fields_2[3].split('_')[0]
            while ll_2 != labelID2 and w2<len(lines_2)-1:
                w2 += 1
                fields_2 = lines_2[w2].split()
                ll_2 = fields_2[3].split('_')[0]

            if ll_2 == labelID2:
                s2 = int(fields_2[1])
                e2 = int(fields_2[2])-100
            else:
                print('end of file 2')
                break

            # if overlap, count
            if not(s1 > e2 or e1 < s2):
                overlap += min(e1, e2) - max(s1, s2)

            # the one with smaller end moves fw
            if e1 < e2:
                w1 +=1
            else:
                w2 +=1

        overlap_mat[i,j] = overlap
        overlap_mat[j,i] = overlap
                        
# get the lenght of each label
res_list_length = np.zeros(len(res_list))
for i,label in enumerate(res_list):
    # get the index
    print(label)
    index = label.split('___')[0]
    labelID = label.split('___')[1]

    annAccession = annAccessionList[int(index)]
    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/' 
    bed = sampleFolderAdd + 'chr19_enh.bed'

    with open(bed, 'r') as f:
        lines = f.readlines()
        for line in lines:
            fields = line.split()
            ll = fields[3].split('_')[0]
            if ll == labelID:
                s = int(fields[1])
                e = int(fields[2])-100
                res_list_length[i] += e - s            


# get the normalized overlap for each pair:
expected_overlap = np.zeros(overlap_mat.shape)
chr_length = 590000
for i in range(len(res_list_length)):
    for j in range(i+1, len(res_list_length)):
        expected_overlap[i, j] = (overlap_mat[i,j]/(chr_length*100))/((res_list_length[i]/(chr_length*100))*(res_list_length[j]/(chr_length*100)))
        expected_overlap[j,i] = expected_overlap[i,j]

#expected_overlap[expected_overlap>4] = 4

log_overlap = np.log10(expected_overlap+1)
dist = 1/(log_overlap+1)
np.fill_diagonal(dist,0)

from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram, linkage

condensed_dist = squareform(dist)
linkresult = linkage(condensed_dist, method='ward')

dn = dendrogram(linkresult)
print(dn['ivl'])
order_list = dn['ivl']

reordered = np.zeros(len(order_list))
for i,ind in enumerate(order_list):
    reordered[i] = int(ind)

reordered = reordered.astype(int)
res_list_reordered = [res_list[i] for i in reordered]

log_overlap2 = np.log10(expected_overlap+1)
overlap_reordered = log_overlap2[reordered,:]
overlap_reordered = overlap_reordered[:, reordered]

df = pd.DataFrame(overlap_reordered, index=res_list_reordered, columns = res_list_reordered)
sns.heatmap(overlap_reordered)
sns.heatmap(df)
sns.set(font_scale=.4)
figfile = plotFolder + 'enhancerlabel_clustering_of_odd_labels.pdf'
plt.savefig(figfile)

# saving all the results:
res = {}
res['idlist_ordered'] = res_list
res['dendro'] = dn
file = dataFolder + 'enhancerlabel_clustering_of_odd_labels.pkl'
with open(file, 'wb') as f:
    pickle.dump(res, f)


# each label, based on its coverage and the coverage of the repeat's coverage, has an odd of having overalp with repeat regions. Whichever label that is behaving oddly, missing or gaining too much, is of interest.
# TODO: catch this. But also do the code for the coverage stuff. 

# basically, I can keep a record of segments with various repeats, and for each segment, its repeat count and also which samples contributed. Then at each level, find the "odd" set of labels that appear together.

# here lets find who is this 5 repeat labels:
print(coverage_ratio[0:10])
sib = norm_lc[:, 5]
kado = np.where(sib > .17)
print(sib[kado])
for i in kado[0]:
    print(label_ID[i])
# NOTE: the 5 repeats are promoters and NOT enhancers.

sib = norm_lc[:, 0] 
kado = np.where(sib > .10)
print(sib[kado])
for i in kado[0]:
    print(label_ID[i])

# another factor is that how much of the label is dwelling in the lower coverage regions, and how much of it in the significant regions - that doesn't necessarily mean the label is wrong, rather than it has sth special going on

print(sum(coverage_ratio[0:10])) # about 76% is covered by this. so lets see how it looks for labels

sib = np.sum(norm_lc[:, 0:3], axis = 1)
halva = np.sort(sib)
plt.plot(halva)
plt.show()

kado = np.where(sib>.2)
for i in kado[0]:
    print(label_ID[i])
    #print(sib[i])
    #print(norm_lc[i,0:4])


# TODO: prune the odd labels based on the coverage

# contribution of each label to the low coverage regions (can we detect off labels, or off samples, or all that is hapening there is just expected distribution)

# 3.3 what is the count of individual enhancer events, and what is the distribution of distance between them
################################################################################

dis_dist = np.zeros(len(chr_coverage))
dis_enh = np.zeros(len(chr_coverage))
dist_counter = 0
enh_counter = 0
dist_length = 0
enh_length = 0
currentState = chr_coverage[0]
for i in range(1,len(chr_coverage)):
    previousState = currentState
    currentState = chr_coverage[i]
    if previousState == 0 and currentState > 0: # getting out of distance
        dis_dist[dist_counter] = dist_length
        dist_counter +=1
        dist_length = 0 # reset
    if currentState == 0: # in the distance
        dist_length +=1

    if previousState > 0 and currentState == 0: # getting out of enhancer
        dis_enh[enh_counter] = enh_length
        enh_counter +=1
        enh_length = 0
    if currentState > 0:
        enh_length +=1

dis_dist = np.sort(dis_dist[0:dist_counter])
dis_enh = np.sort(dis_enh[0:enh_counter])

fig, ax = plt.subplots()
plt.scatter(np.log10(dis_dist[1:]), np.log10(dis_enh))
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.xlim(-.5,5)
plt.ylim(-.5,5)
plt.show()

dist_qq = np.quantile(np.log10(dis_dist), np.linspace(0,1, 100))
enh_qq = np.quantile(np.log10(dis_enh), np.linspace(0,1, 100))

book = np.sort(dis_enh)
book = np.sort(np.log(dis_dist))

import matplotlib.pyplot as plt
plt.scatter(range(len(book)), book, alpha=.3, s=8)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)

plt.show()

fig, ax = plt.subplots()
ax.scatter(dist_qq, enh_qq)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)

ax.set_xlabel('x-axis', fontsize = 12)
ax.set_ylabel('y-axis', fontsize = 10)

fig, ax = plt.subplots()
ax.scatter(np.log10(dis_enh), range(len(dis_enh)), alpha=.3, s=8)
  
plt.show()

########################################
# 4. peak identification
########################################

# coverage obtained, now the peaks. My assumption is that most peaks narrow down significantly before visiting another peak.
book = np.zeros(1000)
book[0] = chr_coverage[0]
for i in range(999):
    book[i+1] = book[i] + chr_coverage[i+1]

plt.plot(book)
plt.show()
# a peak is where the the value is higher or equal in the neighbouring indices. Therefore, some peaks can be wide and some can be narrow.
# a peak region is idenfied by a region starting with the increasing values towards a peak, and decreasing from it. So as soon as we pass 5 overlaps, we start counting. This is not covering the significance level by any means, but rather keeps the overlap. I keep for each peak: max, start, end. We will fill out the labels later.

# my guess is that we will not have more than 50k peaks.
peak_thr = 10
peak_sig_thr = 15
peakInfo = np.zeros((20000, 3))
in_peak_region = False
peak_counter = 0
max_peak = 0
peak_starts = 0
peak_ends = 0
for i in range(1,len(chr_coverage)):
    #print(i)
    #print(peak_counter)
    #pprint('____')
    if chr_coverage[i] < peak_thr:
        if in_peak_region:
            #record the peak
            if max_peak >= peak_sig_thr:
                peak_ends = i-1
                peakInfo[peak_counter, 0] = peak_starts
                peakInfo[peak_counter, 1] = peak_ends
                peakInfo[peak_counter, 2] = max_peak
                peak_counter +=1;
                # reset the peak
                in_peak_region = False
                max_peak = 0
            else:
                in_peak_region = False
                max_peak = 0
        continue

    if chr_coverage[i] >= peak_thr:
        if in_peak_region:
            if chr_coverage[i] > max_peak:
                max_peak = chr_coverage[i]
        else:
            # initiate peak
            in_peak_region = True
            peak_starts = i
            max_peak = chr_coverage[i]

peakInfo = peakInfo[0:peak_counter, ]

peak_length = peakInfo[:, 1] - peakInfo[:, 0]

plt.hist(peak_length)
plt.show()

high_length_ind = np.where(peak_length>20)
print(len(high_length_ind[0]))
print(peakInfo[3386,])
peakInd = high_length_ind[0][0]
ps = int(peakInfo[peakInd, 0])
pe = int(peakInfo[peakInd, 1])
kado = np.mean(chr_coverage[ps:pe])
print(kado)
plt.plot(chr_coverage[ps-400:pe+400])
plt.show()

hl_means = np.zeros(len(high_length_ind[0]))
for i in range(len(high_length_ind[0])):
    peakInd = high_length_ind[0][i]
    ps = int(peakInfo[peakInd, 0])
    pe = int(peakInfo[peakInd, 1])
    hl_means[i] = np.mean(chr_coverage[ps:pe])
    
plt.hist(hl_means)
plt.show()

other_length_ind = np.where(peak_length<=20)
# distribution of peak length:
plt.hist(peak_length[other_length_ind[0]])
plt.show()
ol_means = np.zeros(len(high_length_ind[0]))
for i in range(len(high_length_ind[0])):
    peakInd = high_length_ind[0][i]
    ps = int(peakInfo[peakInd, 0])
    pe = int(peakInfo[peakInd, 1])
    hl_means[i] = np.mean(chr_coverage[ps:pe])
    
plt.hist(hl_means)
plt.show()


# distance between peaks
peak_dis = np.zeros(peakInfo.shape[0])
for i in range(1,len(peak_dis)):
    peak_dis[i-1] = peakInfo[i,0] - peakInfo[i-1, 1]

print(sum(peak_dis<2))
peak_dis[peak_dis > 200] = 200
plt.hist(peak_dis)
plt.show()

sib = np.cumsum(peak_dis)
plt.plot(sib)
plt.show()

plt.plot(chr_coverage[2459:6000])
plt.show()

####
sib = np.sort(chr_coverage)
plt.plot(sib)
plt.show()

plt.plot(chr_coverage[0:10000])
plt.show()

print(len(high_length_ind[0]))
high_length_max = peakInfo[high_length_ind, 2]

plt.hist(high_length_mean[0])
plt.show()

#########################################
# 5. Peak label examination 
#########################################
print(peakInfo.shape)
label_IDs

label_ID_reverse_dic = {}
for label in label_IDs:
    ID = label_IDs[label]
    label_ID_reverse_dic[ID] = label

# peak presence
peakCount = peakInfo.shape[0] 
labelCount = label_coverages.shape[0]
peak_presence = np.zeros((peakCount, labelCount))

peak_presence_sample = np.zeros((peakCount, sample_count))

# walk through the labels and lines, fill them as fit
for annAccession in annAccessionList:

    sampleFolderAdd = dataFolder + dataSubFolder + annAccession + '/'
    print(sampleFolderAdd)

    index = accession_index_map[annAccession]
    print(index)

    if en_label_count[index] == 0:
        continue

    # get the coverage for each of the labels in the sample
    # temp_peak_coverage = np.zeros((int(en_label_count[index]), peakCount))
    # getting the list of the labels for the sample
    # en_label_list = en_label_dict[index]

    # getting the chromosome coverage for each label
    peak_walk = 0
    out = sampleFolderAdd + 'chr19_enh.bed'
    with open(out, 'r') as enFile:

        # we go through lines (segments), as long as the seg_start < peak_end
        # for each step, we check if segment and peak overlap. once overalp, we move fw on the peaks as long as they do.

        # then we move fw on the peak, as long as seg_start < peak_end

        # get the first peak
        current_peak = peakInfo[peak_walk,:]
        peak_start = current_peak[0]
        peak_end = current_peak[1]

        line = enFile.readline()
        line = enFile.readline()
        print(line)
        fields = line.split()
        seg_start = int(fields[1][0:-2])
        seg_end = int(fields[2][0:-2])-1
        label = fields[3].split('_')[0]

        while seg_start > peak_end: # taking the peak fw
            peak_walk +=1
            current_peak = peakInfo[peak_walk,:]
            peak_start = current_peak[0]
            peak_end = current_peak[1]
            

        while ((seg_start < peak_end) or (peak_walk < peakCount-1)) and peak_walk != peakCount-1: # while there is a peak ahead of a segment
            # print(peak_walk)
            if not((peak_start > seg_end) or (peak_end < seg_start)): # if overlap, record the segment, move the peak fw until passed the segment (so seg_start < peak_end)
                # get label index
                label_ID = (index, label)
                l_index = label_ID_reverse_dic[label_ID]
                # add segment to the peak
                peak_presence[peak_walk, l_index] = 1
                peak_walk +=1  # move one peak fw
                #print(peak_walk)
                peak_presence_sample[peak_walk, index] = 1
                current_peak = peakInfo[peak_walk,:]
                peak_start = current_peak[0]
                peak_end = current_peak[1]

                # while overlap: move peak fw, add segment to the peak
                while not((peak_start > seg_end) or (peak_end < seg_start)) and peak_walk < (peakCount-1):

                    peak_presence[peak_walk, l_index] = 1
                    peak_walk +=1  # move one peak fw
                    print(peak_walk)

                    current_peak = peakInfo[peak_walk,:]
                    peak_start = current_peak[0]
                    peak_end = current_peak[1]

                
            else: # not overlap, seg_start < peak_end
                # move one segment fw
                try:
                    line = enFile.readline()
                    fields = line.split()
                    seg_start = int(fields[1][0:-2])
                    seg_end = int(fields[2][0:-2])-1
                    label = fields[3].split('_')[0]
                except IndexError:
                    print('end of file')
                    break
                

                while (seg_start > peak_end) and peak_walk < (peakCount-1):
                    # move one peak fw
                    peak_walk +=1  # move one peak fw
                    current_peak = peakInfo[peak_walk,:]
                    peak_start = current_peak[0]
                    peak_end = current_peak[1]


print(peak_presence.shape)
method = 'ward'
sib = sns.clustermap(peak_presence[1:4000,], cmap='binary', method=method)
plt.show()

method = 'ward'
sib = sns.clustermap(peak_presence_sample[1:4000,], cmap='binary')
plt.show()
reo_samples = sib.dendrogram_col.reordered_ind

# TODO: make decision about mislabeld labels
# TODO: train the classifier and get the new mnemonics
# TODO: do the clustering and get the labels for samples 

# print the indices. 

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



        #line = enFile.readline()
        #print(line)

        #fields = line.split()
        #center = int(fie)
        for line in enFile:
            fields = line.split()
            seg_start = int(fields[1][0:-2])
            seg_end = int(fields[2][0:-2])-1
            label = fields[3].split('_')[0]

            while 
                # while the peak and segment overlap
                # todo: add segmen to the peak,
                # move to the next peak

                
            # get the index of the label, if it is within the peak range, to fill the coverage matrix


            # do the double coverage code 
### TODO: we are here ######


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


