# TODO : this file needs major clean up
# this is to do the experiment and test for checking the probabilities for the enhances
import pickle
import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt

# for each sample, load the prob matrix

dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
dataSubFolder = 'testBatch105/fromAPI/'

# for each sample we want the prob and the labels be generated. See how the probs file was made.

inputFileName = 'all_annInfo_list.pkl'
inputFile = dataFolder + dataSubFolder + inputFileName
with open(inputFile, 'rb') as f:
    annInfo_list = pickle.load(f)


term_list = ['Enhancer', 'Enhancer_low', 'Promoter','Promoter_flanking', 'CTCF', 'K9K36', 'Bivalent', 'Transcribed', 'FacultativeHet', 'ConstitutiveHet', 'Quiescent']

label_color_mapping = {}
label_color_mapping['Promoter'] = [[255, 0, 0]]
label_color_mapping['Promoter_flanking'] = [[255, 69, 0]] # promoter
label_color_mapping['Enhancer'] = [[220, 210, 0], [195, 180, 0]]
label_color_mapping['Enhancer_low'] = [[255, 255, 0],[180, 160, 0]] # enhancer
label_color_mapping['CTCF'] = [[102, 205, 170]] # CTCF
label_color_mapping['Transcribed'] = [[0, 128, 0], [0, 170, 0]]# Transcribed
label_color_mapping['FacultativeHet'] = [[200, 10, 200], [150, 5, 150]] # FacHet
label_color_mapping['ConstitutiveHet'] = [[138, 145, 208], [180, 195, 250], [100, 120, 180]] # consHet
label_color_mapping['Quiescent'] = [[200, 200, 220], [230, 230, 230]] # Quis
label_color_mapping['K9K36'] = [[180, 180, 180]]
label_color_mapping['Bivalent'] = [[0,0,0]]

# for each of the samples, keep the label, and list of probabilities
all_label_lists = {}
med_values = np.zeros(len(annInfo_list))
q1_values = np.zeros(len(annInfo_list))
low_mean_vals = np.zeros(len(annInfo_list))
# median is sorted by the index, index and accession are mapped. So once I resort the median, I can get the accession from the index-accession map
index_accession = {}
for ann in annInfo_list:

    index = ann['index']
    print(index)
    sampleFolder =  ann['accession']

    index_accession[index] = sampleFolder

    sampleFolderAdd = dataFolder + dataSubFolder + sampleFolder + '/'

    mnemonics_file = sampleFolderAdd + 'mnemonics_v02.txt'
    df_m = pd.read_csv(mnemonics_file, sep='\t')

    probs_file = sampleFolderAdd + 'probs_v02.csv'
    df_p = pd.read_csv(probs_file)

    # extracting the maximum prob for each of the labels
    df_p = df_p.drop('Unnamed: 0', axis =1)
    label_names = df_p.idxmax(axis=1)
    label_prob = df_p.max(axis=1)

    lp_array = label_prob.to_numpy()
    med_values[index] = np.median(lp_array)
    q1_values[index] = np.quantile(lp_array, .25)
    vals = np.sort(lp_array)
    low_mean_vals[index] = np.mean(vals[0:5])

    label_list = []
    for i,item in enumerate(label_names):
        color = label_color_mapping[item][0]
        prob = label_prob[i]
        this_label = (i, item, prob, color)
        label_list.append(this_label)

    
    all_label_lists[sampleFolder] = label_list


# plotting the values
# sorting the median's from low to high
med_sorted = np.sort(med_values)
med_sorted_arg = np.argsort(med_values)

low_mean_sorted = np.sort(low_mean_vals)
low_mean_arg = np.argsort(low_mean_vals)

q1_sorted = np.sort(q1_values)
q1_sorted_arg = np.argsort(q1_values)


# size of the figure, grid on, black background, increase the circle sizes, add the sample number
# save the plots for the three values - see which one you wan to pick
# compare the low values with chromhmm Quiescents.
# finish this and get the updated classifier data
fig, axs = plt.subplots(1, 1, figsize=[12, 8])
var_range = .17
for i,index in enumerate(low_mean_arg[78:105]):
    print(i)
    accession = index_accession[index]
    label_list = all_label_lists[accession]

    x = i
    for item in label_list:
        var = random.uniform(-var_range, var_range)
        carr = np.array(item[3])
        plt.scatter(x + var, item[2], s = 45, color = carr/256, alpha=.85)
        plt.axvline(x + .5, color='lightgray', linestyle='--')
        
plt.xticks(range(27), labels = low_mean_arg[78:], fontsize=10, rotation=90)
plt.xlabel('sample ID')
plt.ylabel('probs')
plt.show()


# samples with mislabeled enhancers:

s1 = [4, 45, 47, 24, 12,  30]

for i in s1:
    print(np.where(med_sorted_arg == i))

# samples with likely mislabedl enhancers:

s2 = [46, 98, 12]

for i in s2:
    print(np.where(med_sorted_arg == i))


# samples with odd k9k36 suspect

s3 =  [69, 86, 57, 2, 83, 37, 40]
for i in s3:
    print(np.where(med_sorted_arg == i))

# we have identified bad samples

for i in range(105):
    print(i)


vlist = []
low_mean_list = []
index_list = []
for j, i in enumerate(low_mean_arg):
    if i in list(chmm_quis_coverage.keys()):
        if chmm_quis_coverage[i] > 0:
            vlist.append(chmm_quis_coverage[i])
            low_mean_list.append(low_mean_sorted[j])
            index_list.append(i)

fig, axs = plt.subplots(1, 1, figsize=[18, 6])            
plt.plot(vlist)
plt.xticks(range(101), labels = index_list, fontsize=10, rotation=90)
plt.grid()
plt.ylabel('chmm quiescent coverage')
plt.show()

vlistarr = np.array(vlist)
low_mean_listarr = np.array(low_mean_list)

sib = np.correlate(vlistarr, low_mean_listarr)
from scipy.stats import pearsonr
corr, _ = pearsonr(vlistarr, q1_listarr)

# the quis by itself is not determinent, is there any other label that is determinent? is it similar in ours?

chmm_quis_coverage = {0: .74, 1: .69,2: .66,3: .85,
                      4: .45,5: 0,6: .85,7: .86,8: .58,9: .82,10: .84,
                      11: .51,12: .64,13: .66,14: .85,15: .63,16: .68,17: .71,18: .86,
                      19: 0,20: .81,21: .46,22: .82,23: .80,24: .76,25: .66,26: .82,
                      27: .84,28: .88,29: .63,30: .70,31: .0,32: .69,33: .73,34: .58,
                      35: .82,36: .6,37: .84,38: .81,39: .86,40: .7,41: .78,42: .75,
                      43: .76,44: .58,45: .63,46: .56, 47: .92,48: .85,49: .79,50: .66,
                      51: .57,52: .83,53: .86,54: .86,55: .69,56: .71,57: .85,58: 0,
                      59: .76,60: .85,61: .82,62: .80,63: .75,64: .74,65: .44,66: .75,
                      67: .79,68: .88,69: .86,70: .63,71: .90,72: .55,73: .76,74: .67,
                      75: .58,76: .82,77: .62,78: .73,79: .69,80: .63,81: .61, 82: .78,
                      83: .81,84: .88,85: .69,86: .90,87: .85,88: .39,89: .7,90: .6,
                      91: .76,92: .9,93: .84,94: .78,95: .69,96: .85,97: .83,98: .78,
                      99: .78,100: .60,101: .58,102: .81,103: .90,104: .83}

# 57-2 is ambiguous labels
# let's find good samples ambiguous labels

### Get for each sample, what percent of the genome was done with more than x% confidence

#  we need the label coverages

# from ann info, for each label get the coverage.
# from probs, for each label get the probability.
# set the x as the percentage. 

index_accession = {}
x = .3 # probability
sample_pCoverage = {}
spc_array = np.zeros(len(annInfo_list))
for ann in annInfo_list:

    index = ann['index']
    print(index)
    sampleFolder =  ann['accession']

    index_accession[index] = sampleFolder

    sampleFolderAdd = dataFolder + dataSubFolder + sampleFolder + '/'

    mnemonics_file = sampleFolderAdd + 'mnemonics_v02.txt'
    df_m = pd.read_csv(mnemonics_file, sep='\t')

    probs_file = sampleFolderAdd + 'probs_v02.csv'
    df_p = pd.read_csv(probs_file)
    df_p = df_p.drop('Unnamed: 0', axis =1)
    label_prob = df_p.max(axis=1)

    label_list = list(annInfo_list[index]['segway_anns']['clusters'].keys())
    label_bp_counts = np.zeros(len(label_list))

    for label in label_list:
        label_label = int(annInfo_list[index]['segway_anns']['clusters'][label].cluster)
        label_bp_counts[label_label] = int(annInfo_list[index]['segway_anns']['clusters'][label].bp_count)

    total_bps = sum(label_bp_counts)
    genome_coverage = label_bp_counts/total_bps

    xp_coverage = 0
    for i,p in enumerate(label_prob):
        if p > x:
            xp_coverage += genome_coverage[i]

    sample_pCoverage[index] = xp_coverage
    spc_array[index] = xp_coverage


fig, axs = plt.subplots(1, 1, figsize=[18, 6])            
plt.plot(spc_array[low_mean_arg])
plt.xticks(range(105), labels = index_list, fontsize=10, rotation=90)
plt.grid()
plt.ylabel('percent of the genome predicted with > 0.3 prob')
plt.show()

# categories of low probs:
# enhancer, enhancerlow confusion: if sums greater than 50 with the other
# promoter, promoterflank confusion: if sums greater than 50 with one other partner
# quis, fachet, conshet confusion: with which partner sums greater than 50? 

list1 = [] # the sample
list2 = [] # the label
list3 = [] # the prob
list4 = [] # the term
list5 = [] # the next term
list6 = [] # the added prob

for ann in annInfo_list:

    index = ann['index']
    print(index)
    sampleFolder =  ann['accession']

    index_accession[index] = sampleFolder

    sampleFolderAdd = dataFolder + dataSubFolder + sampleFolder + '/'

    probs_file = sampleFolderAdd + 'probs_v02.csv'
    df_p = pd.read_csv(probs_file)
    df_p = df_p.drop('Unnamed: 0', axis =1)
    label_prob = df_p.max(axis=1)
    label_term = df_p.idxmax(axis=1)

    for i,prob in enumerate(label_prob):
        if prob < .4:
            term = label_term[i]
            book = np.array(df_p.loc[i])
            argsecondp = np.argsort(book)[-2]
            secondp = book[argsecondp]
            secondterm = list(df_p.columns)[argsecondp]

            list1.append(index)
            list2.append(i)
            list3.append(prob)
            list4.append(term)
            list5.append(secondterm)
            list6.append(secondp+prob)

print(len(list3))
sib = np.array(list3)
print(sum(sib>.3))


termIndices = {}
for term in term_list: 
    indices = [i for i, x in enumerate(list4) if x == term]
    termIndices[term] = indices
    print('term %s, %d' %(term, len(indices)))

# for each term, get the count of other terms and how many of them had the sum prob greater than 50 or 65

term1_based = {}
throwAwayCount = 0
for term1 in term_list:
    print('_______term1 ' + term1)
    term2_based = {}
    list5p = [list5[i] for i in termIndices[term1]] # term
    list6p = [list6[i] for i in termIndices[term1]] # second prob
    list3p = [list3[i] for i in termIndices[term1]] # first prob
    list1p = [list1[i] for i in termIndices[term1]] # first prob
    list2p = [list2[i] for i in termIndices[term1]] # first prob
    for term2 in term_list:
        print('__' + term2)
        tempind = [i for i, x in enumerate(list5p) if (x == term2)]
        p3meet = 0
        p4meet = 0
        nomeet = 0
        print(len(tempind))
        p3meet_labels = []
        p4meet_labels = []
        nomeet_labels = []
        total40 = 0
        total30 = 0
        for i in tempind:
            
            # how many have prob < .3, greater than .5 with the second
            sampleQindex = np.where(low_mean_arg == list1p[i])[0][0]
            if sampleQindex < 15:
                throwAwayCount +=1
            else:
                if list3p[i] > .3:
                    total40+=1
                    
                if list3p[i] < .3:
                    total30+=1
        
                if list3p[i] < .3 and list6p[i] > .45:
                    p3meet +=1
                    l = (list1p[i], list2p[i])
                    p3meet_labels.append(l)
                else:
                    if list3p[i] > .3 and list6p[i] > .55:
                        p4meet +=1
                        l = (list1p[i], list2p[i])
                        p4meet_labels.append(l)
                    else:
                        nomeet +=1
                        l = (list1p[i], list2p[i])
                        nomeet_labels.append(l)
        # how many have prob < .4, >.3, greater than .6 with second
        if p3meet > 0 or p4meet>0 or nomeet >0:
            term2_based[term2] = (total40, total30, p3meet, p4meet, nomeet, p3meet_labels, p4meet_labels, nomeet_labels)

    term1_based[term1] = term2_based
    
# how many labels are lower than .3, who is the next candidate and how up it goes

for term1 in term_list:
    print('_____' + term1)
    for term2 in list(term1_based[term1].keys()):
        print('_' + term2)
        print(term1_based[term1][term2])
    print('\n')


# how many labels are lower than .4 and higher than .3 - who is their next candidate and how up it goes

print(np.where(low_mean_arg == 85))

for ann in annInfo_list:

    index = ann['index']
    #print(index)
    sampleFolder =  ann['accession']

    if sampleFolder == 'ENCSR950GQD':
        print(index)


badsampleList1 = [ 86,  20,  22,  69,  57,   9,  92,  37,  54,  39, 103,  82,  23, 27,  25]
for ann in annInfo_list:
    index = ann['index']
    if index in badsampleList1:
        print(ann['accession'])
    
##### Bad samples for the May run
########################################

ID_mnemonics = {}
low_mean_vals = np.zeros(len(runIDs))
for i,runID in enumerate(runIDs):

    print(runID)

    if runID == '5857d68b-e559-4776-9c12-a6e10aea7f76': # this runID doesn't have the signal_dist file
        continue

    sampleFolderAdd = dataFolder + dataSubFolder + runID + '/'
    mnemfile = sampleFolderAdd + 'probs_v03.csv'
    df_p = pd.read_csv(mnemfile)

    
    # extracting the maximum prob for each of the labels
    df_p = df_p.drop('Unnamed: 0', axis =1)
    label_names = df_p.idxmax(axis=1)
    label_prob = df_p.max(axis=1)

    lp_array = label_prob.to_numpy()
    vals = np.sort(lp_array)
    low_mean_vals[i] = np.mean(vals[0:5])



outputFile = dataFolder + dataSubFolder + 'classifier_output.pkl'
with open(outputFile, 'wb') as output:
    pickle.dump(ID_mnemonics, output)



