# this is to pars the segway parameter batch run output. For a given sample folder, the parameter setting can be extracted from the name of the folder and the segment count is the [3?]rd line of segment_sizes.tab file in length_distribution folder.

import pandas as pd
import os
import matplotlib.pyplot as plt

results = pd.DataFrame(columns = ['ps', 'stws', 'tw', 'meanlen' ,'scount'], index = range(1260))


##############################
# write a module to get the folder address and return the result file
##############################

# for a given sample
sampleName = 'DND41'
sampleName = 'H1'
sampleName = 'adrenal_gland'
resultFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/BatchRuns_parameter/%s' %(sampleName)

# fetch the list of directories for that sample
dirList = os.listdir(resultFolder)

# pars the name of directories, fill the results
for dirName in dirList:
    nameSplit = dirName.split('_')

    if nameSplit[0][0:3] == 'seg':
        index = int(nameSplit[0][9:])
        ps = float(nameSplit[4])
        stws = float(nameSplit[5])
        tw = float(nameSplit[6])

        segSizeFile = '%s/%s/length_distribution/segment_sizes.tab' %(resultFolder, dirName)
        try:
            with open(segSizeFile, 'r') as file:
                line = file.readline()
                infoLine = file.readline().split()
                scount = int(infoLine[1])
                meanLen = float(infoLine[2])
                results.loc[index] = [ps, stws, tw, meanLen ,scount]
        except FileNotFoundError:
            results.loc[index] = [-1, -1, -1, -1, -1]


# save results into a tab delimited file with a header
outputFile = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/BatchRuns_parameter/%s.csv' %(sampleName)
results.to_csv(outputFile, sep='\t', header = ['ps', 'stws', 'tw', 'meanLen' ,'scount'])

adrenalRes = results.copy()

dndRes = results.copy()

h1Res = results.copy()

##############################
# plotting 
##############################

plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/'

# first plot: for a given ps, get the tw, and stws , get the length
# get all the rows with ps = 0
# reference 
ps_list = [0, .1, .5, 1, 10, 100, 1000]
# seg transition weight scale
stws_list = [0, .001, .01, .5, 1, 10, 50, 100, 500, 1000, 5000, 10000]
# track weight
tw_list = [.0005, .001, .002, .005, .008, .01, .02, .05, .08, .1, .5, 1, 10, 50, 100]

# >>>>>>>>>>>>>>>>>>>>
# the line plot
# >>>>>>>>>>>>>>>>>>>>

# selecting based on the ps
# subtable = results.loc[results['ps'] == 0]
maxCount = math.log2(max(results['meanlen']))

m1 = math.log2(max(adrenalRes['meanlen'])) # the max is equal, so one is enough

for plotCount in range(len(ps_list)):

    fig = plt.figure(figsize=(12,10))
    ps = ps_list[plotCount]

    x = list(range(len(stws_list)))
    y = list(range(int(m1)+2))

    # >>>>>>>>>> first plot: adrenal gland
    ax = fig.add_subplot(1,3,1) #(111)
    ax.set_ylabel('mean segment length (bp)', fontsize = 18)
    plotRes = adrenalRes
    
    subtable = plotRes.loc[plotRes['ps'] == ps]

    for tw in tw_list:    # for each tw, I want to plot the things
        subsubtable = subtable.loc[subtable['tw'] == tw]
    
    for row in range(len(subtable)):
        stws = subtable.iloc[row][1] # getting stws value
        tw = subtable.iloc[row][2] # getting tw value
        
        i = stws_list.index(stws)
        j = tw_list.index(tw)

        thisCount = math.log2(subtable.iloc[row][3])
        halfLen = ((thisCount/maxCount)/2)*.75
        plt.plot([i-halfLen, i+halfLen], [j, j], color='olive', linewidth=4)


    plt.title('prior_strength = %f' %(ps))    
    plt.xticks([-1] + x + [len(x)], [0] + stws_list + [0])
    plt.yticks([-1] + y + [len(y)], [0] + tw_list + [0])
    plt.grid()
    axt = ax.twinx()
    plt.yticks([-1] + y + [len(y)], [0] + tw_list + [0])
    plt.grid()
    #plt.show()

    # >>>>>>>>>> second plot: H1
    ax = fig.add_subplot(3,1,2) #(111)
    ax.set_ylabel('track_weight', fontsize = 18)
    plotRes = h1Res
    subtable = plotRes.loc[plotRes['ps'] == ps]
    for row in range(len(subtable)):
        stws = subtable.iloc[row][1] # getting stws value
        tw = subtable.iloc[row][2] # getting tw value
        
        i = stws_list.index(stws)
        j = tw_list.index(tw)

        thisCount = math.log2(subtable.iloc[row][3])
        halfLen = ((thisCount/maxCount)/2)*.75
        plt.plot([i-halfLen, i+halfLen], [j, j], color='olive', linewidth=4)


    plt.title('prior_strength = %f' %(ps))    
    plt.xticks([-1] + x + [len(x)], [0] + stws_list + [0])
    plt.yticks([-1] + y + [len(y)], [0] + tw_list + [0])
    plt.grid()
    axt = ax.twinx()
    plt.yticks([-1] + y + [len(y)], [0] + tw_list + [0])
    plt.grid()

    # >>>>>>>>>> second plot: DND41
    ax = fig.add_subplot(3,1,3) #(111)
    ax.set_ylabel('track_weight', fontsize = 18)
    plotRes = dndRes
    subtable = plotRes.loc[plotRes['ps'] == ps]
    for row in range(len(subtable)):
        stws = subtable.iloc[row][1] # getting stws value
        tw = subtable.iloc[row][2] # getting tw value
        
        i = stws_list.index(stws)
        j = tw_list.index(tw)

        thisCount = math.log2(subtable.iloc[row][3])
        halfLen = ((thisCount/maxCount)/2)*.75
        plt.plot([i-halfLen, i+halfLen], [j, j], color='olive', linewidth=4)


    plt.title('prior_strength = %f' %(ps))    
    plt.xticks([-1] + x + [len(x)], [0] + stws_list + [0])
    plt.yticks([-1] + y + [len(y)], [0] + tw_list + [0])
    plt.grid()
    axt = ax.twinx()
    plt.yticks([-1] + y + [len(y)], [0] + tw_list + [0])
    plt.grid()

    plt.savefig(plotFolder + '02threeTissues_ps_%.4f.pdf' %(ps))

# wegiht scale vs seg transition
# >>>>>>>>>>>>>>>>>>>>

for plotCount in range(len(ps_list)):

    fig = plt.figure(figsize=(12,10))
    ps = ps_list[plotCount]

    y = list(range(len(stws_list)))
    x = list(range(len(tw_list)))

    # >>>>>>>>>> first plot: adrenal gland
    ax = fig.add_subplot(3,1,1) #(111)
    ax.set_ylabel('seg_transition', fontsize = 18)
    plotRes = adrenalRes
    subtable = plotRes.loc[plotRes['ps'] == ps]
    for row in range(len(subtable)):
        stws = subtable.iloc[row][1] # getting stws value
        tw = subtable.iloc[row][2] # getting tw value
        
        j = stws_list.index(stws)
        i = tw_list.index(tw)

        thisCount = math.log2(subtable.iloc[row][3])
        halfLen = ((thisCount/maxCount)/2)*.75
        plt.plot([i-halfLen, i+halfLen], [j, j], color='olive', linewidth=4)


    plt.title('prior_strength = %f' %(ps))    
    plt.xticks([-1] + x + [len(x)], [0] + tw_list + [0])
    plt.yticks([-1] + y + [len(y)], [0] + stws_list + [0])
    plt.grid()
    axt = ax.twinx()
    plt.yticks([-1] + y + [len(y)], [0] + stws_list + [0])
    plt.grid()
    #plt.show()

    # >>>>>>>>>> second plot: H1
    ax = fig.add_subplot(3,1,2) #(111)
    ax.set_ylabel('seg_transition', fontsize = 18)
    plotRes = h1Res
    subtable = plotRes.loc[plotRes['ps'] == ps]
    for row in range(len(subtable)):
        stws = subtable.iloc[row][1] # getting stws value
        tw = subtable.iloc[row][2] # getting tw value
        
        j = stws_list.index(stws)
        i = tw_list.index(tw)

        thisCount = math.log2(subtable.iloc[row][3])
        halfLen = ((thisCount/maxCount)/2)*.75
        plt.plot([i-halfLen, i+halfLen], [j, j], color='olive', linewidth=4)


    plt.xticks([-1] + x + [len(x)], [0] + tw_list + [0])
    plt.yticks([-1] + y + [len(y)], [0] + stws_list + [0])
    plt.grid()
    axt = ax.twinx()
    plt.yticks([-1] + y + [len(y)], [0] + stws_list + [0])
    plt.grid()

    # >>>>>>>>>> second plot: DND41
    ax = fig.add_subplot(3,1,3) #(111)
    ax.set_ylabel('seg_transition', fontsize = 18)
    plotRes = dndRes
    subtable = plotRes.loc[plotRes['ps'] == ps]
    for row in range(len(subtable)):
        stws = subtable.iloc[row][1] # getting stws value
        tw = subtable.iloc[row][2] # getting tw value
        
        j = stws_list.index(stws)
        i = tw_list.index(tw)

        thisCount = math.log2(subtable.iloc[row][3])
        halfLen = ((thisCount/maxCount)/2)*.75
        plt.plot([i-halfLen, i+halfLen], [j, j], color='olive', linewidth=4)


    plt.xticks([-1] + x + [len(x)], [0] + tw_list + [0])
    plt.yticks([-1] + y + [len(y)], [0] + stws_list + [0])
    plt.grid()
    axt = ax.twinx()
    plt.yticks([-1] + y + [len(y)], [0] + stws_list + [0])
    plt.grid()
    ax.set_xlabel('track_weight', fontsize = 18)
    plt.savefig(plotFolder + '03threeTissues_ps_%.4f.pdf' %(ps))

# >>>>>>>>>>>>>>>>>>>>
# same thing with lines per length
# >>>>>>>>>>>>>>>>>>>>

# selecting based on the ps
# subtable = results.loc[results['ps'] == 0]
maxMean = math.log10(max(results['meanlen']))

for plotCount in range(len(ps_list)):

    fig = plt.figure(figsize=(16,10))
    ps = ps_list[plotCount]

    x = list(range(len(stws_list)))
    y = list(range(len(tw_list)))
    # y = list(range(len(tw_list)))

    # >>>>>>>>>> first plot: adrenal gland
    ax = fig.add_subplot(3,1,1) #(111)
    ax.set_ylabel('track_weight', fontsize = 18)
    plotRes = adrenalRes
    subtable = plotRes.loc[plotRes['ps'] == ps]
    
    for row in range(len(subtable)):
        stws = subtable.iloc[row][1] # getting stws value
        tw = subtable.iloc[row][2] # getting tw value
        
        i = stws_list.index(stws)
        j = tw_list.index(tw)

        thisCount = math.log2(subtable.iloc[row][3])
        halfLen = ((thisCount/maxCount)/2)*.75
        plt.plot([i-halfLen, i+halfLen], [j, j], color='olive', linewidth=4)


    plt.title('prior_strength = %f' %(ps))    
    plt.xticks([-1] + x + [len(x)], [0] + stws_list + [0])
    plt.yticks([-1] + y + [len(y)], [0] + tw_list + [0])
    plt.grid()
    axt = ax.twinx()
    plt.yticks([-1] + y + [len(y)], [0] + tw_list + [0])
    plt.grid()
    #plt.show()

    # >>>>>>>>>> second plot: H1
    ax = fig.add_subplot(3,1,2) #(111)
    ax.set_ylabel('track_weight', fontsize = 18)
    plotRes = h1Res
    subtable = plotRes.loc[plotRes['ps'] == ps]
    for row in range(len(subtable)):
        stws = subtable.iloc[row][1] # getting stws value
        tw = subtable.iloc[row][2] # getting tw value
        
        i = stws_list.index(stws)
        j = tw_list.index(tw)

        thisCount = math.log2(subtable.iloc[row][3])
        halfLen = ((thisCount/maxCount)/2)*.75
        plt.plot([i-halfLen, i+halfLen], [j, j], color='olive', linewidth=4)


    plt.title('prior_strength = %f' %(ps))    
    plt.xticks([-1] + x + [len(x)], [0] + stws_list + [0])
    plt.yticks([-1] + y + [len(y)], [0] + tw_list + [0])
    plt.grid()
    axt = ax.twinx()
    plt.yticks([-1] + y + [len(y)], [0] + tw_list + [0])
    plt.grid()

    # >>>>>>>>>> second plot: DND41
    ax = fig.add_subplot(3,1,3) #(111)
    ax.set_ylabel('track_weight', fontsize = 18)
    plotRes = dndRes
    subtable = plotRes.loc[plotRes['ps'] == ps]
    for row in range(len(subtable)):
        stws = subtable.iloc[row][1] # getting stws value
        tw = subtable.iloc[row][2] # getting tw value
        
        i = stws_list.index(stws)
        j = tw_list.index(tw)

        thisCount = math.log2(subtable.iloc[row][3])
        halfLen = ((thisCount/maxCount)/2)*.75
        plt.plot([i-halfLen, i+halfLen], [j, j], color='olive', linewidth=4)


    plt.title('prior_strength = %f' %(ps))    
    plt.xticks([-1] + x + [len(x)], [0] + stws_list + [0])
    plt.yticks([-1] + y + [len(y)], [0] + tw_list + [0])
    plt.grid()
    axt = ax.twinx()
    plt.yticks([-1] + y + [len(y)], [0] + tw_list + [0])
    plt.grid()

    plt.savefig(plotFolder + '02threeTissues_ps_%.4f.pdf' %(ps))


# old ones
#>>>>>>>>>>>>>>>>>>>>

for plotCount in range(len(ps_list)):

    fig = plt.figure(figsize=(12,10))
    ps = ps_list[plotCount]

    x = list(range(len(stws_list)))
    y = list(range(len(tw_list)))

    # >>>>>>>>>> first plot: adrenal gland
    ax = fig.add_subplot(3,1,1) #(111)
    ax.set_ylabel('track_weight', fontsize = 18)
    plotRes = adrenalRes
    subtable = plotRes.loc[plotRes['ps'] == ps]
    row = 0
    for i in range(len(stws_list)):
        for j in range(len(tw_list)):
            thisCount = math.log2(subtable.iloc[row][3])
            halfLen = ((thisCount/maxCount)/2)*.75
            plt.plot([i-halfLen, i+halfLen], [j, j], color='olive', linewidth=4)
            row +=1

    plt.title('prior_strength = %f' %(ps))    
    plt.xticks([-1] + x + [len(x)], [0] + stws_list + [0])
    plt.yticks([-1] + y + [len(y)], [0] + tw_list + [0])
    plt.grid()
    axt = ax.twinx()
    plt.yticks([-1] + y + [len(y)], [0] + tw_list + [0])
    plt.grid()
    #plt.show()

    # >>>>>>>>>> second plot: H1
    ax1 = fig.add_subplot(3,1, 2) #(111)
    ax1.set_ylabel('track_weight', fontsize = 18)

    plotRes = h1Res
    subtable = plotRes.loc[plotRes['ps'] == ps]
    row = 0
    for i in range(len(stws_list)):
        for j in range(len(tw_list)):
            thisCount = math.log2(subtable.iloc[row][3])
            halfLen = ((thisCount/maxCount)/2)*.75
            plt.plot([i-halfLen, i+halfLen], [j, j], color='olive', linewidth=4)
            row +=1

    plt.xticks([-1] + x + [len(x)], [0] + stws_list + [0])
    plt.yticks([-1] + y + [len(y)], [0] + tw_list + [0])
    plt.grid()
    ax1t = ax1.twinx()
    plt.yticks([-1] + y + [len(y)], [0] + tw_list + [0])
    plt.grid()

    #ax1.set_xlabel('seg_transition_weight_scale', fontsize = 18)
    #plt.show()

    # >>>>>>>>>> second plot: DND41
    ax2 = fig.add_subplot(3, 1, 3) #(111)
    ax2.set_ylabel('track_weight', fontsize = 18)

    plotRes = dndRes
    subtable = plotRes.loc[plotRes['ps'] == ps]
    row = 0
    for i in range(len(stws_list)):
        for j in range(len(tw_list)):
            thisCount = math.log2(subtable.iloc[row][3])
            halfLen = ((thisCount/maxCount)/2)*.75
            plt.plot([i-halfLen, i+halfLen], [j, j], color='olive', linewidth=4)
            row +=1

    plt.xticks([-1] + x + [len(x)], [0] + stws_list + [0])
    plt.yticks([-1] + y + [len(y)], [0] + tw_list + [0])
    plt.grid()
    ax2t = ax2.twinx()
    plt.yticks([-1] + y + [len(y)], [0] + tw_list + [0])
    plt.grid()
    
    ax2.set_xlabel('seg_transition_weight_scale', fontsize = 18)
#    plt.show()

    plt.savefig(plotFolder + 'threeTissues_ps_%.4f.pdf' %(ps))

# tw vs stws
# >>>>>>>>>>>>>>>>>>>>

for plotCount in range(len(ps_list)):

    fig = plt.figure(figsize=(12,10))
    ps = ps_list[plotCount]

    y = list(range(len(stws_list)))
    x = list(range(len(tw_list)))

    # >>>>>>>>>> first plot: adrenal gland
    ax = fig.add_subplot(3,1,1) #(111)
    ax.set_ylabel('seg_transition', fontsize = 18)
    plotRes = adrenalRes
    subtable = plotRes.loc[plotRes['ps'] == ps]
    row = 0
    for i in range(len(stws_list)):
        for j in range(len(tw_list)):
            thisCount = math.log2(subtable.iloc[row][3])
            halfLen = ((thisCount/maxCount)/2)*.75
            plt.plot([i-halfLen, i+halfLen], [j, j], color='olive', linewidth=4)
            row +=1

    plt.title('prior_strength = %f' %(ps))    
    plt.xticks([-1] + x + [len(x)], [0] + stws_list + [0])
    plt.yticks([-1] + y + [len(y)], [0] + tw_list + [0])
    plt.grid()
    axt = ax.twinx()
    plt.yticks([-1] + y + [len(y)], [0] + tw_list + [0])
    plt.grid()
    #plt.show()

    # >>>>>>>>>> second plot: H1
    ax1 = fig.add_subplot(3,1, 2) #(111)
    ax1.set_ylabel('track_weight', fontsize = 18)

    plotRes = h1Res
    subtable = plotRes.loc[plotRes['ps'] == ps]
    row = 0
    for i in range(len(stws_list)):
        for j in range(len(tw_list)):
            thisCount = math.log2(subtable.iloc[row][3])
            halfLen = ((thisCount/maxCount)/2)*.75
            plt.plot([i-halfLen, i+halfLen], [j, j], color='olive', linewidth=4)
            row +=1

    plt.xticks([-1] + x + [len(x)], [0] + stws_list + [0])
    plt.yticks([-1] + y + [len(y)], [0] + tw_list + [0])
    plt.grid()
    ax1t = ax1.twinx()
    plt.yticks([-1] + y + [len(y)], [0] + tw_list + [0])
    plt.grid()

    #ax1.set_xlabel('seg_transition_weight_scale', fontsize = 18)
    #plt.show()

    # >>>>>>>>>> second plot: DND41
    ax2 = fig.add_subplot(3, 1, 3) #(111)
    ax2.set_ylabel('track_weight', fontsize = 18)

    plotRes = dndRes
    subtable = plotRes.loc[plotRes['ps'] == ps]
    row = 0
    for i in range(len(stws_list)):
        for j in range(len(tw_list)):
            thisCount = math.log2(subtable.iloc[row][3])
            halfLen = ((thisCount/maxCount)/2)*.75
            plt.plot([i-halfLen, i+halfLen], [j, j], color='olive', linewidth=4)
            row +=1

    plt.xticks([-1] + x + [len(x)], [0] + stws_list + [0])
    plt.yticks([-1] + y + [len(y)], [0] + tw_list + [0])
    plt.grid()
    ax2t = ax2.twinx()
    plt.yticks([-1] + y + [len(y)], [0] + tw_list + [0])
    plt.grid()
    
    ax2.set_xlabel('seg_transition_weight_scale', fontsize = 18)
#    plt.show()

    plt.savefig(plotFolder + 'threeTissues_ps_%.4f.pdf' %(ps))    


# draft section
####################
book = results[results['ps'] ==0]

# axis labels
ax.set_xlabel('kbp distance from the SNPs', fontsize = 16)
ax.set_ylabel('SNPs', fontsize = 16)
plt.title('count of overlapping encyclopedia segments', fontsize = 20)

plt.savefig(figureFolder + 'heatmapSegCountShuff01_sorted.pdf')
plt.show()

plt.savefig(figureFolder + 'heatmapSegCountMain.pdf')
book = sum(segCount)
plt.plot(range(200), book, '-', color='black')
plt.show()

# 4. line plot showing the density of the regions as we get further away (curve fitting here)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# the main dots/line for 1Mbp
fig = plt.figure(figsize=(14,8))
ax = fig.add_subplot(111)
x,y = range(1000), book
s = interpolate.UnivariateSpline(x, y, s = 5e9, k = 5)
xs = np.linspace(0, 999, 3000)
ys = s(xs)
ax.plot(x, y, color='olive', alpha=1)
ax.plot(xs, ys, lw=3)
ax.set_xlabel('kbp distance from the SNPs', fontsize = 18)
ax.set_ylabel('overlapping segments', fontsize = 18)
plt.savefig(figureFolder + 'encycMain1Mbp_COPD139.pdf')
plt.show()

# the main dots/line for 200kbp
x,y = range(200), book
s = interpolate.UnivariateSpline(x, y, s = 5e9, k = 5)
xs = np.linspace(0, 199, 600)
ys = s(xs)

s2 = interpolate.UnivariateSpline(x, meanSegCount, s = 5e9, k = 5)
y2 = s2(xs)

plt.plot(x, meanSegCount)

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)

for i in range(9):
    plt.plot(x, segCountPerDis[:,i], c = 'k', alpha = .3)

shuffleLines, = ax.plot(x, segCountPerDis[:,9], c = 'k', alpha = .3)
shuffleMean, = ax.plot(xs, y2, lw = 3)
mainDots, = ax.plot(x, y, '.')
mainLine, = ax.plot(xs, ys, lw = 3)

ax.legend((mainLine, mainDots, shuffleMean, shuffleLines), ('encyclopedia fitted','encyclopedia overlap','shuffle fitted','shuffled overlap'))

# axis labels
ax.set_xlabel('kbp distance from the SNPs', fontsize = 18)
ax.set_ylabel('overlapping segments', fontsize = 18)

plt.savefig(figureFolder + 'shuffleVSmain02.pdf')
plt.show()

# >>>>>> the count to go with the heatmap


# do the plots
from matplotlib import colors
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)
cmap = colors.ListedColormap(['gray', 'red', 'blue', 'yellow', 'white'])
bounds = [-.5, .5, 1.5, 2.5, 3.4, 4.5]
norm = colors.BoundaryNorm(bounds, cmap.N)
heatmap = plt.pcolor(segCountMain[sindexMain,], cmap=cmap, norm=norm)
heatmap = plt.pcolor(segCountsh01[sindexShuff,], cmap=cmap, norm=norm)
plt.colorbar(heatmap, ticks=[0, 1, 2, 3, 4])



dirName = dirList[0]
nameSplit = dirName.split('_')
print(nameSplit)

index = int(nameSplit[0][9:])
ps = float(nameSplit[4])
stws = float(nameSplit[5])
tw = float(nameSplit[6])

segSizeFile = '%s/%s/length_distribution/segmentSizes.tab' %(resultFolder, dirName)

####################






