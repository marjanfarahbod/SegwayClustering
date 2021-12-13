# for batch run on cedar on three parameters: track-weight, STWS, prior-strength
# changing the code so I can delete files after each run

# final run, all is working
# select genomedata
# test run, fraction .1, 

import math
import os
import shutil

# where the genomedata is
#genomedataFile =  'data/DND41.genomedata'
genomedataFile =  'data/H1.genomedata'

# parameter setting
parameterIndex = int(os.environ['SLURM_ARRAY_TASK_ID'])
#parameterIndex = 602

ps_list = [0, .1, .5, 1, 10, 100, 1000]
# seg transition weight scale
stws_list = [0, .001, .01, .5, 1, 10, 50, 100, 500, 1000, 5000, 10000]
# track weight
tw_list = [.0005, .001, .002, .005, .008, .01, .02, .05, .08, .1, .5, 1, 10, 50, 100]

den1 = len(tw_list) * len(stws_list)
den2 = len(tw_list)

ps_ind = int(math.modf(parameterIndex/den1)[1])
temp = parameterIndex%den1
stws_ind = int(math.modf(temp/den2)[1])
tw_ind = temp%den2

ps = ps_list[ps_ind]
stws = stws_list[stws_ind]
tw = tw_list[tw_ind]

# make the output folder for the run: segwayRun[runID]_ps_stws_tw_[psvalue]_[stwsvalue]_[twvalue]
outputFolder = 'segwayRun%d_ps_stws_tw_%.1f_%.3f_%.4f' %(parameterIndex, ps, stws, tw)
os.system('mkdir %s' %(outputFolder))

# segway train command
trainFolder = '%s/train' %(outputFolder)
trainCommand = 'PYTHONDONTWRITEBYTECODE= segway train --include-coords=data/encodePilotRegions.hg19.bed --num-labels=15 --resolution=100 --minibatch-fraction=0.1 --num-instances=10 --prior-strength=%.1f --segtransition-weight-scale=%.3f --track-weight=%.4f --ruler-scale=100 --max-train-rounds=25 %s %s' %(ps, stws, tw, genomedataFile, trainFolder)

os.system(trainCommand)

# segway posterior command
postFolder = '%s/posterior' %(outputFolder)
postCommand = 'PYTHONDONTWRITEBYTECODE= segway posterior --include-coords=data/encodePilotRegions.hg19.bed %s %s %s' %(genomedataFile, trainFolder, postFolder)

os.system(postCommand)

# # collect results
# os.system('cp %s/*.gz %s' %(postFolder, outputFolder))
# os.system('cp %s/log/segway.sh %s/segway.sh' %(trainFolder, outputFolder))
# os.system('cp %s/params/params.params %s/params.params' %(trainFolder, outputFolder))

# collect results
os.system('cp %s/segway.bed.gz %s/' %(postFolder, outputFolder))
os.system('cp %s/log/segway.sh %s/' %(trainFolder, outputFolder))
os.system('cp %s/segway.str %s/' %(trainFolder, outputFolder))
os.system('cp %s/params/params.params %s' %(trainFolder, outputFolder))

# # # # doing the clean up
shutil.rmtree(trainFolder)
shutil.rmtree(postFolder)

# # # cd to the segwayRun folder (see above) command
os.chdir(outputFolder)


# # # unzip the .bed file command
os.system('gunzip segway.bed.gz')


# # # segtools-length-distribution command
os.system('segtools-length-distribution segway.bed')

# # # segtools-gmtk-parameters command
os.system('segtools-gmtk-parameters params.params')
