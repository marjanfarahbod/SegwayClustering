# for batch run on cedar on three parameters: track-weight, STWS, prior-strength

import math
import os
import shutil

# FIX FRACTION

# where the genomedata is
genomedataFile =  'data/DND-41.genomedata'


# parameter setting
# parameterIndex = int(os.environ['SLURM_ARRAY_TASK_ID'])
parameterIndex = 602

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
# tempTrainFolder = '%s/tempTrain' %(outputFolder)
outputTrainFolder = '%s/train' %(outputFolder)
trainCommand = 'segway train --include-coords=data/encodePilotRegions.hg19.bed --num-labels=15 --resolution=100 --minibatch-fraction=0.01 --num-instances=10 --prior-strength=%.1f --segtransition-weight-scale=%.3f --track-weight=%.4f --ruler-scale=100 --max-train-rounds=25 %s %s' %(ps, stws, tw, genomedataFile, outputTrainFolder)

os.system(trainCommand)

# clean up
# os.system('mkdir %s' %trainFolder)
# shutil.rmtree('%s/output' %(outputTrainFolder))
# shutil.rmtree('%s/cmdline' %(outputTrainFolder))
# shutil.rmtree('%s/accumulators' %(outputTrainFolder))
# shutil.rmtree('%s/auxiliary' %(outputTrainFolder))
# shutil.rmtree('%s/intermediate' %(outputTrainFolder))
# shutil.rmtree('%s/observations' %(outputTrainFolder))

# trainFolder = '%s/train' %(outputFolder)
# os.system('mkdir %s' %trainFolder)
# os.rmdir('%s/output' %(outputTrainFolder))
# os.rmdir('%s/cmdline' %(outputTrainFolder))
# os.rmdir('%s/accumulators' %(outputTrainFolder))
# os.rmdir('%s/auxiliary' %(outputTrainFolder))
# os.rmdir('%s/intermediate' %(outputTrainFolder))
# os.rmdir('%s/observations' %(outputTrainFolder))

# segway posterior command
outputPostFolder = '%s/posterior' %(outputFolder)
postCommand = 'segway posterior --include-coords=data/encodePilotRegions.hg19.bed data/DND-41.genomedata %s %s' %(outputTrainFolder, outputPostFolder)

os.system(postCommand)

# unzip the .bed file command
os.system('gunzip %s/segway.bed.gz %s/segway.bed' %(outputPostFolder, outputPostFolder))

# cd to the segwayRun folder (see above) command
os.chdir(outputFolder)

# segtools-length-distribution command
os.system('segtools-length-distribution posterior/segway.bed')

# segtools-gmtk-parameters command
os.system('segtools-gmtk-parameters train/params/params.params')


