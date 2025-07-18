import shutil

#runID_file = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/testBatch105/' + 'runID_accession_map_105run.pkl'
runID_file = '/home/mfarahbo/scratch/runID_accession_map_105run.pkl'

with open(runID_file, 'rb') as pickledFile:
    id_info = pickle.load(pickledFile)

runIDList = list(id_info.keys())

for runID in runIDList:
    source = '/home/mfarahbo/projects/rrg-maxwl/mfarahbo/encode-segway/' + runID + '/call-segtools/'
    dest = '/home/mfarahbo/scratch/all_segtools/' + runID

    shutil.copy2(source, dest)

