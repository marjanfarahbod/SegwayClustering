# just do that the name says, this is for the 112 batch runs

### local add
dataFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/the112batch/'

sampleFileLists = 'SegwayFileSets.tsv'
# for each sample, make the folder with the accession name
inputFile = dataFolder + sampleFileLists
segAccessionList = []
with open(inputFile, 'r') as f:
    header = f.readline()
    for line in f:
        fields = line.strip().split('\t')
        print(line)

        # with fields 0 make the folder
        segwayAccession = fields[0]
        segAccessionList.append(segwayAccession)

outputFile = dataFolder + 'accessionList.pkl'
with open(outputFile, 'wb') as f:
    pickle.dump(segAccessionList,f)
        
