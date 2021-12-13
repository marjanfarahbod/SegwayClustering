# to get trackIDs-histone mark mappings from a list of samples, the mappings should come in the .txt, with first column IDs and second column histone marks or identifiers.

from urllib.request import urlopen
import json

# fetch the list of samples
sampleFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/'
contents = os.listdir(sampleFolder)

outputFile = sampleFolder + 'trackname_assay.txt'

with open(outputFile, 'w') as output:
# for a given sample, fetch the feature: they are in the header of the signal_distribution.csv file.
    for item in contents:
        if os.path.isdir(sampleFolder + item):
            csvfile = sampleFolder + item + '/Oct22/call-segtools/glob-2f0453ff805b02858906905b93a1b15d/signal_distribution.csv'

            f = open(csvfile, 'r')
            header = f.readline()
            IDs = header.strip().split(',')
            IDs.pop(0)

            for i, id in enumerate(IDs):
                IDs[i] = id.replace('"', '')
    
            # for each ID, get the jason file
            # get the mark, print the ID and the track
            for id in IDs:
                url = "https://www.encodeproject.org/files/%s/?format=json" % (id)
                response = urlopen(url)
        
                data_json = json.loads(response.read())
                if 'target' in data_json:
                    label = data_json['target']['label']
                else:
                    label = data_json['assay_title']
            

                output.write('%s\t%s\n' %(id,label))
                # get the error

