# This is to get the data files (hard coded for bigwig) from a given list of ENCODE file IDs. For this file, I intend to get the IDs from the old segway runs, so I use the file signal_distribution.csv to fetch the IDs. It needs to be modified to read IDs from an input file

from urllib
import json



"""""""""""""""""""""""
To get the list of file IDs that we want to download in IDs

"""""""""""""""""""""""

sampleFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/adrenal-gland/'
csvfile = sampleFolder + 'Oct22/call-segtools/glob-2f0453ff805b02858906905b93a1b15d/signal_distribution.csv'

f = open(csvfile, 'r')
header = f.readline()
IDs = header.strip().split(',')
IDs.pop(0)

"""""""""""""""""""""""
running the wget command

"""""""""""""""""""""""

for id in (IDs):
    thisID = id.replace('"', '')
    command = 'wget https://www.encodeproject.org/files/' + thisID + '/@@download/' + thisID + '.bigWig'
    
            # for each ID, make the gwet command
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
