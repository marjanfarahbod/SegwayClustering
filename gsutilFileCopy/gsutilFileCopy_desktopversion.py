
import os

# inputFile = 'fileList.txt'

inputFile = 'fileList_filteredForDesktop.txt'
with open(inputFile, 'r') as fileList:
    line = fileList.readline()
    
    for line in fileList:
        #    line = fileList.readline()
        stripLine = line.strip()
        #    print(stripLine)
        resultPath = '/'.join(stripLine.split('/')[4:-1])
        #    print(resultPath)
        os.makedirs(resultPath, exist_ok=True)
        os.system('/Users/marjanfarahbod/google-cloud-sdk/bin/gsutil cp %s %s' %(stripLine, resultPath))
        #os.system('gsutil cp %s %s' %(stripLine, resultPath))


