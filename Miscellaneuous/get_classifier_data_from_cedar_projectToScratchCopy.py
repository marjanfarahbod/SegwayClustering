
import shutil

class_files = 'classifier_data_file_list.txt'
# class_files = dataFolder + 'classifier_data_file_list.txt'

c = 0
with open(class_files, 'r') as inputFile:
    for line in inputFile:

        source = '/home/mfarahbo/projects/rrg-maxwl/mfarahbo/encode-segway' + line.strip()
        dest = '/home/mfarahbo/scratch/' + line.split('/')[1] + '_classifier_data.tab'
        
        shutil.copy2(source, dest)

