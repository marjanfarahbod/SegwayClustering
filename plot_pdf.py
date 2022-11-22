import sys
import os
import pickle

# load the plot address for each of the samples.

# a function for, given the index, give me the length of the segment


index = 98
index = 79 2wjp
index = 56 DZH
index = 77 CCR
index = 12

inputFileName = 'all_annInfo_list.pkl'
inputFile = dataFolder + dataSubFolder + inputFileName
with open(inputFile, 'rb') as f:
    ann_info_list = pickle.load(f)

    
runID_file = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/testBatch105/' + 'runID_accession_map_105run.pkl'

with open(runID_file, 'rb') as pickledFile:
    runID_accession = pickle.load(pickledFile)
accession_runID = {}
for runID in list(runID_accession.keys()):
    ac = runID_accession[runID]
    accession_runID[ac] = runID

plotFolder = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/plots/testBatch105/'


html_file = plotFolder + 'plot_directory.html'
init_text = '''
<html>
  <head>
    <title>My First Web Page!</title>
  </head>
  <body>
'''
end_text = '''
  </body>
</html>
'''

with open(html_file, 'w') as output:
    output.write(init_text)

    for index in range(105):

        ann = ann_info_list[index]
        annAccession = ann['accession']
        print(annAccession)
    
        runID = accession_runID[annAccession]
        segtools_folder_add = dataFolder + 'testBatch105/all_segtools/' + runID + '/'


        sample_plot_folder = plotFolder + annAccession 

        plot_file = sample_plot_folder + '/the_panel_03.pdf'
        text_sample_plot = str(index) + '___' + annAccession + '__sample'
        mid_text = '<p><a href="%s">%s</a></p>' %(plot_file, text_sample_plot)
        output.write(mid_text)

        plot_file = sample_plot_folder + '/GMTK_emmission.pdf'
        text_sample_plot = str(index) + '___' + annAccession + '__GMTK'
        mid_text = '<p><a href="%s">%s</a></p>' %(plot_file, text_sample_plot)
        output.write(mid_text)

        plot_file = sample_plot_folder + '/expression_transcript_02.pdf'
        if os.path.isfile(plot_file):
            text_sample_plot = str(index) + '___' + annAccession + '__transcript'
            mid_text = '<p><a href="%s">%s</a></p>' %(plot_file, text_sample_plot)
            output.write(mid_text)

        plot_file = sample_plot_folder + '/length_dist_hist.pdf'
        if os.path.isfile(plot_file):
            text_sample_plot = str(index) + '___' + annAccession + '__length_dist'
            mid_text = '<p><a href="%s">%s</a></p>' %(plot_file, text_sample_plot)
            output.write(mid_text)

        plot_file = sample_plot_folder + '/genome_coverage.pdf'
        if os.path.isfile(plot_file):
            text_sample_plot = str(index) + '___' + annAccession + '__genome_coverage'
            mid_text = '<p><a href="%s">%s</a></p>' %(plot_file, text_sample_plot)
            output.write(mid_text)


        #plot_file = segtools_folder_add + 'glob-a287da44f32bf6a3fc6d7c51c52ddafa/' + 'length_distribution.pdf'
        #text_sample_plot = str(index) + '___' + annAccession + '__length_dist'
        #mid_text = '<p><a href="%s">%s</a></p>' %(plot_file, text_sample_plot)
        #output.write(mid_text)
        

    # printing the html file
    output.write(end_text)
    
#########################################        
#########################################    
#########################################    
index = 26
label = '7_RegPermissive'

ann = ann_info_list[index]
annAccession = ann['accession']
print(annAccession)

runID = accession_runID[annAccession]

segtools_add = dataFolder + 'testBatch105/all_segtools/' + runID
print(segtools_add)

plots_add = plotFolder + annAccession
print(plots_add)

total_bp = 0
cluster_list = ann_info_list[index]['segway_anns']['clusters']
for cluster in cluster_list:
    total_bp += ann_info_list[index]['segway_anns']['clusters'][cluster].bp_count

print('bp_ratio')
print(ann_info_list[index]['segway_anns']['clusters'][label].bp_count/total_bp)


