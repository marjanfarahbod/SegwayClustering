# for each sample, get the gmtk file,

# get the track mappings,

inputFileName = 'all_annInfo_list.pkl'
inputFile = dataFolder + dataSubFolder + outputFileName
with open(inputFile, 'rb') as f:
    ann_info_list = pickle.load(f)


runID_file = '/Users/marjanfarahbod/Documents/projects/segwayLabeling/data/testBatch105/' + 'runID_accession_map_105run.pkl'

with open(runID_file, 'rb') as pickledFile:
    runID_accession = pickle.load(pickledFile)

accession_runID = {}
for runID in list(runID_accession.keys()):
    ac = runID_accession[runID]
    accession_runID[ac] = runID

#########################################    
index = 79

for index in range(105):
    ann = ann_info_list[index]
    annAccession = ann['accession']
    print(annAccession)

    runID = accession_runID[annAccession]

    segtools_add = dataFolder + 'testBatch105/all_segtools/' + runID
    print(segtools_add)

    gmtk_file = segtools_add + '/glob-03b7332b8fdb9a1ca33a23093d5878d5' + '/gmtk_parameters.stats.csv'
    df = pd.read_csv(gmtk_file)

    df = df[df.columns[1:]]

    track_accession = list(df)
    labels = list(df.index)

    # get column list
    track_names = []
    for track in track_accession:
        track_names.append(ann['track_assay_map'][track])

    # get row list
    ordered_labels = list(df.index)
    segway_labels = list(ann['segway_anns']['clusters'].keys())
    for label in segway_labels:
        number = int(label.split('_')[0])
        ordered_labels[number] = label
    

    df.columns = track_names
    df.index = ordered_labels

    kado = df.transpose()

    fig,axs = plt.subplots(1, figsize=(4,3))
    sns.heatmap(kado, cmap='vlag')
    sns.set(font_scale=.7)
    plt.tight_layout()
    plt.title(annAccession)

    figFile = plotFolder + annAccession + '/GMTK_emmission.pdf'
    print(figFile)
    plt.savefig(figFile)
    plt.close('all')


plt.show()

plots_add = plotFolder + annAccession
print(plots_add)


# read .csv to a dataframe


