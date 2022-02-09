## Content
Utility and pipeline code for Segway classifier labels.

## Items

### QC with transcriptomic data

For this task we need three input files:

- The gene coordinates. We use the GTF file used for running some segtool commands for this ENCODE runs: https://www.encodeproject.org/files/gencode.v29.primary_assembly.annotation_UCSC_names/ genomic coordinates used for Segtools
- The RNA-seq data for the sample(s) - each sample has one RNAseq file
- The annotations with class labels (the initial sample with the transcript data name:H1 ID:ENCSR938GXK) - each sample has one .bed file

### Parsing and exploring the GTF file


#### Products:
- geneTypesFromGTF.pkl: list of gene types and their gene counts in the file 
- coorsFile.pkl: a dataframe with gene coordinates


### RNAseq info and pre-processing

Which file to download from protal

#### code 

#### Products

### Transcript data comparison

code:

#### Products


### Getting the files I want for a list of samples from GCP
TODO (so far I have downloaded them using command line) 

### Getting the trackname_assay.txt for a list of samples
For the Oct2021 test runs, I got it for each sample from the signal_distribution.csv, based on the headers I got the track IDs and from json files I got the assay type - see code in general/trackIDmappings.py

### Given the sample ID, obtain .bed files from the GCP, merge the class and cluster labels
There are two bed files in the GCP storage:

- gs://segway-testing/segway/[GCP dir ID]/call-segway_annotate/segway.bed.gz (cluster labels)
- gs://segway-testing/segway/[GCP dir ID]/call-recolor_bed/recolored.bed.gz (class labels)

The GCP ID was given to us with experiment IDs along with tissue names. In local I use tissue names rather than the IDs to name folders for different samples. While the IDs stay unique for tissue names, the <dir> in the GCP changes for each segway run. In the code I had to do some matching and the matching part might vary for different runs based on the meta file format.
the code is in getbedsfromGCP.py




