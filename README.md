## Content
Utility and pipeline code for Segway classifier labels.

## Items

### given the sample ID, obtain .bed files from the GCP
Last update: Jan 10 - 2022
There are two bed files in the GCP storage:
gs://segway-testing/segway/[GCP ID]/call-segway_annotate/segway.bed.gz
gs://segway-testing/segway/[GCP ID]/call-recolor_bed/recolored.bed.gz

The GCP ID was given to us with experiment IDs. I used a bedtools command to merge the cluster/classifier labels in one file.


