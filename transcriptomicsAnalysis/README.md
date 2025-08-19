# Transcriptomics analysis for annotations, comparing ChromHHM and Segway 

## Overview

This project contains Python scripts comparing ChromHMM and Segway annotations based on transcriptomics data.

All data are downloaded from ENCODE portal. In the code, local addresses for  Segway annotation, ChromHMM annotation and Transcriptomic data for each sample are included in annMeta, in metaInfo.pkl. 

The primary goal is to assess the relationship between chromatin states (for both Segway and ChromHMM) and gene expression levels (transcriptomic data). 

The analysis includes:
- Quality control (QC)  
- Regression modeling  
- Plotting and visualization  

---

## Scripts

Here is a breakdown of the scripts in this directory and their functions:

- `chromHMMSegway_comparison_transcription_01.py`: For comparing the coverage of ChromHMM and Segway labels in genomic regions, based on the transcriptomic data. Calls functions from SegwayInterpretation/transcription_overlap.py
- `chromHMMSegway_comparison_transcription.py`: An alternate version of the above code, with the loop instead of the function call
- `chromHMMSegway_comparison_transcription_regression.py`: Performs regression analysis to model the relationship between chromatin features and gene expression, for both Segway and ChromHMM
- `chromHMMSegway_comparison_transcription_regression_otherSamples.py`: Extension of the regression analysis, where predictions were based on annotations from different samples. 
- `QC_transcriptionComparison_main.py`: old script for running various quality control for transcriptomic data analyses - main script is in SegwayInterpretation/QC_transcriptionComparision_03.py
- `QC_transcriptionComparison_01.py`: old script for running various quality control for transcriptomic data analyses - main script is in SegwayInterpretation/QC_transcriptionComparision_03.py
- `QC_transcriptionComparison_02.py`: old script for running various quality control for transcriptomic data analyses - main script is in SegwayInterpretation/QC_transcriptionComparision_03.py
- `transcript_plot.py`: Generating plots and visualizations from the analysis results.  

---

## Dependencies

These scripts require **Python 3** and the following libraries:

- `pandas`  
- `numpy`  
- `matplotlib`  
- `seaborn`  
- `scipy`  
- `scikit-learn`  
- `statsmodels`

Note: You need to modify the scripts to specify hte correct input file paths. 

