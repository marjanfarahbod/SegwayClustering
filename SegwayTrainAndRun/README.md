# SegwayTrainAndRun

This directory contains scripts to automate the process of running **Segway** and **Segtools** for genome annotation. The workflow is designed to be executed on a high-performance computing (HPC) cluster using the **SLURM** job scheduler. The scripts are written in Python and Bash.

For more details on Segway itself, refer to the official documentation:  
[Segway Documentation](https://www.pmgenomics.ca/hoffmanlab/proj/segway/doc/0.2.1/)

## Workflow

The overall workflow consists of the following steps:

### Data Preparation
- An accession list of samples is loaded from a pickle file (`accessionList.pkl`).
- For each sample, the corresponding **bedGraph** files are used to create a **genomedata** file using the `genomedata-load` command.  
  This step is handled by `wrapperForGenomeData.py`.

### Segway Execution
- The `wrapperForSegway.py` script runs Segway on each sample's genomedata file.  
  - It first runs `segway train` to train a model.  
  - Then, it runs `segway posterior` to generate a genome annotation in the form of a `.bed.gz` file.  
- The script is designed to be run as a **job array** on a cluster, where each job processes one sample.

### Segtools Analysis
- The `wrapperForSegtools.py` script takes the output of Segway and runs various Segtools analyses on it.  
- These analyses include:  
  - `segtools-signal-distribution`  
  - `segtools-aggregation`  
  - `segtools-length-distribution`  
  - `segtools-gmtk-parameters`  
- This script is also designed to be run as a job array.

### Parameter Tuning
- The `batch_run05.py` script is used for **hyperparameter tuning** of the Segway training process.  
- It iterates through different values for `prior-strength`, `segtransition-weight-scale`, and `track-weight`.

## How to Use

The scripts in this directory are designed to be run on the **Cedar HPC cluster**. The following steps outline the general usage:

### Prerequisites
- Ensure that **Segway**, **Segtools**, and **Python** are installed and available in your environment.  
- The scripts expect a specific directory structure and the presence of an `accessionList.pkl` file containing the list of samples to be processed.

### Data Preparation
- Run the `wrapperForGenomeData.py` script as a **SLURM job array** to create the genomedata files for each sample.

### Running Segway
- Run the `wrapperForSegway.py` script as a **SLURM job array** to run Segway on each sample.

### Running Segtools
- Run the `wrapperForSegtools.py` script as a **SLURM job array** to perform Segtools analysis on the Segway output.

### Parameter Tuning
- Use the `batch_run05.py` script to perform hyperparameter tuning for Segway.  
- This script should also be run as a **SLURM job array**.

## Scripts

- `wrapperForGenomeData.py`: A wrapper script to create genomedata files from bedGraph files.  
- `wrapperForSegway.py`: A wrapper script to run the Segway training and posterior commands.  
- `wrapperForSegtools.py`: A wrapper script to run various Segtools analyses.  
- `batch_run05.py`: A script for hyperparameter tuning of Segway.  
- `fileDownloadFromENCODE.py`: A utility script to download files from ENCODE.  
- `gettingBedGraph_bash.sh`, `gettingGenomedata_bash.sh`, `gettingSegway_bash.sh`, `gettingSegtools_bash.sh`: Bash scripts for performing various steps of the workflow.  
- Other scripts in this directory are helper scripts and utilities for tasks such as file management, data selection, and plotting.
