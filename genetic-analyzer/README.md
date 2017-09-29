<img src="../foundry-logo.png" height="120px"/>

# Characterize and debug genetic circuits

This directory contains scripts and example files for the characterization and debugging of genetic circuits using RNA-seq data. For further information see Gorochowski _et al._ "Genetic circuit characterization and debugging using RNA-seq", _Molecular Systems Biology_, 2017 (_in press_).

## Overview

The code within this directory implements an analysis pipeline for characterizing genetic circuits using RNA-seq data. The `bin` directory contains the core set of scripts that perform the analyses. An example analysis pipeline and associated data are contained in the `circuit_example` directory. The pipeline requires the following software to be available in the local path:
- Python version 2.7.9
- R version 3.2.1
- edgeR version 3.8.6
- BWA version 0.7.4
- Samtools version 1.4
- HTSeq version 0.9.1

## Analysis Pipeline 

The pipeline is composed of a number of steps that call functions contained within the `genetic_analyzer.py` script and several separate R scripts. To simplify this process, additional wrapper scripts are also provided that can be called from the command-line and take arguments to specify the data files to use and parameters for the analyses to perform. All of these scripts are contained within the `bin` directory.

To coordinate the calls needed to run the entire pipeline, these wrappers are executed by a set of numbered shell scripts that should be executed in order. These shell scripts are customised for each specific experiment and have the following functions (a set of examples are provided in the `circuit_example` directory):
- `00_setup.sh` – Create all required directories to hold the results, logs and for temporary storage.
- `01_map_reads.sh` – Map the raw RNA-seq reads to the reference sequences (using Samtools and BWA).
- `02_count_reads.sh` – Count the reads per gene in the host chromosome and any circuit plasmids (using HTSeq).
- `03_fragment_distributions.sh` – Calculate the fragment size distribution.
- `04_read_analysis.sh` – Collate mapped read information from all samples into a single matrix and calculate TMM normalization factors (using R and edgeR)
- `05_transcription_profiles.sh` – Generate the normalized transcription profiles.
- `06_de_analysis.sh` – Perform differential gene expression analysis (using R and edgeR).
- `07_part_analysis.sh` – Calculate the performance of all promoters, ribozymes and terminators in the circuit.
- `08_promoter_fitting.sh` – Fit the performance data of the input and output promoters to each NOT-gate to a hill function describing its response function.
- `09_clean_up` - Remove any temporary files and logs.

Data required by the pipeline is recommended to be contained with a `data` directory that includes a `settings.txt` tab-delimited file. The first row labeled "None" gives the location of where the results should be stored. All subsequent rows correspond to a particular sample and state for the circuit and the relevant location of data files needed for the analysis. To simplify the storage of data files we recommend the following subdirectories are used:
- `bed/` - BED files defining the sequence ranges to create the transcription profile.
- `fasta/` - FASTA files containing sequences of the host and circuit.
- `gff/` - GFF files giving annotations of the host and circuit (see Circuit Description below for further details).
- `fastq/` - FASTQ files containing the raw sequencing reads.

## Circuit Description

For the pipeline to be able to understand what genetic parts are present, their location, and other information specific to their function (e.g. cut site location of a ribozyme), we use a number of custom GFF annotation types to store this information: 
- transcript – this defines the start and end of a transcript in the design
- promoter – a promoter part
- ribozyme – a ribozyme part. This has the attribute "cut_site" that defines the base at which cleavage occurs.
- terminator – a terminator part
- promoter_unit – an entire promoter unit. This spans from the start of the first promoter to the cur site of the ribozyme downstream (as this is where we calculate the combined expression). A number of attributes are then used to tell the pipeline what this promoter unit contains and how the internal promoters are regulated (an example is provided in the `gff` directory). Specifically:
  - promoter_names – list of internal promoter names separated by commas. 
  - promoter_types – list of the internal promoter types separated by commas. Types can be "induced" or "repressed".
  - chrom_inputs – list of chromosomes for the input promoter units with names separated by commas.
  - promoter_unit_inputs – These are the promoter units that feed as input to each of the internal promoters. For example, if one of the internal promoters was P_Lit, then the input promoter unit would be the one that expresses LitR. The ordering of the inputs is the same as the ordering used for all the other attributes and the list should be separated by commas. If the promoter is induced, then this is where the truth table for the various samples is entered. This takes the format: "sampleA1>0:sampleA2>1", where >0 and >1 mean that for that sample the input is ON or OFF respectively. Note that all elements should have a "Name" attribute that must be unique.
