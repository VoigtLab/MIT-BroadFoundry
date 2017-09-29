<img src="../foundry-logo.png" height="120px"/>

# Dialout pooled composite parts

This directory contains scripts and example files for the dialout process of pooled composite parts. For further information see [Woodruff _et al._ "Registry in a tube: multiplexed pools of retrievable parts for genetic design space exploration", _Nucleic Acids Research_ **45**, 1553-1565, 2017](http://nar.oxfordjournals.org/content/early/2016/12/21/nar.gkw1226.full).

## 1. Linker design

To choose a set of compatible and orthogonal linkers (scars) the `linker_design.py` script it used. An example of it's use is:

    python ./bin/linker_design.py 4 2 -1 2 ./example/linker-design/seed_scars.txt None out_scars.txt

This will generate scars of length 4 bp with a maximum homology of 2 bp. It will enumerate all possible linkers and use `seed_scars.txt` as a set of initial selected scars before outputing the results to `out_scars.txt`. To see full options for the script use `python ./bin/linker_design.py -h`.

## 2. Dialout barcodes for composite parts

Once the pool has been sequenced the raw FASTQ file is interrogated to assess the designs that are present and the associated barcodes for each. This functionality is contained within the `dialout_barcodes.py` script. This takes as input regular expressions of each design with variable `N` regions of fixed lengths that correspond to potential barcode regions (see the `./example/regular-exprs/` folder for examples), two FASTQ files containing the raw sequencing data for each paired read (order of reads in each file must match), the length of the primer regions, the indexes of the two barcodes (for retrieval), plus the indexes of other variable regions to extract (given as comma separated list of indexes with no spaces), and finally the output path for generated file to be placed (see script comments for further information).

    python ./bin/dialout_barcodes.py ./example/regular-exprs/AmeR_2NOR_regexs.txt seq_data_read_1.fastq seq_data_read_2.fastq 20 20 1 2 3 ./

After execution, four files are generated:

- `dialout_designs.csv` - information about the barcodes found for each of the possible designs.

- `dialout_barcodes.csv` - details of all barcoded designs. Not all will be suitable for retrieval, but data is useful for analyzing barcoding process and troubleshooting issues.

- `dialout_design_unique_barcodes.csv` - details of all the uniquely barcoded designs that are suitable for retrieval.

- `dialout_summary.txt` - a summary of the dialout process including information regarding the run-time, total reads, matched reads, uniquely barcoded designs, etc.

## 3. Analyze the dialout process

Analysis of the dialout process can be performed to give statistics on the success of barcoding designs containing different combinations of parts. The `dialout_part_analysis.py` script analyzes the `dialout_designs.csv` file generated in the previous step.

    python ./bin/dialout_part_analysis.py dialout_barcodes.csv output_analysis.csv

## 4. Primer design for composite part barcodes

The final step in the process is to automatically generate suitable primers for the barcodes found from the sequencing data. We recommend the optimized Primer3 settings in the `dialout_primer3_settings.txt` file. This can be used as input for the `dialout_primer_design.py` script with the input barcode file generated from the dialout process, the indexes of the forward and reverse barcode to use and any short sequences that will be present in the backbone of the circularized vector. An example of calling this script is as follows, with the generated primers contained in the `output_primers.csv` file:

    python ./bin/dialout_primer_design.py ./bin/dialout_primer3_settings.txt ./example/input_barcodes.csv 3 2 output_primers.csv GC GC
