<img src="../assets/foundry-logo.png" height="120px"/>

# Dialout composite part barcodes



## 1. Make the pools of composite parts


## 2. Dialout barcodes for composite parts


- `dialout_designs.csv`

- `dialout_barcodes.csv`

- `dialout_summary.txt`


## 3. Analyse the dialout process

Analysis of the dialout process can be performed to give statistics on the success of barcoding designs containing different combinations of parts. 

## 4. Primer design for composite part barcodes

The final step in the process is to automatically generate suitable primers for the barcodes found from the sequencing data. We recommend the standard Primer3 settings in the `dialout_primer3_settings.txt` file. This can be used as input for the `dialout_primer_design.py` script with the input barcode file generated from the dialout process, the indexes of the forward and reverse barcode to use and any short sequences that will be present in the backbone of the circularised vector. An example of calling this script is as follows, with the generated primers contained in the `output_primers.csv` file:

	python ./bin/dialout_primer_design.py ./bin/dialout_primer3_settings.txt ./example/input_barcodes.csv 3 2 output_primers.csv GC GC
