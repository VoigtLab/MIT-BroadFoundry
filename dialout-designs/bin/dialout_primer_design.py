#!/usr/bin/env python
"""
    Use Primer3 to automatically design primers for all barcoded designs.
    Prerequisites: primer3-py (can install with "pip install primer3-py")
"""
import primer3 as p3
import sys

__author__  = 'Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT'
__license__ = 'MIT'
__version__ = '1.0'

# Parse command line parameters
if len(sys.argv) != 6 and len(sys.argv) != 8:
    print("Usage: python dialout_primer_design.py <primer3_settings> <barcode file> <5'barcode index> <3'barcode index> <output file> [add 5'seq] [add 3'seq]")
    sys.exit()

# Extract parameters
p3_settings_filename = sys.argv[1].strip()
barcodes_filename = sys.argv[2].strip()
seq_bc_idx_5 = int(sys.argv[3].strip())
seq_bc_idx_3 = int(sys.argv[4].strip())
output_filename = sys.argv[5].strip()
seq_add_5 = ''
seq_add_3 = ''
if len(sys.argv) == 8:
    seq_add_5 = sys.argv[6].strip()
    seq_add_3 = sys.argv[7].strip()

# Used to convert setting values to the correct type
def convert_type(s):
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return s

# Load the settings file
p3_settings = {}
f = open(p3_settings_filename, 'r')
for line in f:
    parts = line.split('=')
    if len(parts) == 2:
        parts = [x.strip() for x in parts]
        if parts[0] == 'PRIMER_PRODUCT_SIZE_RANGE':
            range_parts = parts[1].split('-')
            if len(range_parts) == 2:
                range_parts = [x.strip() for x in range_parts]
                p3_range = [ [int(range_parts[0]), int(range_parts[1])] ]
            else:
                print('ERROR: Missing PRIMER_PRODUCT_SIZE_RANGE setting.')
                exit()
            p3_settings[parts[0]] = p3_range
        else:
            p3_settings[parts[0]] = convert_type(parts[1])

# Load the barcodes file
barcode_data = []
f = open(barcodes_filename, 'r')
# Ignore the header
f.readline()
for line in f:
    parts = line.split(',')
    if len(parts) == 3:
        parts = [x.strip() for x in parts]
        barcode_data.append(parts)

# For each barcoded design generate the primer sequences
for bc_idx in range(len(barcode_data)):
    print('Processing design '+str(bc_idx)+' of '+str(len(barcode_data)))
    bc_data = barcode_data[bc_idx]
    # Extract the correct parts of the barcode
    bc_parts = bc_data[1].split('-')
    bc_parts = [x.strip() for x in bc_parts]
    bc_seq_5 = seq_add_5 + bc_parts[seq_bc_idx_5-1]
    bc_seq_3 = bc_parts[seq_bc_idx_3-1] + seq_add_3
    # Create the Primer3 input parameters and sequence
    p3_seq = {}
    p3_seq['SEQUENCE_ID'] = 'SEQ'
    # Generate the template
    template_seq = bc_seq_5 + '[NNN]' + bc_seq_3
    p3_seq['SEQUENCE_TEMPLATE'] = template_seq
    p3_seq['SEQUENCE_EXCLUDED_REGION'] = [len(bc_seq_5),5] # This is the [NNN] part
    # Design the primers
    p3_designs = p3.designPrimers(p3_seq, p3_settings)
    # Check one was found
    p3_primer_left = ''
    p3_primer_right = ''
    if 'PRIMER_LEFT_0_SEQUENCE' in p3_designs.keys() and 'PRIMER_RIGHT_0_SEQUENCE' in p3_designs.keys():
        p3_primer_left = p3_designs['PRIMER_LEFT_0_SEQUENCE']
        p3_primer_right = p3_designs['PRIMER_RIGHT_0_SEQUENCE']
    # Update the barcode data ready for outputing later
    bc_data = bc_data + [template_seq, p3_primer_left, p3_primer_right]
    barcode_data[bc_idx] = bc_data

# Save the updated barcoded designs file with primer sequences
header = 'Design,Barcode,# Reads,Primer3 Seq,5\' Primer,3\' Primer\n' 
f_out = open(output_filename, 'w')
f_out.write(header)
for bc_idx in range(len(barcode_data)):
    f_out.write(','.join(barcode_data[bc_idx])+'\n')
f_out.close()
