#!/usr/bin/env python
"""
Part read analysis from dialout
===============================
    Generates summary of read statistics for each position in the design.
    Allows for potential biases in part representation to be found.
"""
from __future__ import print_function, division
import os
import sys
import string
import timeit

__author__  = 'Thomas E. Gorochowski, Voigt Lab, MIT'
__license__ = 'MIT'
__version__ = '1.0'

## HELPERS ====================================================================

# Reverse complement
def revcomp(seq, trans=string.maketrans("ACGT", "TGCA")):
    return "".join(reversed(seq.translate(trans)))

## MAIN =======================================================================

start_time = timeit.default_timer()

# Parse command line parameters
if len(sys.argv) != 3:
    print("Usage: python {} <dialout designs> <output file>".format(sys.argv[0]), file=sys.stderr)
    sys.exit()
design_filename, output_filename = sys.argv[1:]
output_filename = output_filename.strip()

# Extract part use
part_use = []
pair_scars = {}
pair_pros = {}
double_pro = False
with open(design_filename, "rU") as design_file:
    # Ignore header
    design_file.next()
    for row in design_file:
        row_parts = row.split(",")
        if len(row_parts) >= 3:
            # Break-up 
            parts = row_parts[0].split("_")
            if part_use == []:
                # Initialize if necessary
                for idx in range(len(parts)):
                    part_use.append({})
            else:
                # Update counts at each position
                for idx in range(len(parts)):
                    assert len(parts) == len(part_use)

                    if parts[idx] not in part_use[idx].keys():
                        part_use[idx][parts[idx]] = [int(row_parts[1]), int(row_parts[2])]
                    else:
                        part_use[idx][parts[idx]][0] += int(row_parts[1])
                        part_use[idx][parts[idx]][1] += int(row_parts[2])
            # Process the possible part pairs
            scar_pair = "("+parts[0]+" : "+parts[-1]+")"
            if scar_pair not in pair_scars.keys():
                pair_scars[scar_pair] = [int(row_parts[1]), int(row_parts[2])]
            else:
                pair_scars[scar_pair][0] += int(row_parts[1])
                pair_scars[scar_pair][1] += int(row_parts[2])
            # Check is promoter pairs possible
            if len(parts) == 5:
                double_pro = True
                pro_pair = "("+parts[1]+" : "+parts[2]+")"
                if pro_pair not in pair_pros.keys():
                    pair_pros[pro_pair] = [int(row_parts[1]), int(row_parts[2])]
                else:
                    pair_pros[pro_pair][0] += int(row_parts[1])
                    pair_pros[pro_pair][1] += int(row_parts[2])

# Output part statistics data
# Barcode, # Designs, Design Names...
print("Writing part statistics...", file=sys.stdout)
parts_file = open(output_filename, "w")
parts_file.write("Position,Part,# Barcode Reads,# Unique Barcode Reads\n")
for idx in range(len(part_use)):
    for part in sorted(part_use[idx].keys()):
        out_list = [str(idx+1), part] + [str(x) for x in part_use[idx][part]]
        parts_file.write(",".join(out_list) + "\n")
parts_file.write("\n")
parts_file.write("Scar Pair,# Barcode Reads,# Unique Barcode Reads\n")
for p in pair_scars.keys():
    out_list = [p] + [str(x) for x in pair_scars[p]]
    parts_file.write(",".join(out_list) + "\n")
if double_pro == True:
    parts_file.write("\n")
    parts_file.write("Promoter Pair,# Barcode Reads,# Unique Barcode Reads\n")
    for p in pair_pros.keys():
        out_list = [p] + [str(x) for x in pair_pros[p]]
        parts_file.write(",".join(out_list) + "\n")
parts_file.close()

stop_time = timeit.default_timer()
print("Done ({0:.2f} seconds)".format(stop_time-start_time), file=sys.stdout)
