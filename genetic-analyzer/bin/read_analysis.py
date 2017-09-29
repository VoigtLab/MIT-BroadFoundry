#!/usr/bin/env python

# Copyright (C) 2017 by
# Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# All rights reserved
# Released under MIT license (see LICENSE.txt)

# Supporting modules
import argparse
import genetic_analyzer as ga

def main():
	# Parse the command line inputs
	parser = argparse.ArgumentParser(description="read_analysis")
	parser.add_argument("-settings",  dest="settings",  required=True,  help="settings.txt", metavar="string")
	parser.add_argument("-bin_path",  dest="bin_path",  required=True,  help="/usr/bin/", metavar="string")
	args = parser.parse_args()
	# Run the command
	settings = ga.load_settings(args.settings)
	cur_bin_path = args.bin_path
	# Collate read counts and mapped reads for all samples
	counts = {}
	mapped_reads = {}
	sample_names = []
	gene_lengths = {}
	for s in sorted(settings.keys()):
		if s != 'None':
			sample_names.append(s)
			counts[s] = ga.read_count_file(ga.count_filename(settings, s))
			mapped_reads[s] = ga.load_mapped_reads(settings, s)
			gene_lengths[s] = ga.load_gene_lengths(settings, s)
	count_matrix = ga.combine_counts(counts, sample_names)
	ga.save_count_matrix(count_matrix, sample_names, ga.count_matrix_filename(settings))
	ga.save_mapped_reads_matrix(mapped_reads, sample_names, ga.mapped_reads_matrix_filename(settings))
	ga.save_gene_length_matrix(gene_lengths, ga.gene_length_matrix_filename(settings))
	# Generate normalisation factors (edgeR) and output RPKM/FPKM values
	ga.norm_fpkm(settings, bin_path=cur_bin_path)

if __name__ == "__main__":
	main()
