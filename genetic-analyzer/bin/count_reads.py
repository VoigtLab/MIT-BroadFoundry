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
	parser = argparse.ArgumentParser(description="count_reads")
	parser.add_argument("-settings",  dest="settings",  required=True,  help="settings.txt", metavar="string")
	parser.add_argument("-samples",  dest="samples",  required=True,  help="1,2", metavar="string")
	parser.add_argument("-feature",  dest="feature",  required=True,  help="gene", metavar="string")
	parser.add_argument("-attribute",  dest="attribute",  required=True,  help="Name", metavar="string")
	parser.add_argument("-strand_opt",  dest="strand_opt",  required=True,  help="no/yes/reverse", metavar="string")
	args = parser.parse_args()
	# Run the command
	samples = args.samples.split(',')
	settings = ga.load_settings(args.settings)
	f = args.feature
	a = args.attribute
	s_opt = args.strand_opt
	for s in samples:
		ga.count_reads(settings, s, feature=f, attribute=a, strand_opt=s_opt)
		ga.gene_lengths(settings, s, feature=f, attribute=a)
		ga.mapped_reads(settings, s)

if __name__ == "__main__":
	main()
