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
	parser = argparse.ArgumentParser(description="de_analysis")
	parser.add_argument("-settings",  dest="settings",  required=True,  help="settings.txt", metavar="string")
	parser.add_argument("-group1",  dest="group1",  required=True,  help="1,2,6,7", metavar="string")
	parser.add_argument("-group2",  dest="group2",  required=True,  help="3,4,8,9", metavar="string")
	parser.add_argument("-output_prefix",  dest="output_prefix",  required=True,  help="WT_vs_Sample", metavar="string")
	parser.add_argument("-bin_path",  dest="bin_path",  required=True,  help="/usr/bin/", metavar="string")
	args = parser.parse_args()
	# Run the command
	settings = ga.load_settings(args.settings)
	cur_group1 = args.group1
	cur_group2 = args.group2
	cur_output_prefix = args.output_prefix
	cur_bin_path = args.bin_path
	# Run the DE analysis
	ga.de_analysis(settings, cur_group1, cur_group2, cur_output_prefix, bin_path=cur_bin_path)

if __name__ == "__main__":
	main()
