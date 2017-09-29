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
	parser = argparse.ArgumentParser(description="part_profile_analysis")
	parser.add_argument("-settings",  dest="settings",  required=True,  help="settings.txt", metavar="string")
	args = parser.parse_args()
	settings = ga.load_settings(args.settings)
	samples = []
	for s in settings.keys():
		if s != 'None':
			samples.append(s)
			chr_promoters = ga.characterize_promoter_units(settings, s, upstream_bp=10, downstream_skip_bp=0, downstream_bp=10, normed=True)
			ga.save_characterization_data(settings, s, chr_promoters, part_type='promoter')
			chr_terminators =  ga.characterize_terminators(settings, s, upstream_bp=10, upstream_skip_bp=0, downstream_bp=10, normed=True)
			ga.save_characterization_data(settings, s, chr_terminators, part_type='terminator')
			chr_ribozymes = ga.characterize_ribozymes(settings, s, upstream_promoter_bp=10, upstream_bp=10, downstream_skip_bp=0, downstream_bp=10, normed=True)
			ga.save_characterization_data(settings, s, chr_ribozymes, part_type='ribozyme')
	ga.combine_promoter_characterizations(settings, samples)
	ga.combine_terminator_characterizations(settings, samples)
	ga.combine_ribozyme_characterizations(settings, samples)

if __name__ == "__main__":
	main()
