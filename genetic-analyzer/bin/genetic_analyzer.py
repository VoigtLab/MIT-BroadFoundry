#!/usr/bin/env python

# Copyright (C) 2017 by
# Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# All rights reserved
# Released under MIT license (see LICENSE.txt)

# Requires following software is available and in current path:
#   - BWA
#   - SAMTools
#   - HTSeq
#   - BEDTools
#   - R + edgeR package (RScript)

# Required modules
import csv
import subprocess
import numpy as np
import re
import math
import random
import scipy.optimize

def bwa_index_filename (settings, sample):
	return settings[sample]['temp_path']+sample

def sam_filename (settings,  sample):
	return settings[sample]['temp_path']+sample+'.sam'

def bam_filename (settings,  sample, extension=True):
	bam_filename = settings[sample]['temp_path']+sample
	if extension == True:
		bam_filename += '.bam'
	return bam_filename

def fragment_dist_filename (settings, sample):
	return settings[sample]['output_path']+sample+'.fragment.distribution.txt'

def count_filename (settings, sample):
	return settings[sample]['output_path']+sample+'.counts.txt'

def mapped_reads_filename (settings, sample):
	return settings[sample]['output_path']+sample+'.mapped.reads.txt'

def mapped_reads_matrix_filename (settings):
	return settings['None']['output_path']+'mapped.reads.matrix.txt'

def gene_length_filename (settings, sample):
	return settings[sample]['output_path']+sample+'.gene.lengths.txt'

def profile_filename (settings, sample):
	return settings[sample]['output_path']+sample+'.profiles.txt'

def profile_fwd_filename (settings, sample):
	return settings[sample]['output_path']+sample+'.fwd.profiles.txt'

def profile_rev_filename (settings, sample):
	return settings[sample]['output_path']+sample+'.rev.profiles.txt'

def profile_norm_fwd_filename (settings, sample):
	return settings[sample]['output_path']+sample+'.fwd.norm.profiles.txt'

def profile_norm_rev_filename (settings, sample):
	return settings[sample]['output_path']+sample+'.rev.norm.profiles.txt'

def position_fraglen_fwd_filename (settings, sample):
	return settings[sample]['output_path']+sample+'.fwd.position.fraglen.txt'

def position_fraglen_rev_filename (settings, sample):
	return settings[sample]['output_path']+sample+'.rev.position.fraglen.txt'

def count_matrix_filename (settings):
	return settings['None']['output_path']+'counts.matrix.txt'

def gene_length_matrix_filename (settings):
	return settings['None']['output_path']+'gene.lengths.matrix.txt'

def promoter_profile_perf_filename (settings, sample):
	return settings[sample]['output_path']+sample+'.promoter.profile.perf.txt'

def terminator_profile_perf_filename (settings, sample):
	return settings[sample]['output_path']+sample+'.terminator.profile.perf.txt'

def ribozyme_profile_perf_filename (settings, sample):
	return settings[sample]['output_path']+sample+'.ribozyme.profile.perf.txt'

def combined_promoter_profile_perf_filename (settings):
	return settings['None']['output_path']+'promoter.profile.perf.txt'

def combined_terminator_profile_perf_filename (settings):
	return settings['None']['output_path']+'terminator.profile.perf.txt'

def combined_ribozyme_profile_perf_filename (settings):
	return settings['None']['output_path']+'ribozyme.profile.perf.txt'

def combined_fitted_promoter_perf_filename (settings, output_name):
	return settings['None']['output_path']+'fitted.promoter.perf.'+output_name+'.txt'

def load_settings (filename):
	"""Load the settings file
	"""
	settings = {}
	data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
	# Ignore header
	header = next(data_reader)
	# Process each line
	for row in data_reader:
		if len(row) == len(header):
			sample = row[0]
			sample_data = {}
			for el_idx, el in enumerate(header[1:]):
				sample_data[el] = row[el_idx+1]
			settings[sample] = sample_data
	return settings

def map_reads (settings, sample):
	"""Map reads using BWA-MEM
	"""
	# Make the indexes
	cmd_index = 'bwa index' + \
				' -p ' + bwa_index_filename(settings, sample) + \
				' ' + settings[sample]['fasta_file']
	print("Making index: "+cmd_index)
	subprocess.call(cmd_index, shell=True)
	# Perform the mapping
	sam_file = sam_filename(settings, sample)
	cmd_mapping = ''
	if settings[sample]['R2_fastq_file'] == '':
		cmd_mapping = 'bwa mem' + \
					  ' ' + bwa_index_filename(settings, sample) + \
					  ' ' + settings[sample]['R1_fastq_file'] + \
					  ' > ' + sam_file
	else:
		cmd_mapping = 'bwa mem' + \
					  ' ' + bwa_index_filename(settings, sample) + \
					  ' ' + settings[sample]['R1_fastq_file'] + \
					  ' ' + settings[sample]['R2_fastq_file'] + \
					  ' > ' + sam_file
	print("Mapping Reads (BWM-MEM): "+cmd_mapping)
	subprocess.call(cmd_mapping, shell=True)
	# Convert to BAM for some tools
	cmd_to_bam = 'samtools view -bS' + \
				 ' ' + sam_file + \
				 ' | samtools sort' + \
				 ' -o ' + bam_filename(settings, sample, extension=True) + \
				 ' -T ' + bam_filename(settings, sample, extension=False) + \
				 ' -' + \
				 ' && samtools index ' + bam_filename(settings, sample, extension=True)
	print("Converting SAM to position sorted BAM: "+cmd_to_bam)
	subprocess.call(cmd_to_bam, shell=True)

def count_reads (settings, sample, feature='gene', attribute='Name', strand_opt='reverse'):
	"""Count reads falling in a specific feature type, group on an attribute
	"""
	# Use HTSeq to count the reads in specific features
	if settings[sample]['R2_fastq_file'] == '' and strand_opt == 'reverse':
		strand_opt = 'yes'
	cmd_count = 'htseq-count' + \
				' -f bam' + \
				' -s ' + strand_opt + \
				' -a 10' + \
				' -m union' + \
				' -r pos' + \
				' -t ' + feature + \
				' -i ' + attribute + \
				' ' + bam_filename(settings, sample, extension=True) + \
				' ' + settings[sample]['gff_file'] + \
				' > ' + count_filename(settings, sample)
	print("Counting reads: "+cmd_count)
	subprocess.call(cmd_count, shell=True)

def mapped_reads (settings, sample):
	cmd_total = 'samtools view -c -F 4' + \
				' ' + bam_filename(settings, sample) + \
				' > ' + mapped_reads_filename(settings, sample)
	print("Total mapped reads: "+cmd_total)
	subprocess.call(cmd_total, shell=True)

def load_mapped_reads (settings, sample):
	file_in = open(mapped_reads_filename(settings, sample), 'rU')
	file_data = file_in.readlines()
	if len(file_data) > 0:
		return int(file_data[0])
	else:
		return 0

def load_gene_lengths (settings, sample):
	gene_lengths = {}
	data_reader = csv.reader(open(gene_length_filename(settings, sample), 'rU'), delimiter='\t')
	header = next(data_reader)
	for row in data_reader:
		if len(row) == 2:
			gene_lengths[row[0]] = int(row[1])
	return gene_lengths

def read_count_file (filename):
	""" Read the count file generated by HTSeq
		count_data is a dict (tag -> count)
	"""
	count_data = {}
	data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
	for row in data_reader:
		# Check that data exists and is not reporting by HTSeq
		if len(row) == 2:
			if len(row[0])>2 and row[0][0:2] != '__':
				count_data[row[0]] = int(row[1])
	return count_data

def combine_counts (counts, sample_names):
	""" Combine a set of count dictionaries
		counts is a dictorinary of count_data where key is sample name
	"""
	full_tag_list = []
	num_of_samples = len(sample_names)
	# Generate the complete tag list (some samples won't have some tags)
	for sample in sample_names:
		full_tag_list = full_tag_list + list(counts[sample].keys())
	full_tag_list = list(set(full_tag_list))
	# Generate matrix zero matrix
	count_matrix = {}
	for tag in full_tag_list:
		count_matrix[tag] = [0]*num_of_samples
	# Update where count exists
	for sample_idx in range(num_of_samples):
		sample = sample_names[sample_idx]
		for tag in counts[sample].keys():
			count_matrix[tag][sample_idx] = counts[sample][tag]
	return count_matrix

def save_count_matrix (count_matrix, sample_names, filename):
	""" Save a count_matrix with the sample_names to file
	"""
	f_out = open(filename, 'w')
	# Write the header
	f_out.write( '\t'.join(['gene_name']+sample_names)+'\n' )
	for tag in sorted(count_matrix):
		count_strs = [str(x) for x in count_matrix[tag]]
		f_out.write( '\t'.join([tag]+count_strs)+'\n' )
	f_out.close()

def save_mapped_reads_matrix (mapped_reads, sample_names, filename):
	""" Save a mapped_reads_matrix with the sample_names to file
	"""
	f_out = open(filename, 'w')
	f_out.write( 'sample\ttotal_mapped_reads\n' )
	for s in sample_names:
		f_out.write( s+'\t'+str(mapped_reads[s])+'\n' )
	f_out.close()

def save_gene_length_matrix (gene_lengths, filename):
	""" Save a gene_length_matrix of all samples to file
	"""
	f_out = open(filename, 'w')
	f_out.write( 'gene\tlength\n' )
	seen = []
	for s in gene_lengths.keys():
		for gene in gene_lengths[s].keys():
			if gene not in seen:
				f_out.write( gene+'\t'+str(gene_lengths[s][gene])+'\n' )
				seen.append(gene)
	f_out.close()

def count_matrix (settings):
	""" Create and save count matrix
	"""
	counts = {}
	for sample in settings.keys():
		if sample != 'None':
			counts[sample] = read_count_file(count_filename(settings, sample))
	sample_names = counts.keys()
	count_matrix = combine_counts(counts, sample_names)
	save_count_matrix(count_matrix, sample_names, 
					  settings['None']['output_path']+'read_count.matrix')

def gene_lengths (settings, sample, feature='gene', attribute='Name'):
	""" Calculate the gene lengths from set of GTF references
	"""
	len_file = gene_length_filename(settings, sample)
	f_out = open(len_file, 'w')
	f_out.write('gene_name\tlength\n')
	seen = []
	data_reader = csv.reader(open(settings[sample]['gff_file'], 'rU'), delimiter='\t')
	for row in data_reader:
		if len(row) == 9 and row[2] == feature:
			attribs = row[8].split(';')
			for el in attribs:
				key = el.split('=')[0]
				value = el.split('=')[1]
				if key == attribute and value not in seen:
					gene_length = int(row[4])-int(row[3])+1
					f_out.write(value+'\t'+str(gene_length)+'\n')
					seen.append(value)
	f_out.close()

def make_profile (settings, sample, stranded=True):
	""" Calculate transcription profile for given regions in BED file
		http://seqanswers.com/forums/showthread.php?t=29399
	"""
	if settings[sample]['R2_fastq_file'] == '':
		# https://www.biostars.org/p/14378/
		fwd_filename = bam_filename(settings, sample, extension=False) + '.fwd.bam'
		cmd_fwd_coverage = 'samtools view -b -F 20 '+bam_filename(settings, sample, extension=True)+' > '+fwd_filename+' && '+\
			'samtools index '+fwd_filename+' && '+\
			'bedtools coverage -d -b '+fwd_filename+' -a '+settings[sample]['bed_file'] + \
			' > '+profile_fwd_filename(settings, sample)
		print("Making forward profile: "+cmd_fwd_coverage)
		subprocess.call(cmd_fwd_coverage, shell=True)
		rev_filename = bam_filename(settings, sample, extension=False) + '.rev.bam'
		cmd_rev_coverage = 'samtools view -b -f 16 '+bam_filename(settings, sample, extension=True)+' > '+rev_filename+' && '+\
			'samtools index '+rev_filename+' && '+\
			'bedtools coverage -d -b '+rev_filename+' -a '+settings[sample]['bed_file'] + \
			' > '+profile_rev_filename(settings, sample)
		print("Making reverse profile: "+cmd_rev_coverage)
		subprocess.call(cmd_rev_coverage, shell=True)
	else:
		if stranded == True:
			fwd_filename = bam_filename(settings, sample, extension=False) + '.fwd.bam'
			fwd1_filename = bam_filename(settings, sample, extension=False) + '.fwd.1.bam'
			fwd2_filename = bam_filename(settings, sample, extension=False) + '.fwd.2.bam'
			cmd_fwd_coverage = 'samtools view -b -f 83 '+bam_filename(settings, sample, extension=True)+' > '+fwd1_filename+' && '+\
				'samtools index '+fwd1_filename+' && '+\
				'samtools view -b -f 163 '+bam_filename(settings, sample, extension=True)+' > '+fwd2_filename+' && '+\
				'samtools index '+fwd2_filename+' && '+\
				'samtools merge -f '+fwd_filename+' '+fwd1_filename+' '+fwd2_filename+' && '+\
				'samtools index '+fwd_filename+' && '+\
				'bedtools coverage -d -b '+fwd_filename+' -a '+settings[sample]['bed_file'] + \
				' > '+profile_fwd_filename(settings, sample)
			print("Making forward profile: "+cmd_fwd_coverage)
			subprocess.call(cmd_fwd_coverage, shell=True)
			rev_filename = bam_filename(settings, sample, extension=False) + '.rev.bam'
			rev1_filename = bam_filename(settings, sample, extension=False) + '.rev.1.bam'
			rev2_filename = bam_filename(settings, sample, extension=False) + '.rev.2.bam'
			cmd_rev_coverage = 'samtools view -b -f 99 '+bam_filename(settings, sample, extension=True)+' > '+rev1_filename+' && '+\
				'samtools index '+rev1_filename+' && '+\
				'samtools view -b -f 147 '+bam_filename(settings, sample, extension=True)+' > '+rev2_filename+' && '+\
				'samtools index '+rev2_filename+' && '+\
				'samtools merge -f '+rev_filename+' '+rev1_filename+' '+rev2_filename+' && '+\
				'samtools index '+rev_filename+' && '+\
				'bedtools coverage -d -b '+rev_filename+' -a '+settings[sample]['bed_file'] + \
				' > '+profile_rev_filename(settings, sample)
			print("Making reverse profile: "+cmd_rev_coverage)
			subprocess.call(cmd_rev_coverage, shell=True)
		else:
			# Non-stranded library (so combine paired-end reads)
			fwd_filename = bam_filename(settings, sample, extension=False) + '.fwd.bam'
			fwd1_filename = bam_filename(settings, sample, extension=False) + '.fwd.1.bam'
			fwd2_filename = bam_filename(settings, sample, extension=False) + '.fwd.2.bam'
			rev_filename = bam_filename(settings, sample, extension=False) + '.rev.bam'
			rev1_filename = bam_filename(settings, sample, extension=False) + '.rev.1.bam'
			rev2_filename = bam_filename(settings, sample, extension=False) + '.rev.2.bam'
			cmd_fwd_coverage = 'samtools view -b -f 83 '+bam_filename(settings, sample, extension=True)+' > '+fwd1_filename+' && '+\
				'samtools index '+fwd1_filename+' && '+\
				'samtools view -b -f 163 '+bam_filename(settings, sample, extension=True)+' > '+fwd2_filename+' && '+\
				'samtools index '+fwd2_filename+' && '+\
				'samtools view -b -f 99 '+bam_filename(settings, sample, extension=True)+' > '+rev1_filename+' && '+\
				'samtools index '+rev1_filename+' && '+\
				'samtools view -b -f 147 '+bam_filename(settings, sample, extension=True)+' > '+rev2_filename+' && '+\
				'samtools index '+rev2_filename+' && '+\
				'samtools merge -f '+fwd_filename+' '+fwd1_filename+' '+fwd2_filename+' '+rev1_filename+' '+rev2_filename+' && '+\
				'samtools index '+fwd_filename+' && '+\
				'bedtools coverage -d -b '+fwd_filename+' -a '+settings[sample]['bed_file'] + \
				' > '+profile_fwd_filename(settings, sample)
			print("Making forward profile: "+cmd_fwd_coverage)
			subprocess.call(cmd_fwd_coverage, shell=True)
			cmd_rev_coverage = 'samtools view -b -f 99 '+bam_filename(settings, sample, extension=True)+' > '+rev1_filename+' && '+\
				'samtools index '+fwd1_filename+' && '+\
				'samtools view -b -f 163 '+bam_filename(settings, sample, extension=True)+' > '+fwd2_filename+' && '+\
				'samtools index '+fwd2_filename+' && '+\
				'samtools view -b -f 99 '+bam_filename(settings, sample, extension=True)+' > '+rev1_filename+' && '+\
				'samtools index '+rev1_filename+' && '+\
				'samtools view -b -f 147 '+bam_filename(settings, sample, extension=True)+' > '+rev2_filename+' && '+\
				'samtools index '+rev2_filename+' && '+\
				'samtools merge -f '+rev_filename+' '+fwd1_filename+' '+fwd2_filename+' '+rev1_filename+' '+rev2_filename+' && '+\
				'samtools index '+rev_filename+' && '+\
				'bedtools coverage -d -b '+rev_filename+' -a '+settings[sample]['bed_file'] + \
				' > '+profile_rev_filename(settings, sample)
			print("Making reverse profile: "+cmd_rev_coverage)
			subprocess.call(cmd_rev_coverage, shell=True)

def norm_fpkm (settings, bin_path=''):
	""" Calculate normalised RPKM/FPKM expression levels
	"""
	cmd_fpkm = 'Rscript '+bin_path+'norm_fpkm.r ' + \
			   ' ' + count_matrix_filename(settings) + \
			   ' ' + mapped_reads_matrix_filename(settings) + \
			   ' ' + gene_length_matrix_filename (settings) + \
			   ' ' + settings['None']['output_path']
	print("Generating normalised RPKM/FPKMs: "+cmd_fpkm)
	subprocess.call(cmd_fpkm, shell=True)

def de_analysis (settings, group1, group2, output_prefix, bin_path=''):
	""" Calculate DEG analysis between two groups (column numbers 1-indexed)
	"""
	cmd_deg = 'Rscript '+bin_path+'de_analysis.r ' + \
			   ' ' + count_matrix_filename (settings) + \
			   ' ' + group1 + \
			   ' ' + group2 + \
			   ' ' + mapped_reads_matrix_filename(settings) + \
			   ' ' + settings['None']['output_path'] + output_prefix
	print("Generating normalised RPKM/FPKMs: "+cmd_deg)
	subprocess.call(cmd_deg, shell=True)

def load_gff (settings, sample):
	""" Load all annotation data from a GFF file
	"""
	gff = {}
	data_reader = csv.reader(open(settings[sample]['gff_file'], 'rU'), delimiter='\t')
	# Process each line
	for row in data_reader:
		if len(row) == 9:
			chromo = row[0]
			part_type = row[2]
			start_bp = int(row[3])
			end_bp = int(row[4])
			part_dir = row[6]
			part_attribs = {}
			split_attribs = row[8].split(';')
			part_name = None
			for attrib in split_attribs:
				key_value = attrib.split('=')
				if len(key_value) == 2:
					if key_value[0] == 'Name':
						part_name = key_value[1]
					else:
						part_attribs[key_value[0]] = key_value[1]
			if part_name != None:
				if chromo not in gff.keys():
					gff[chromo] = {}
				gff[chromo][part_name] = [part_type, part_dir, start_bp, end_bp, part_attribs]
	return gff

def load_profiles (settings, sample, normed=False):
	""" Profiles have the form of a list chr: [start_bp, end_bp, [profile_fwd],[profile_rev]]
	"""
	profiles = {}
	fwd_profile_filename = profile_fwd_filename(settings, sample)
	if normed == True:
		fwd_profile_filename = profile_norm_fwd_filename(settings, sample)
	data_reader = csv.reader(open(fwd_profile_filename, 'rU'), delimiter='\t')
	# Process each line in fwd profile
	for row in data_reader:
		if len(row) == 5:
			cur_chrom = row[0]
			if cur_chrom not in profiles.keys():
				profiles[cur_chrom] = []
			cur_start_bp = int(row[1])
			cur_end_bp = int(row[2])
			cur_profile = find_profile(profiles, cur_chrom, cur_start_bp, cur_end_bp)
			if cur_profile == None:
				new_profile = [cur_start_bp, cur_end_bp, np.zeros(cur_end_bp-cur_start_bp), np.zeros(cur_end_bp-cur_start_bp)]
				new_profile[2][int(row[3])-1] = float(row[4])
				profiles[cur_chrom].append(new_profile)
			else:
				cur_profile[0][int(row[3])-1] = float(row[4])
	rev_profile_filename = profile_rev_filename(settings, sample)
	if normed == True:
		rev_profile_filename = profile_norm_rev_filename(settings, sample)
	data_reader = csv.reader(open(rev_profile_filename, 'rU'), delimiter='\t')
	# Process each line in fwd profile
	for row in data_reader:
		if len(row) == 5:
			cur_chrom = row[0]
			if cur_chrom not in profiles.keys():
				profiles[cur_chrom] = []
			cur_start_bp = int(row[1])
			cur_end_bp = int(row[2])
			cur_profile = find_profile(profiles, cur_chrom, cur_start_bp, cur_end_bp)
			if cur_profile != None:
				cur_profile[1][int(row[3])-1] = float(row[4])
	return profiles

def find_profile (profiles, chrom, start_bp, end_bp):
	""" Find a profile for a given chromosome that spans a given range of bp positions
	"""
	if chrom in profiles.keys():
		for el in profiles[chrom]:
			if el[0] == start_bp and el[1] == end_bp:
				return [el[2], el[3]]
	return None

def extract_profile_region (profiles, chrom, start_bp, end_bp):
	""" Extract a region of a transcription profile for a given chromosome and region
	"""
	region = None
	if chrom in profiles.keys():
		for profile in profiles[chrom]:
			full_chrom = False
			if profile[0] == 0 and profile[1] == len(profile[2]):
				full_chrom = True
			if full_chrom == True:
				fwd_profile = list(profile[2])
				rev_profile = list(profile[3])
				profile_len = len(fwd_profile)
				ext_start_fwd = []
				ext_end_fwd = []
				ext_start_rev = []
				ext_end_rev = []
				# The region will exist
				if start_bp < 0:
					# extend the profile at start
					ext_start_fwd = fwd_profile[start_bp:]
					ext_start_rev = rev_profile[start_bp:]
				if end_bp > profile_len:
					# extend the profile at end
					ext_end_fwd = fwd_profile[:(end_bp-profile_len)]
					ext_end_rev = rev_profile[:(end_bp-profile_len)]
				new_start_bp = start_bp
				new_end_bp = end_bp
				if ext_start_fwd != []:
					new_start_bp = 0
					new_end_bp = end_bp+len(ext_start_fwd)
				new_fwd_profile = ext_start_fwd+fwd_profile+ext_end_fwd
				new_rev_profile = ext_start_rev+rev_profile+ext_end_rev
				region = [new_fwd_profile[new_start_bp:new_end_bp], 
						  new_rev_profile[new_start_bp:new_end_bp]]
				break
			else:
				if start_bp >= profile[0] and end_bp <= profile[1]:
					fwd_profile = list(profile[2])
					rev_profile = list(profile[3])
					profile_len = len(fwd_profile)
					region = [fwd_profile[start_bp-profile[0]:end_bp-profile[0]], 
					          rev_profile[start_bp-profile[0]:end_bp-profile[0]]]
					break
	return region

def reverse_region (region):
	""" Reverse a given region
	"""
	return [region[1][::-1], region[0][::-1]]

def avg_fn (data):
	""" The average function to use
	"""
	return np.mean(data)

def characterize_promoter_units (settings, sample, upstream_bp=10, downstream_skip_bp=0, downstream_bp=10, normed=False):
	""" Characterize all promoter units for a given sample
	"""
	profiles = load_profiles(settings, sample, normed=normed)
	char_data = []
	raw_region = None
	gff = load_gff (settings, sample)
	for chrom in gff.keys():
		for part_name in gff[chrom].keys():
			part_data = gff[chrom][part_name]
			if part_data[0] == 'promoter_unit':
				if part_data[1] == '+': 
					raw_region = extract_profile_region(profiles, chrom, 
									 (part_data[2]-1)-upstream_bp, part_data[3]+downstream_skip_bp+downstream_bp)
				else:
					raw_region = reverse_region(extract_profile_region(profiles, chrom, 
									 (part_data[2]-1)-downstream_skip_bp-downstream_bp, part_data[3]+upstream_bp))
				# Calculate performance
				avg_us = avg_fn(raw_region[0][0:upstream_bp])
				avg_ds = avg_fn(raw_region[0][-downstream_bp:])
				perf = avg_ds-avg_us
				char_data.append([chrom, part_name, avg_us, avg_ds, perf])
	return char_data

def characterize_promoters (settings, sample, upstream_bp=10, downstream_skip_bp=0, downstream_bp=10, normed=False):
	""" Characterize all promoters for a given sample
	"""
	profiles = load_profiles(settings, sample, normed=normed)
	char_data = []
	raw_region = None
	gff = load_gff (settings, sample)
	for chrom in gff.keys():
		for part_name in gff[chrom].keys():
			part_data = gff[chrom][part_name]
			if part_data[0] == 'promoter':
				if part_data[1] == '+': 
					raw_region = extract_profile_region(profiles, chrom, 
									 (part_data[2]-1)-upstream_bp, part_data[3]+downstream_skip_bp+downstream_bp)
				else:
					raw_region = reverse_region(extract_profile_region(profiles, chrom, 
									 (part_data[2]-1)-downstream_skip_bp-downstream_bp, part_data[3]+upstream_bp))
				# Calculate performance
				avg_us = avg_fn(raw_region[0][0:upstream_bp])
				avg_ds = avg_fn(raw_region[0][-downstream_bp:])
				perf = avg_ds-avg_us
				char_data.append([chrom, part_name, avg_us, avg_ds, perf])
	return char_data

def characterize_terminators (settings, sample, upstream_bp=10, upstream_skip_bp=0, downstream_bp=10, normed=False):
	""" Characterize all terminators for a given sample
	"""
	profiles = load_profiles(settings, sample, normed=normed)
	char_data = []
	raw_region = None
	gff = load_gff (settings, sample)
	for chrom in gff.keys():
		for part_name in gff[chrom].keys():
			part_data = gff[chrom][part_name]
			if part_data[0] == 'terminator':
				if part_data[1] == '+': 
					raw_region = extract_profile_region(profiles, chrom, 
									 (part_data[2]-1)-upstream_skip_bp-upstream_bp, part_data[3]+downstream_bp)
				else:
					raw_region = reverse_region(extract_profile_region(profiles, chrom, 
									 (part_data[2]-1)-downstream_bp, part_data[3]+upstream_skip_bp+upstream_bp))
				# Calculate performance
				if raw_region != None:
					avg_us = avg_fn(raw_region[0][0:upstream_bp])
					avg_ds = avg_fn(raw_region[0][-downstream_bp:])
					max_term = 'N'
					t_e = 0.0
					if avg_us == 0.0:
						max_term = 'Y'
						t_e = 0.0
					else:
						if avg_ds < 1.0:
							t_e = 1.0/float(avg_us)
							max_term = 'Y'
						else:
							t_e = float(avg_ds)/float(avg_us)
					if t_e != 0.0:
						t_s = 1.0/t_e
					else:
						t_s = -1.0
					char_data.append([chrom, part_name, avg_us, avg_ds, t_e, t_s, max_term])
	return char_data

def characterize_ribozymes (settings, sample, upstream_promoter_bp=10, upstream_bp=10, downstream_skip_bp=0, downstream_bp=10, normed=False):
	""" Characterize all ribozymes for a given sample
	"""
	profiles = load_profiles(settings, sample, normed=normed)
	char_data = []
	raw_region = None
	promoter_region = None
	gff = load_gff (settings, sample)
	for chrom in gff.keys():
		for part_name in gff[chrom].keys():
			part_data = gff[chrom][part_name]
			if part_data[0] == 'ribozyme':
				cut_site = 0
				promoter_start = int(part_data[4]['upstream_promoter_start'])
				if 'cut_site' in part_data[4].keys():
					cut_site = int(part_data[4]['cut_site'])
				if part_data[1] == '+':
					cur_site_bp = (part_data[2]-1)+cut_site
					raw_region = extract_profile_region(profiles, chrom, 
									 cur_site_bp-upstream_bp, cur_site_bp+downstream_skip_bp+downstream_bp)
					promoter_region = extract_profile_region(profiles, chrom, 
									 promoter_start-upstream_promoter_bp, promoter_start)
				else:
					cur_site_bp = (part_data[3])-cut_site
					raw_region = reverse_region(extract_profile_region(profiles, chrom, 
									 cur_site_bp-downstream_skip_bp-downstream_bp, cur_site_bp+upstream_bp))
					promoter_region = reverse_region(extract_profile_region(profiles, chrom, 
									 promoter_start, promoter_start+upstream_promoter_bp))
				# Calculate performance
				avg_promoter = avg_fn(promoter_region[0])
				avg_us = avg_fn(raw_region[0][0:upstream_bp])
				avg_ds = avg_fn(raw_region[0][-downstream_bp:])
				# Correct for input transcription to promoter
				avg_us = avg_us - avg_promoter
				if avg_us < 0.0:
					avg_us = 0.0
				avg_ds = avg_ds - avg_promoter
				if avg_ds < 0.0:
					avg_ds = 0.0
				max_cut = 'N'
				c_e = 0.0
				if avg_ds <= 0.0:
					c_e = 0.0
					max_cut = 'Y'
				else:
					if avg_us <= 0:
						if avg_ds < avg_us:
							c_e = 0.0
						else:
							c_e = 1.0-(1.0/float(avg_ds))
						max_cut = 'Y'
					else:
						if avg_ds < avg_us:
							c_e = 0.0
						else:
							c_e = 1.0-(float(avg_us)/float(avg_ds))
				char_data.append([chrom, part_name, avg_us, avg_ds, c_e, max_cut, cut_site])
	return char_data

def save_characterization_data (settings, sample, data, part_type=None):
	""" Save all characterisation data (promoters, terminators, ribozymes) for a given sample
	"""
	if part_type == 'promoter':
		filename = promoter_profile_perf_filename(settings, sample)
		f_out = open(filename, 'w')
		f_out.write( 'sample\tchromosome\tpart_name\treads_us\treads_ds\treads_strength\n' )
		for d in data:
			f_out.write( sample+'\t'+'\t'.join([str(x) for x in d])+'\n' )
		f_out.close()
	if part_type == 'terminator':
		filename = terminator_profile_perf_filename(settings, sample)
		f_out = open(filename, 'w')
		f_out.write( 'sample\tchromosome\tpart_name\treads_us\treads_ds\tt_e\tt_s\tmax_term\n' )
		for d in data:
			f_out.write( sample+'\t'+'\t'.join([str(x) for x in d])+'\n' )
		f_out.close()
	if part_type == 'ribozyme':
		filename = ribozyme_profile_perf_filename(settings, sample)
		f_out = open(filename, 'w')
		f_out.write( 'sample\tchromosome\tpart_name\treads_us\treads_ds\tc_e\tmax_cut\tcut_site\n' )
		for d in data:
			f_out.write( sample+'\t'+'\t'.join([str(x) for x in d])+'\n' )
		f_out.close()

def combine_promoter_characterizations (settings, samples):
	""" Combine all promoter characterization data across a set of samples
	"""
	data = {}
	for s in samples:
		filename = promoter_profile_perf_filename(settings, s)
		data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
		header = next(data_reader)
		for row in data_reader:
			if row[1] not in data.keys():
				data[row[1]] = {}
			chrom_data = data[row[1]]
			if row[2] not in chrom_data.keys():
				chrom_data[row[2]] = []
			chrom_part_data = chrom_data[row[2]]
			chrom_part_data.append([row[0]]+row[3:])
	f_out = open(combined_promoter_profile_perf_filename(settings), 'w')
	f_out.write('chromosome\tpart_name\tsample\treads_us\treads_ds\treads_strength\n')
	for chrom in sorted(data.keys()):
		chrom_data = data[chrom]
		for part in sorted(chrom_data.keys()):
			chrom_part_data = chrom_data[part]
			for data_rec in chrom_part_data:
				f_out.write( '\t'.join([chrom, part]+data_rec)+'\n' )
	f_out.close()

def combine_terminator_characterizations (settings, samples):
	""" Combine all terminator characterization data across a set of samples
	"""
	data = {}
	for s in samples:
		filename = terminator_profile_perf_filename(settings, s)
		data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
		header = next(data_reader)
		for row in data_reader:
			if row[1] not in data.keys():
				data[row[1]] = {}
			chrom_data = data[row[1]]
			if row[2] not in chrom_data.keys():
				chrom_data[row[2]] = []
			chrom_part_data = chrom_data[row[2]]
			chrom_part_data.append([row[0]]+row[3:])
	f_out = open(combined_terminator_profile_perf_filename(settings), 'w')
	f_out.write('chromosome\tpart_name\tsample\treads_us\treads_ds\tt_e\tt_s\tmax_term\n')
	for chrom in sorted(data.keys()):
		chrom_data = data[chrom]
		for part in sorted(chrom_data.keys()):
			chrom_part_data = chrom_data[part]
			for data_rec in chrom_part_data:
				f_out.write( '\t'.join([chrom, part]+data_rec)+'\n' )
	f_out.close()

def combine_ribozyme_characterizations (settings, samples):
	""" Combine all ribozyme characterization data across a set of samples
	"""
	data = {}
	for s in samples:
		filename = ribozyme_profile_perf_filename(settings, s)
		data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
		header = next(data_reader)
		for row in data_reader:
			if row[1] not in data.keys():
				data[row[1]] = {}
			chrom_data = data[row[1]]
			if row[2] not in chrom_data.keys():
				chrom_data[row[2]] = []
			chrom_part_data = chrom_data[row[2]]
			chrom_part_data.append([row[0]]+row[3:])
	f_out = open(combined_ribozyme_profile_perf_filename(settings), 'w')
	f_out.write('chromosome\tpart_name\tsample\treads_us\treads_ds\tc_e\tmax_cut\tcut_site\n')
	for chrom in sorted(data.keys()):
		chrom_data = data[chrom]
		for part in sorted(chrom_data.keys()):
			chrom_part_data = chrom_data[part]
			for data_rec in chrom_part_data:
				f_out.write( '\t'.join([chrom, part]+data_rec)+'\n' )
	f_out.close()

def fragment_length_dists (settings, sample, reads_to_sample=1000000):
	""" Generate the fragment length distribution for a sample (adapted from get_insert_size.py (Wei Li))
	"""
	frag_file = fragment_dist_filename(settings, sample)
	sam_file = sam_filename(settings, sample)
	f_in = open(sam_file, 'rU')
	plrdlen={}
	plrdspan={}
	objmrl=re.compile('([0-9]+)M$')
	objmtj=re.compile('NH:i:(\d+)')
	nline=0
	with open(sam_file, 'rU') as ins:
		for lines in ins:
			field=lines.strip().split()
			nline=nline+1
			if nline >= reads_to_sample:
				break
			if len(field)<12:
				continue
			try:
				mrl=objmrl.match(field[5])
				if mrl==None: # ignore non-perfect reads
					continue
				readlen=int(mrl.group(1))
				if readlen in plrdlen.keys():
					plrdlen[readlen]=plrdlen[readlen]+1
				else:
					plrdlen[readlen]=1
				if field[6]!='=':
					continue
				dist=int(field[8])
				if dist<=0: # ignore neg dist
					continue
				mtj=objmtj.search(lines)
				if dist in plrdspan.keys():
					plrdspan[dist]=plrdspan[dist]+1
				else:
					plrdspan[dist]=1
			except ValueError:
				continue
	f_out = open(frag_file, 'w')
	for k in sorted(plrdspan.keys()):
		f_out.write(str(k)+'\t'+str(plrdspan[k])+'\n')

def load_norm_factor (settings, sample):
	""" Load edgeR normalization factors from file
	"""
	norm_facs = {}
	norm_fac_file = settings['None']['output_path']+'norm.factors.matrix.txt'
	data_reader = csv.reader(open(norm_fac_file, 'rU'), delimiter='\t')
	# Ignore the header
	next(data_reader)
	# Process each line
	for row in data_reader:
		if len(row) == 3:
			norm_facs[row[0]] = [float(row[1]), float(row[2])]
	return (norm_facs[sample][0]*norm_facs[sample][1])

def load_fragmentation_dist (settings, sample, max_frag_len=1000):
	""" Load the fragmentation distribution from file
	"""
	frag_dist = np.zeros(max_frag_len+1)
	frag_file = fragment_dist_filename(settings, sample)
	data_reader = csv.reader(open(frag_file, 'rU'), delimiter='\t')
	# Process each line
	for row in data_reader:
		frag_len = int(row[0])
		frag_count = int(row[1])
		if frag_len <= max_frag_len:
			frag_dist[frag_len] = frag_count
	return frag_dist

def identify_internal_and_overlapped_fragments(read_dictionary,transcript_unit,plasmid_length,scale):
	""" Identify fragments that are internal to a transcription unit and those that overlap with the edges
	"""
	internal_profile = np.ones(plasmid_length)*scale
	overlapped_profile = np.ones(plasmid_length)*scale
	internal_read = {}
	for i in range(plasmid_length):
		internal_read[i] = []
	overlapped_read = {}
	for i in range(plasmid_length):
		overlapped_read[i] = []
	for rr in range(plasmid_length):
		temp_locations = [0]+[item[0] for item in transcript_unit]+[item[1] for item in transcript_unit]+[plasmid_length-1]
		left_pos  = [item-rr for item in temp_locations if item-rr <= 0]
		right_pos = [item-rr for item in temp_locations if item-rr >= 0]
		closest_left  = max(left_pos)  + rr
		closest_right = min(right_pos) + rr
		if closest_left in [item[0] for item in transcript_unit] and closest_right in [item[1] for item in transcript_unit]: 
			for frag_len in read_dictionary[rr]:
				if frag_len > 1 and (rr + frag_len -1) <= closest_right:
					internal_profile[rr:rr + frag_len] += np.ones(frag_len)
					internal_read[rr].append(frag_len)
				elif frag_len > 1 and (rr + frag_len -1) > closest_right:
					overlapped_profile[rr:rr + frag_len] += np.ones(frag_len)
					overlapped_read[rr].append(frag_len)
		else:
			for frag_len in read_dictionary[rr]:
				if frag_len > 1:
					overlapped_profile[rr:rr + frag_len] += np.ones(frag_len)
					overlapped_read[rr].append(frag_len)  
	return internal_profile,internal_read,overlapped_profile,overlapped_read

def generate_corrected_profiles (corr_factor, fwd_border_dicts, fwd_profile, plasmid_length, scale):
	""" Correct the edges of the transcription profiles
	"""
	# New profile to hold the normalised data
	Normalized_fwd_profile = np.ones(plasmid_length)*scale
	# Cycle through each transcript (no idea why this is a dict)
	for item in sorted(fwd_border_dicts):
		# For each transcript
		for h in range(len(fwd_border_dicts[item])):			
			transcript_len = fwd_border_dicts[item][h][1]-fwd_border_dicts[item][h][0]+1
			correction = np.zeros(transcript_len)
			# Extract the relevant part of the correction profile
			for kk in range(fwd_border_dicts[item][h][0], fwd_border_dicts[item][h][1]+1):
				dist_5 = kk - fwd_border_dicts[item][h][0]
				dist_3 = fwd_border_dicts[item][h][1] - kk
				pp_5 = corr_factor[dist_5]
				pp_3 = corr_factor[dist_3]
				correction[kk-fwd_border_dicts[item][h][0]] = min(pp_5,pp_3)
			Normalized_fwd_profile[fwd_border_dicts[item][h][0]:fwd_border_dicts[item][h][1]+1] += fwd_profile[fwd_border_dicts[item][h][0]:fwd_border_dicts[item][h][1]+1] / correction
	return Normalized_fwd_profile

def get_chromosome_details(settings, sample):
	""" Find the lengths of all chromosomes in a sample
	"""
	chrom_data = {}
	chroms_file = settings[sample]['fasta_file']
	cur_chrom = ''
	cur_seq = ''
	with open(chroms_file, 'rU') as f:
	    for line in f:
	        cur_line = line.strip()
	        if cur_chrom != '':
	        	if len(cur_line) > 1 and cur_line[0:1] == '>':
	        		chrom_data[cur_chrom] = len(cur_seq)
	        		cur_seq = ''
	        		cur_chrom = cur_line[1:].strip()
	        	else:
	        		cur_seq = cur_seq+cur_line
	        else:
	        	if len(cur_line) > 1 and cur_line[0:1] == '>':
	        		cur_chrom = cur_line[1:].strip()
	        		cur_seq = ''
	chrom_data[cur_chrom] = len(cur_seq)
	return chrom_data

def load_all_reads (settings, sample, chrom_data, min_length=10, stranded=True):
	""" Load each qualifying read into data structures for later analysis
	"""
	# Load necessary details
	sam_fn = sam_filename(settings,  sample)
	chroms = chrom_data.keys()
	# Generate basic profiles for all chromosomes
	reads_fwd = {}
	reads_rev = {}
	for c in chroms:
		reads_fwd[c] = {}
		reads_rev[c] = {}
		for l in range(chrom_data[c]):
			reads_fwd[c][l] = [1]
			reads_rev[c][l] = [1]
	if stranded == True:
		# Go through all mapped reads check mapped to chromosome and quality (split by strand)
		nline = 0
		print(sam_fn)
		field = None
		with open(sam_fn, 'rU') as ins:
			for lines in ins:
				field=lines.strip().split()
				# Skip header
				if nline > 1:
					if field[2] in chroms and int(field[4]) >= 10:
						if int(field[1]) == 83 or int(field[1]) == 163:
							if abs(int(field[8])) >= min_length:
								reads_fwd[field[2]][min(int(field[3]),int(field[7]))].append(abs(int(field[8]))-1)
						elif int(field[1]) == 99 or int(field[1]) == 147:
							if abs(int(field[8])) >= min_length:
								reads_rev[field[2]][min(int(field[3]),int(field[7]))].append(abs(int(field[8]))-1)
				nline += 1
	else:
		# Go through all mapped reads check mapped to chromosome and quality (add to both strands profiles - strandless)
		# https://www.biostars.org/p/56246/
		nline = 0
		print(sam_fn)
		field = None
		with open(sam_fn, 'rU') as ins:
			for lines in ins:
				field=lines.strip().split()
				# Skip header
				if nline > 1:
					if field[2] in chroms and int(field[4]) >= 10:
						if int(field[1]) != 4:
							if abs(int(field[8])) >= min_length:
								reads_fwd[field[2]][min(int(field[3]),int(field[7]))].append(abs(int(field[8]))-1)
								reads_rev[field[2]][min(int(field[3]),int(field[7]))].append(abs(int(field[8]))-1)
				nline += 1
	return reads_fwd, reads_rev

def extract_TUs (gff, chrom):
	""" Extract all transcription units from the GFF data
	"""
	TUs_fwd = []
	TUs_rev = []
	for part_name in gff[chrom].keys():
		part_data = gff[chrom][part_name]
		if part_data[0] == 'transcript':
			if part_data[1] == '+':
				TUs_fwd.append((part_data[2], part_data[3]))
			else:
				TUs_rev.append((part_data[2], part_data[3]))
	return TUs_fwd, TUs_rev

def calc_corr_profile (settings, sample, transcript_len, num_to_sample=10000):
	""" Generate the expected profile used for correction
	"""
	profile = np.ones(transcript_len)
	ones = np.ones(transcript_len)
	frag_dist = load_fragmentation_dist(settings, sample, max_frag_len=transcript_len+100)
	# Make into probability distribution (sum=1)
	frag_dist_norm = frag_dist[0:transcript_len+1]/sum(frag_dist[0:transcript_len+1])
	# Sample a large number
	lengths = np.arange(0, transcript_len+1)
	for i in range(int(num_to_sample)):
		# Randomly choose a fragment length
		frag_len = int(np.random.choice(lengths, p=frag_dist_norm))
		# Randomly place within boundaries
		pos_5 = int(random.uniform(0, transcript_len-frag_len))
		profile[pos_5:(pos_5+frag_len)] += ones[pos_5:(pos_5+frag_len)]
	return profile

def create_profiles (settings, sample, chroms, scale=10e-5, stranded=True):
	""" Create and save to file corrected profiles for a given sample
	"""
	norm_fac = load_norm_factor(settings, sample)
	gff = load_gff (settings, sample)
	chrom_data = get_chromosome_details(settings, sample)
	reads_fwd, reads_rev = load_all_reads(settings, sample, chrom_data, min_length=10, stranded=stranded)
	profiles = {}
	for chrom in chroms:
		# Generate the transcript structure as list of (start, end) tuples.
		TUs_fwd, TUs_rev = extract_TUs(gff, chrom)
		# Check that there are any edge effects to correct for
		if TUs_fwd != [] or TUs_rev != []:
			# identify_internal_and_overlapped_fragments
			(internal_fwd_profile,internal_fwd_read,overlapped_fwd_profile,overlapped_fwd_read) = identify_internal_and_overlapped_fragments(reads_fwd[chrom], TUs_fwd, chrom_data[chrom], scale)
			(internal_rev_profile,internal_rev_read,overlapped_rev_profile,overlapped_rev_read) = identify_internal_and_overlapped_fragments(reads_rev[chrom], TUs_rev, chrom_data[chrom], scale)
			# Generate correction profile
			chrom_len = chrom_data[chrom]
			# Fixed length transcript to calculate edge correction for
			transcript_len = 2000
			corr_profile = calc_corr_profile(settings, sample, transcript_len, num_to_sample=100000)
			ratio_corr_profile = corr_profile/max(corr_profile)
			chrom_len = chrom_data[chrom]
			# Correct 500 bp from edge
			corr_factor = np.ones(chrom_len)
			corr_factor[0:500] = ratio_corr_profile[0:500]
			# Forward profiles
			internal_fwd_borders_dict = {}
			for n in range(len(TUs_fwd)):
				internal_fwd_borders_dict[TUs_fwd[n]] = [TUs_fwd[n]]
			overlapped_fwd_borders_dict = {}
			overlapped_fwd_borders_dict[(0,chrom_len-1)] = [(0,chrom_len-1)]
			Normalized_internal_fwd_profile = generate_corrected_profiles(corr_factor,internal_fwd_borders_dict,internal_fwd_profile,chrom_len,scale)
			Normalized_overlapped_fwd_profile = generate_corrected_profiles(corr_factor,overlapped_fwd_borders_dict,overlapped_fwd_profile,chrom_len,scale)
			Final_Normalized_fwd_profile = Normalized_internal_fwd_profile + Normalized_overlapped_fwd_profile
			# Reverse profiles
			internal_rev_borders_dict = {}
			for n in range(len(TUs_rev)):
				internal_rev_borders_dict[TUs_rev[n]] = [TUs_rev[n]]
			overlapped_rev_borders_dict = {}
			overlapped_rev_borders_dict[(0,chrom_len-1)] = [(0,chrom_len-1)]
			Normalized_internal_rev_profile = generate_corrected_profiles(corr_factor,internal_rev_borders_dict,internal_rev_profile,chrom_len,scale)
			Normalized_overlapped_rev_profile = generate_corrected_profiles(corr_factor,overlapped_rev_borders_dict,overlapped_rev_profile,chrom_len,scale)
			Final_Normalized_rev_profile = Normalized_internal_rev_profile + Normalized_overlapped_rev_profile
			profiles[chrom] = [0, chrom_len, Final_Normalized_fwd_profile, Final_Normalized_rev_profile]
	# Finally save profiles to file
	fwd_file = profile_fwd_filename(settings, sample)
	rev_file = profile_rev_filename(settings, sample)
	f_out_fwd = open(fwd_file, 'w')
	f_out_rev = open(rev_file, 'w')
	for chrom in profiles.keys():
		region = profiles[chrom]
		cur_start_bp = region[0]
		cur_end_bp = region[1]
		cur_fwd_profile = region[2]
		cur_rev_profile = region[3]
		for idx in range(len(cur_fwd_profile)):
			cur_fwd_data = [chrom, str(cur_start_bp), str(cur_end_bp), str(idx+1), str(cur_fwd_profile[idx])]
			cur_rev_data = [chrom, str(cur_start_bp), str(cur_end_bp), str(idx+1), str(cur_rev_profile[idx])]
			f_out_fwd.write('\t'.join(cur_fwd_data)+'\n')
			f_out_rev.write('\t'.join(cur_rev_data)+'\n')
	f_out_fwd.close()
	f_out_rev.close()
	# Save the normalised profiles
	fwd_norm_file = profile_norm_fwd_filename(settings, sample)
	rev_norm_file = profile_norm_rev_filename(settings, sample)
	f_out_norm_fwd = open(fwd_norm_file, 'w')
	f_out_norm_rev = open(rev_norm_file, 'w')
	for chrom in profiles.keys():
		region = profiles[chrom]
		cur_start_bp = region[0]
		cur_end_bp = region[1]
		cur_fwd_profile = region[2]
		cur_rev_profile = region[3]
		for idx in range(len(cur_fwd_profile)):
			cur_fwd_data = [chrom, str(cur_start_bp), str(cur_end_bp), str(idx+1), str((cur_fwd_profile[idx]/norm_fac)*1e6)]
			cur_rev_data = [chrom, str(cur_start_bp), str(cur_end_bp), str(idx+1), str((cur_rev_profile[idx]/norm_fac)*1e6)]
			f_out_norm_fwd.write('\t'.join(cur_fwd_data)+'\n')
			f_out_norm_rev.write('\t'.join(cur_rev_data)+'\n')
	f_out_norm_fwd.close()
	f_out_norm_rev.close()

def load_promoter_characterizations (filename, samples):
	""" Load promoter characterization data from file
	"""
	pro_data_perf = {}
	pro_data_full = {}
	data_reader = csv.reader(open(filename, 'rU'), delimiter='\t')
	# Ignore header
	next(data_reader)
	# Process each line
	for row in data_reader:
		cur_chrom = row[0]
		cur_pu = row[1]
		cur_sample = row[2]
		cur_perf = float(row[5])
		cur_all = float(row[4])
		if cur_sample in samples:
			# Characterised increase from promoter
			if cur_sample not in pro_data_perf.keys():
				pro_data_perf[cur_sample] = {}
			if cur_chrom not in pro_data_perf[cur_sample].keys():
				pro_data_perf[cur_sample][cur_chrom] = {}
			pro_data_perf[cur_sample][cur_chrom][cur_pu] = cur_perf
			# Total downstream reads
			if cur_sample not in pro_data_full.keys():
				pro_data_full[cur_sample] = {}
			if cur_chrom not in pro_data_full[cur_sample].keys():
				pro_data_full[cur_sample][cur_chrom] = {}
			pro_data_full[cur_sample][cur_chrom][cur_pu] = cur_all
	return pro_data_perf, pro_data_full

def hill_func (x, Pmin, Pmin_inc, K, n, repress=False):
	""" Hill function used for fitting
	"""
	if x < 0.0:
		x = 0.0
	if repress == True:
		return Pmin + Pmin_inc*( math.pow(K,n) / (math.pow(K,n)+math.pow(x,n)) )
	else: 
		return Pmin + Pmin_inc*( math.pow(x,n) / (math.pow(K,n)+math.pow(x,n)) )

def extract_fit_params (x, P_names, P_types):
	""" Extract the parameters into more user friendly dict
	"""
	fit_params = {}
	cur_idx = 0
	for p_idx in range(len(P_names)):
		p_name = P_names[p_idx]
		p_type = P_types[p_idx]
		if p_type == 'induced':
			# If an induced promoter then 2 parameters
			fit_params[p_name] = {}
			fit_params[p_name]['Pmin'] = x[cur_idx]
			fit_params[p_name]['Pmin_inc'] = x[cur_idx+1]
			cur_idx += 2
		else:
			# If a repessed promoter then 4 parameters
			fit_params[p_name] = {}
			fit_params[p_name]['Pmin'] = x[cur_idx]
			fit_params[p_name]['Pmin_inc'] = x[cur_idx+1]
			fit_params[p_name]['K'] = x[cur_idx+2]
			fit_params[p_name]['n'] = x[cur_idx+3]
			cur_idx += 4
	return fit_params

def promoter_unit_err_func (x, exp_data, exp_data_full_ds, chrom_to_fit, PU_to_fit, chrom_inputs, PU_inputs, P_names, P_types, fac_reduce_size):		
	""" Error function for the promoter units
	"""
	# Extract the parameters into more user friendly dict
	fit_params = extract_fit_params(x, P_names, P_types)
	# Only calculate error based on samples specified
	err_diffs = []
	for sample in exp_data.keys():
		sample_data = exp_data[sample]
		exp_val = exp_data[sample][chrom_to_fit][PU_to_fit]
		# Calculate fitted output assuming additive fluxes
		fit_outs = []
		cur_param_idx = 0
		for p_idx in range(len(P_names)):
			cur_p = P_names[p_idx]
			# Check the type of promoter and calculate 
			if P_types[p_idx] == 'induced':
				# Inducible promoter
				input_val = extract_key_vals(PU_inputs[p_idx])[sample]
				if input_val == 1.0:
					fit_outs.append(fit_params[cur_p]['Pmin']+fit_params[cur_p]['Pmin_inc'])
				else:
					fit_outs.append(fit_params[cur_p]['Pmin'])
			else:
				# Repressor promoter
				input_val = exp_data_full_ds[sample][chrom_inputs[p_idx]][PU_inputs[p_idx]]
				fit_outs.append( hill_func(input_val, 
					                       fit_params[cur_p]['Pmin'], 
					                       fit_params[cur_p]['Pmin_inc'],
					                       fit_params[cur_p]['K'],
					                       fit_params[cur_p]['n'],
					                       repress=True) )
		fit_val = np.sum(fit_outs)
		# To ensure low values are fitted use a log fitting
		# This part can be adapted to tune the fitting performed
		if exp_val > 0.1:
			exp_val = np.log10(exp_val)
		else:
			exp_val = -2.0
		if fit_val > 0.1:
			fit_val = np.log10(fit_val)
		else:
			fit_val = -2.0
		err_diffs.append(((exp_val-fit_val)*10.0))
	# Return SSE
	return np.sum(np.power(err_diffs, 2.0))

def extract_key_vals (data_str):
	""" Extract key -> value pairs into a dictionary from a string
	"""
	split_str = data_str.split(':')
	key_vals = {}
	for el in split_str:
		key_val_split = el.split('>')
		if len(key_val_split) == 2:
			key_vals[key_val_split[0]] = float(key_val_split[1])
	return key_vals

def extract_list (data_str, to_float=False):
	""" Extract a list of floating point values from a string
	"""
	split_str = data_str.split(',')
	if to_float == True:
		split_str = [float(x) for x in split_str]
	return split_str
	
def fit_promoter_response_functions (settings, samples, output_name, fac_reduce_size=1000.0):
	""" Perform a fitting of all gates to the samples and save results to file
	"""
	# All the data needed for the fitting
	promoter_filename = combined_promoter_profile_perf_filename(settings)
	pro_data, pro_data_full_ds = load_promoter_characterizations(promoter_filename, samples)
	gff = load_gff(settings, samples[0])
	# Somewhere to save the results
	fitted_pro_params = {}
	# Parameters required for each fiting
	chrom_to_fit = ''
	PU_to_fit = ''
	chrom_inputs = []
	PU_inputs = []
	P_names = []
	P_types = []
	# Cycle through each promoter unit and fit the internal promoter functions
	for chrom in gff.keys():
		for part_name in gff[chrom].keys():
			part_data = gff[chrom][part_name]
			if part_data[0] == 'promoter_unit':
				part_attribs = part_data[-1]
				# Populate all the parameters of the promoter unit
				chrom_to_fit = chrom
				PU_to_fit = part_name
				chrom_inputs = extract_list(part_attribs['chrom_inputs'])
				P_names = extract_list(part_attribs['promoter_names'])
				P_types = extract_list(part_attribs['promoter_types'])
				PU_inputs = extract_list(part_attribs['promoter_unit_inputs'])
				# Calculate the number of parameters we have
				num_of_params = 0
				for p_idx in range(len(P_names)):
					if P_types[p_idx] == 'induced':
						num_of_params += 2
					else:
						num_of_params += 4
				# Parameters for fit
				x0 = np.zeros(num_of_params)
				# Some contraints (keep everything positive)
				bnds = []
				cur_param_idx = 0
				for p_idx in range(len(P_names)):
					if P_types[p_idx] == 'induced':
						# Initial conditions (start realistic to improve fitting)
						x0[cur_param_idx] = 1.0
						x0[cur_param_idx+1] = 500.0
						# Set bounds
						bnds.append((0.1, None))
						bnds.append((0.1, None))
						cur_param_idx += 2
					else:
						# Initial conditions (start realistic to improve fitting)
						x0[cur_param_idx] = 1.0
						x0[cur_param_idx+1] = 1000.0
						x0[cur_param_idx+2] = 25.0
						x0[cur_param_idx+3] = 1.0
						# Set bounds
						bnds.append((0.01, None)) # Pmin
						bnds.append((10.0, None)) # Pmax
						bnds.append((10.0, None)) # K
						bnds.append((1.0, 4.0))  # n (limited to 4, but can be removed if parts able to accomodate higher cooperativity)
						cur_param_idx += 4
				# methods = BFGS, nelder-mead, Powell, TNC, SLSQP,  L-BFGS-B
				res = scipy.optimize.minimize(promoter_unit_err_func, x0, args=(pro_data, pro_data_full_ds, chrom_to_fit, PU_to_fit, chrom_inputs, PU_inputs, P_names, P_types, fac_reduce_size),
											  bounds=bnds,
					                          method='SLSQP', jac=False,
					                          options={'disp': True, 'maxiter': 100000})
				# Save the fitted parameters
				cur_param_idx = 0
				for p_idx in range(len(P_names)):
					p_name = P_names[p_idx]
					if P_types[p_idx] == 'induced':
						fitted_pro_params[p_name] = {}
						fitted_pro_params[p_name]['Pmin'] = res.x[cur_param_idx]
						fitted_pro_params[p_name]['Pmin_inc'] = res.x[cur_param_idx+1]
						cur_param_idx += 2
					else:
						fitted_pro_params[p_name] = {}
						fitted_pro_params[p_name]['Pmin'] = res.x[cur_param_idx]
						fitted_pro_params[p_name]['Pmin_inc'] = res.x[cur_param_idx+1]
						fitted_pro_params[p_name]['K'] = res.x[cur_param_idx+2]
						fitted_pro_params[p_name]['n'] = res.x[cur_param_idx+3]
						cur_param_idx += 4
	# Save the results to file	
	out_filename = combined_fitted_promoter_perf_filename(settings, output_name)
	f_out = open(out_filename, 'w')
	for p_name in fitted_pro_params.keys():
		f_out.write(p_name)
		for param in fitted_pro_params[p_name].keys():
			f_out.write('\t'+param+'\t'+str(fitted_pro_params[p_name][param]))
		f_out.write('\n')
	f_out.close()
