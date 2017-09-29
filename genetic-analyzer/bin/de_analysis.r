#!/usr/bin/env Rscript

# Copyright (C) 2017 by
# Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# All rights reserved
# Released under MIT license (see LICENSE.txt)

# Arguments:
# 1. count_matrix_filename
# 2. group1, e.g., 1,2,5,6
# 3. group2, e.g., 3,4,7,8
# 4. library_size_matrix
# 5. output_file_prefix
args <- commandArgs(TRUE)
arg_count_matrix <- args[1]
arg_group1 <- as.integer(unlist(strsplit(args[2], ",")))
arg_group2 <- as.integer(unlist(strsplit(args[3], ",")))
arg_library_size_matrix <- args[4]
arg_output_file_prefix <- args[5]

# We use edgeR for DE analysis
library(edgeR)

# Load the data file
count_matrix <- read.table(arg_count_matrix, header=T, row.names=1, com='')

# Extract those columns we are interested in
col_ordering <- c(arg_group1, arg_group2)
count_matrix <- count_matrix[,col_ordering]
conditions <- factor(c(rep("group_1", length(arg_group1)), rep("group_2", length(arg_group2))))

# Load the actual mapped reads
library_size_matrix <- read.table(arg_library_size_matrix, header=T, row.names=1, com='')

# Get the library sizes (total counts in annotated genes)
lib_sizes <- as.vector(library_size_matrix$total_mapped_reads)[col_ordering]

# Calculate normalisation factors using TMM
expr <- DGEList(counts=count_matrix, group=conditions, lib.size=lib_sizes)
expr <- calcNormFactors(expr)
expr <- estimateCommonDisp(expr)
expr <- estimateTagwiseDisp(expr)
de_data <- exactTest(expr)
de_results <- topTags(de_data, n=length(count_matrix[,1]))

# write the output to a text file
write.table(as.matrix(de_results$table), col.names=NA, quote=FALSE, file=paste0(arg_output_file_prefix, ".de.analysis.txt"), sep="\t")
