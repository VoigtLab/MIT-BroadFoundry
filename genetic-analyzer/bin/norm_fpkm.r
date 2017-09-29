#!/usr/bin/env Rscript

# Copyright (C) 2017 by
# Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# All rights reserved
# Released under MIT license (see LICENSE.txt)

# Arguments:
# 1. count_matrix_filename
# 3. library_size_matrix
# 4. gene_length_matrix
# 5. output_file_prefix
args <- commandArgs(TRUE)
arg_count_matrix <- args[1]
arg_library_size_matrix <- args[2]
arg_gene_length_matrix <- args[3]
arg_output_file_prefix <- args[4]

# We use edgeR for between sample normalisation factors
library(edgeR)

# Load the data file
count_matrix <- read.table(arg_count_matrix, header=T, row.names=1, com='')

# Load the gene lengths (for RPKM calculation) - reorder
gene_length_matrix <- read.table(arg_gene_length_matrix, header=T, row.names=1, com='')
gene_lengths <- gene_length_matrix[order(match(rownames(gene_length_matrix),rownames(count_matrix))),1]

# Load the actual mapped reads
library_size_matrix <- read.table(arg_library_size_matrix, header=T, row.names=1, com='')

# Get the library sizes (total counts in annotated genes)
print(library_size_matrix)
print(library_size_matrix$total_mapped_reads)
print(colSums(count_matrix))
lib_sizes <- as.vector(library_size_matrix$total_mapped_reads)

# Calculate normalisation factors using TMM
expr <- DGEList(counts=count_matrix, lib.size=lib_sizes)
expr <- calcNormFactors(expr)

# Calculate RPKM value (normalisation factors are included by default) log=FALSE
expr_norm <- rpkm(expr, gene.length=gene_lengths)

# Write the RPKM values to a text file
write.table(expr_norm, col.names=NA, quote=FALSE, file=paste0(arg_output_file_prefix, "/fpkm.normed.matrix.txt"), sep="\t")

# Write normalisation factors used
write.table(data.frame("design"=colnames(count_matrix), "factor"=expr$samples$norm.factors, "lib_size"=expr$samples$lib.size), file=paste0(arg_output_file_prefix, "/norm.factors.matrix.txt"), sep='\t', quote=FALSE, row.names=FALSE)
