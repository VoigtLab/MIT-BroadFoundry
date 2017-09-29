
# Copyright (C) 2017 by
# Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# All rights reserved
# Released under MIT license (see LICENSE.txt)

# Differential gene expression analysis studies are performed by selecting the groups for comparison (these are the sample column indexes from the combined mapped reads matrix: ./results/counts.matrix.txt)

BIN_PATH=../bin

# Comparison of all flask vs. tube samples (output file is placed in ./results with prefix flask_vs_tube)
python $BIN_PATH/de_analysis.py -settings ./data/settings.txt -group1 9,10,11,12,13,14,15,16 -group2 1,2,3,4,5,6,7,8 -output_prefix flask_vs_tube -bin_path $BIN_PATH/
