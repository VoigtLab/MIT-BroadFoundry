
# Copyright (C) 2017 by
# Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# All rights reserved
# Released under MIT license (see LICENSE.txt)

# Analyse and collate the reads for each sample into a single table and calculate TMM normalisation factors

BIN_PATH=../bin

python $BIN_PATH/read_analysis.py -settings ./data/settings.txt -bin_path $BIN_PATH/
