
# Copyright (C) 2017 by
# Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# All rights reserved
# Released under MIT license (see LICENSE.txt)

# Perform part analysis (calculate individual strengths and performance) for all promoters, terminators and ribozymes specified in the GFF file.

BIN_PATH=../bin

python $BIN_PATH/part_profile_analysis.py -settings ./data/settings.txt
