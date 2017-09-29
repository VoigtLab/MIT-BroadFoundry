
# Copyright (C) 2017 by
# Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# All rights reserved
# Released under MIT license (see LICENSE.txt)

# Attempt to fit gate response functions to promoter performance values across all states (how inputs and outputs are connected are defined in the GFF file).

BIN_PATH=../bin

# Perform fitting for tube samples
python $BIN_PATH/promoter_fitting.py -settings ./data/settings.txt -output_name tube -samples tube_1,tube_2,tube_3,tube_4,tube_5,tube_6,tube_7,tube_8

# Perform fitting for flask samples
python $BIN_PATH/promoter_fitting.py -settings ./data/settings.txt -output_name flask -samples flask_1,flask_2,flask_3,flask_4,flask_5,flask_6,flask_7,flask_8
