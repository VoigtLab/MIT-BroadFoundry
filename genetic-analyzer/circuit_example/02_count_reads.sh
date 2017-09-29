
# Copyright (C) 2017 by
# Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# All rights reserved
# Released under MIT license (see LICENSE.txt)

# Count mapped reads for each "gene" feature

BIN_PATH=../bin

python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples tube_1 -feature gene -attribute Name -strand_opt reverse
python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples tube_2 -feature gene -attribute Name -strand_opt reverse
python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples tube_3 -feature gene -attribute Name -strand_opt reverse
python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples tube_4 -feature gene -attribute Name -strand_opt reverse
python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples tube_5 -feature gene -attribute Name -strand_opt reverse
python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples tube_6 -feature gene -attribute Name -strand_opt reverse
python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples tube_7 -feature gene -attribute Name -strand_opt reverse
python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples tube_8 -feature gene -attribute Name -strand_opt reverse

python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples flask_1 -feature gene -attribute Name -strand_opt reverse
python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples flask_2 -feature gene -attribute Name -strand_opt reverse
python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples flask_3 -feature gene -attribute Name -strand_opt reverse
python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples flask_4 -feature gene -attribute Name -strand_opt reverse
python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples flask_5 -feature gene -attribute Name -strand_opt reverse
python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples flask_6 -feature gene -attribute Name -strand_opt reverse
python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples flask_7 -feature gene -attribute Name -strand_opt reverse
python $BIN_PATH/count_reads.py -settings ./data/settings.txt -samples flask_8 -feature gene -attribute Name -strand_opt reverse
