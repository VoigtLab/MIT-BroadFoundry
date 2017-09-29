
# Copyright (C) 2017 by
# Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# All rights reserved
# Released under MIT license (see LICENSE.txt)

# Map the raw reads running each sample separately

BIN_PATH=../bin

python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples tube_1
python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples tube_2
python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples tube_3
python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples tube_4
python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples tube_5
python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples tube_6
python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples tube_7
python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples tube_8

python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples flask_1
python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples flask_2
python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples flask_3
python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples flask_4
python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples flask_5
python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples flask_6
python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples flask_7
python $BIN_PATH/map_reads.py -settings ./data/settings.txt -samples flask_8
