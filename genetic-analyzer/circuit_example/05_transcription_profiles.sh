
# Copyright (C) 2017 by
# Thomas E. Gorochowski <tom@chofski.co.uk>, Voigt Lab, MIT
# All rights reserved
# Released under MIT license (see LICENSE.txt)

# Generate the transcription profile for each sample separately (-chroms defines the chromosome to create the profile for)

BIN_PATH=../bin

python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples tube_1 -chroms 0x58v50
python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples tube_2 -chroms 0x58v50
python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples tube_3 -chroms 0x58v50
python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples tube_4 -chroms 0x58v50
python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples tube_5 -chroms 0x58v50
python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples tube_6 -chroms 0x58v50
python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples tube_7 -chroms 0x58v50
python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples tube_8 -chroms 0x58v50

python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples flask_1 -chroms 0x58v50
python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples flask_2 -chroms 0x58v50
python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples flask_3 -chroms 0x58v50
python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples flask_4 -chroms 0x58v50
python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples flask_5 -chroms 0x58v50
python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples flask_6 -chroms 0x58v50
python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples flask_7 -chroms 0x58v50
python $BIN_PATH/transcription_profile.py -settings ./data/settings.txt -samples flask_8 -chroms 0x58v50
