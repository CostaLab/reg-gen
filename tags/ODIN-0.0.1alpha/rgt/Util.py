import os
import sys
import ConfigParser

# Reading config file directory
root_dir = os.path.dirname(os.path.abspath(__file__))
data_config_path_file_name = os.path.join(root_dir, 'data.config.path')
data_config_path_file = open(data_config_path_file_name,"r")
data_config_file_name = data_config_path_file.readline().strip()
data_config_path_file.close()

# Parsing config file
config = ConfigParser.ConfigParser()
config.read(data_config_file_name)

# Reading data directory
DATA_DIR = "/".join(data_config_file_name.split("/")[:-1])


class GenomeData:
    ORGANISM = config.get('GenomeData','organism')
    GENOME = os.path.join(DATA_DIR,ORGANISM,config.get('GenomeData','genome'))
    CHROMOSOME_SIZES = os.path.join(DATA_DIR,ORGANISM,config.get('GenomeData','chromosome_sizes'))
    ASSOCIATION_FILE = os.path.join(DATA_DIR,ORGANISM,config.get('GenomeData','association_file'))

class MotifData:
    PWM_DATASET = os.path.join(DATA_DIR,config.get('MotifData','pwm_dataset'))
    LOGO_DATASET = os.path.join(DATA_DIR,config.get('MotifData','logo_dataset'))

class OverlapType:
    OVERLAP = 0 
    ORIGINAL = 1
    COMP_INCL = 2
