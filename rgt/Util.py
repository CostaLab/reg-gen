import os
import sys
import ConfigParser

root_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
data_dir = os.path.join(root_dir,'data')
config_file = os.path.join(root_dir, 'data.config')
print root_dir
print data_dir
print config_file

config = ConfigParser.ConfigParser()
config.read(config_file)

class GenomeData:
    ORGANISM = config.get('GenomeData','organism')
    GENOME = os.path.join(data_dir,ORGANISM,config.get('GenomeData','genome'))
    CHROMOSOME_SIZES = os.path.join(data_dir,ORGANISM,config.get('GenomeData','chromosome_sizes'))
    ASSOCIATION_FILE = os.path.join(data_dir,ORGANISM,config.get('GenomeData','association_file'))

class MotifData:
    PWM_DATASET = os.path.join(data_dir,config.get('MotifData','pwm_dataset'))
    LOGO_DATASET = os.path.join(data_dir,config.get('MotifData','logo_dataset'))

class OverlapType:
    OVERLAP = 0 
    ORIGINAL = 1
    COMP_INCL = 2
