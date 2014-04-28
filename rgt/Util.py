import os
import sys
import ConfigParser
from pkg_resources import Requirement, resource_filename

"""
The Util classes contains many utilities needed by other classes
such as the paths to input files.

Authors: Eduardo G. Gusmao, Manuel Allhoff, Joseph Kuo and Ivan G. Costa.

"""

class ConfigurationFile:
    """
    Represent the data path configuration file (data.config).

    Authors: Eduardo G. Gusmao.

    It serves as a superclass to classes that will contain
    default variables (such as paths, parameters to tools, etc.)
    for a certain purpose (genomic data, motif data, etc.).

    Variables:

    * self.config: Represents the configuration file.
    * self.data_dir: Represents the root path to data files.

    """

    def __init__(self):
        """
        Initializes ConfigurationFile:
        1. Reads the auxiliary data.config.path file that contains the path to the rgt data folder.
        2. Creates the ConfigParser self.config.
        3. Creates the str self.data_dir.
        """
        
        # Reading config file directory
        root_dir = os.path.dirname(os.path.abspath(__file__))
        data_config_path_file_name = os.path.join(root_dir, 'data.config.path')
        data_config_path_file = open(data_config_path_file_name,"r")
        data_config_file_name = data_config_path_file.readline().strip()
        data_config_path_file.close()

        # Parsing config file
        self.config = ConfigParser.ConfigParser()
        self.config.read(data_config_file_name)

        # Reading data directory
        self.data_dir = os.path.split(data_config_file_name)[0]

class GenomeData(ConfigurationFile):
    """
    Represent genomic data.

    Authors: Eduardo G. Gusmao.

    Inherits ConfigurationFile.

    Methods:

    get_organism():
    Returns the current organism.

    get_genome():
    Returns the current path to the genome fasta file.

    get_chromosome_sizes():
    Returns the current path to the chromosome sizes text file.

    get_association_file():
    Returns the current path to the gene association text file.

    """

    def __init__(self,organism):
        """
        Initializes GenomeData.
        """
        ConfigurationFile.__init__(self)
        self.organism = organism
        self.genome = os.path.join(self.data_dir,self.organism,self.config.get('GenomeData','genome'))
        self.chromosome_sizes = os.path.join(self.data_dir,self.organism,self.config.get('GenomeData','chromosome_sizes'))
        self.association_file = os.path.join(self.data_dir,self.organism,self.config.get('GenomeData','association_file'))

    def get_organism(self):
        """
        Returns the current organism.
        """
        return self.organism

    def get_genome(self):
        """
        Returns the current path to the genome fasta file.
        """
        return self.genome

    def get_chromosome_sizes(self):
        """
        Returns the current path to the chromosome sizes text file.
        """
        return self.chromosome_sizes

    def get_association_file(self):
        """
        Returns the current path to the gene association text file.
        """
        return self.association_file


class MotifData(ConfigurationFile):
    """
    Represent motif (PWM) data.

    Authors: Eduardo G. Gusmao.

    Inherits ConfigurationFile.

    Methods:

    get_repositories_list():
    Returns the current repository list.

    get_pwm_list():
    Returns the list of current paths to the PWM repositories.

    get_logo_list():
    Returns the list of current paths to the logo images of PWMs
    in the given repositories.

    """

    def __init__(self):
        """
        Initializes MotifData.
        """
        ConfigurationFile.__init__(self)
        self.repositories_list = self.config.get('MotifData','repositories').split(",")
        self.pwm_list = []
        self.logo_list = []
        for current_repository in self.repositories_list:
            self.pwm_list.append(os.path.join(self.data_dir,self.config.get('MotifData','pwm_dataset'),current_repository))
            self.logo_list.append(os.path.join(self.data_dir,self.config.get('MotifData','logo_dataset'),current_repository))

    def get_repositories_list(self):
        """
        Returns the current repository list.
        """
        return self.repositories_list

    def get_pwm_list(self):
        """
        Returns the list of current paths to the PWM repositories.
        """
        return self.pwm_list

    def get_logo_list(self):
        """
        Returns the list of current paths to the logo images of PWMs
        in the given repositories.
        """
        return self.logo_list

class OverlapType:
    """
    Class of overlap type constants.

    Authors: Joseph Kuo.

    Constants:

    OVERLAP: Represents #TODO
    ORIGINAL: Represents #TODO
    COMP_INCL: Represents # TODO

    """

    OVERLAP = 0 
    ORIGINAL = 1
    COMP_INCL = 2


