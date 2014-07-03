from __future__ import print_function
import os
import sys
import ConfigParser
from optparse import OptionParser,BadOptionError,AmbiguousOptionError

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
        self.gencode_annotation = os.path.join(self.data_dir,self.organism,self.config.get('GenomeData','gencode_annotation'))

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

    def get_gencode_annotation(self):
        """
        Returns the current path to the gencode annotation gtf file.
        """
        return self.gencode_annotation


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
        self.mtf_list = []
        for current_repository in self.repositories_list:
            self.pwm_list.append(os.path.join(self.data_dir,self.config.get('MotifData','pwm_dataset'),current_repository))
            self.logo_list.append(os.path.join(self.data_dir,self.config.get('MotifData','logo_dataset'),current_repository))
            self.mtf_list.append(os.path.join(self.data_dir,self.config.get('MotifData','pwm_dataset'),current_repository+".mtf"))

    def get_repositories_list(self):
        """
        Returns the current repository list.
        """
        return self.repositories_list

    def get_pwm_path(self, current_repository):
        """
        Returns the path to a specific motif repository.
        """
        return os.path.join(self.data_dir,self.config.get('MotifData','pwm_dataset'),current_repository)

    def get_pwm_list(self):
        """
        Returns the list of current paths to the PWM repositories.
        """
        return self.pwm_list

    def get_logo_file(self, current_repository):
        """
        Returns the path to a specific logo repository.
        """
        return os.path.join(self.data_dir,self.config.get('MotifData','logo_dataset'),current_repository)

    def get_logo_list(self):
        """
        Returns the list of current paths to the logo images of PWMs
        in the given repositories.
        """
        return self.logo_list

    def get_mtf_path(self, current_repository):
        """
        Returns the path to a specific mtf file.
        """
        return os.path.join(self.data_dir,self.config.get('MotifData','pwm_dataset'),current_repository+".mtf")

    def get_mtf_list(self):
        """
        Returns the list of current paths to the mtf files.
        """
        return self.mtf_list

class HmmData(ConfigurationFile):
    """
    Represent HMM data.

    Authors: Eduardo G. Gusmao.

    Inherits ConfigurationFile.

    Methods:

    get_default_hmm():
    Returns the current repository list.

    """

    def __init__(self):
        """
        Initializes HmmData.
        """
        ConfigurationFile.__init__(self)
        self.default_hmm = os.path.join(self.data_dir,self.config.get('HmmData','default_hmm'))

    def get_default_hmm(self):
        """
        Returns the current default hmm.
        """
        return self.default_hmm

class OverlapType:
    """
    Class of overlap type constants.

    Authors: Joseph Kuo.

    Constants:

    OVERLAP:
        Return new GenomicRegionSet including only the overlapping regions.
    ORIGINAL: 
        Return the regions of original GenomicRegionSet which have any intersections.
    COMP_INCL: 
        Return region(s) of the GenomicRegionSet which are 'completely' included.

    """

    OVERLAP = 0 
    ORIGINAL = 1
    COMP_INCL = 2

class HelpfulOptionParser(OptionParser):
    """An OptionParser that prints full help on errors."""
    def error(self, msg):
        self.print_help(sys.stderr)
        self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))

class PassThroughOptionParser(HelpfulOptionParser):
    """
    An unknown option pass-through implementation of OptionParser.
    When unknown arguments are encountered, bundle with largs and try again,
    until rargs is depleted.
    sys.exit(status) will still be called if a known argument is passed
    incorrectly (e.g. missing arguments or bad argument types, etc.)
    """
    def _process_args(self, largs, rargs, values):
        while rargs:
            try:
                HelpfulOptionParser._process_args(self,largs,rargs,values)
            except (BadOptionError,AmbiguousOptionError), e:
                largs.append(e.opt_str)

class ErrorHandler():
    """
    Handles errors in a standardized way.
    """

    def __init__(self):
        """
        Initializes ErrorHandler.
        """

        self.program_name = os.path.basename(sys.argv[0])

        """
        Error Dictionary Standard:
        Each entry consists of a key+list in the form X:[Y,Z,W] where:
        * X: The key representing the internal error name.
        * Y: Error number.
        * Z: Exit status.
        * W: Error message to be print.
        """
        self.error_dictionary = {
            "DEFAULT_ERROR": [0,0,"Undefined error. Program terminated with exit status 0."],

            "MOTIF_ANALYSIS_OPTION_ERROR": [1,0,"You must define one specific analysis. Run '"+self.program_name+" -h' for help."],

            "FP_WRONG_ARGUMENT": [2,0,"You must provide at least one and no more than one experimental matrix as input argument."],
            "FP_WRONG_EXPMAT": [3,0,"The experimental matrix could not be loaded. Check if it is correctly formatted."],
            "FP_ONE_REGION": [4,0,"You must provide one 'regions' bed file in the experiment matrix."],
            "FP_NO_DNASE": [5,0,"You must provide one 'reads' file termed DNASE in the experiment matrix"],
            "FP_NO_HISTONE": [6,0,"You must provide at least one 'reads' file not termed DNASE (histone modification) in the experiment matrix."],
            "FP_NB_HMMS": [7,0,"You must provide one HMM file or X HMM files where X is the number of histone tracks detected in the experiment matrix."],
            "FP_HMM_FILES": [8,0,"Your HMM file could not be read. If you did not provide any HMM file, you may not have installed scikit correctly."],
            "FP_BB_CREATION": [9,0,"Big Bed file (.bb) could not be created. Check if you have the bedToBigBed script in $PATH."],
            "XXXXXXX": [10,0,"Xxxxxx"]
        }
        self.error_number = 0
        self.exit_status = 1
        self.error_message = 2

        """
        Warning Dictionary Standard:
        Each entry consists of a key+list in the form X:[Y,Z] where:
        * X: The key representing the internal warning name.
        * Y: Warning number.
        * Z: Warning message to be print.
        """
        self.warning_dictionary = {
            "DEFAULT_WARNING": [0,"Undefined warning."],
            "WARNING1": [1,"Warning 1 Test"],
            "WARNING2": [2,"Warning 2 Test"],

            "FP_ONE_REGION": [3,"There are more than one 'regions' file in the experiment matrix. Only the first will be used."],
            "FP_MANY_DNASE": [4,"There are more than one DNASE 'reads' file. Only the first one will be used."],
            "FP_MANY_HISTONE": [5,"It is recomended that no more than three histone modifications should be used."],
            "FP_DNASE_PROC": [6,"The DNase file could not be processed"],
            "FP_HISTONE_PROC": [7,"The Histone file could not be processed"],
            "FP_SEQ_FORMAT": [8,"The DNase+Histone sequence could not be formatted to be input for scikit"],
            "FP_HMM_APPLIC": [9,"The scikit HMM encountered errors when applied"],
            "XXXXXXX": [10,"Xxxxxx"]
        }
        self.warning_number = 0
        self.warning_message = 1

    def throw_error(self, error_type, add_msg = ""):
        """
        Throws the specified error type. If the error type does not
        exist, throws a default error message and exits.
        """

        # Fetching error type
        try:
            error_number = self.error_dictionary[error_type][self.error_number]
            exit_status = self.error_dictionary[error_type][self.exit_status]
            error_message = self.error_dictionary[error_type][self.error_message]
        except KeyError, IndexError:
            error_number = self.error_dictionary["DEFAULT_ERROR"][self.error_number]
            exit_status = self.error_dictionary["DEFAULT_ERROR"][self.exit_status]
            error_message = self.error_dictionary["DEFAULT_ERROR"][self.error_message]

        # Handling error
        complete_error_message = ("--------------------------------------------------\n"
                                  "Error Number: "+str(error_number)+".\n"
                                  "Program: "+self.program_name+".\n"
                                  "Report: "+error_message+" "+add_msg+"\n"
                                  "Behaviour: The program will quit with exit status "+str(exit_status)+".\n"
                                  "--------------------------------------------------")
        print(complete_error_message, file=sys.stderr)
        sys.exit(exit_status)

    def throw_warning(self, warning_type, add_msg = ""):
        """
        Throws the specified warning type. If the warning type does not
        exist, throws a default warning message and exits.
        """

        # Fetching warning type
        try:
            warning_number = self.warning_dictionary[warning_type][self.warning_number]
            warning_message = self.warning_dictionary[warning_type][self.warning_message]
        except KeyError, IndexError:
            warning_number = self.warning_dictionary["DEFAULT_WARNING"][self.warning_number]
            warning_message = self.warning_dictionary["DEFAULT_WARNING"][self.warning_message]

        # Handling warning
        complete_warning_message = ("--------------------------------------------------\n"
                                    "Warning Number: "+str(warning_number)+".\n"
                                    "Program: "+self.program_name+".\n"
                                    "Report: "+warning_message+" "+add_msg+"\n"
                                    "--------------------------------------------------")
        print(complete_warning_message, file=sys.stderr)

class AuxiliaryFunctions:
    """
    Class of auxiliary functions.

    Authors: Eduardo G. Gusmao.

    Methods:

    #TODO
    """

    @staticmethod
    def string_is_int(s):
        """ Verifies if a string is a numeric integer """
        try:
            int(s)
            return True
        except ValueError:
            return False

    @staticmethod
    def string_is_float(s):
        """ Verifies if a string is a numeric float """
        try:
            float(s)
            return True
        except ValueError:
            return False

    @staticmethod
    def correct_standard_bed_score(score):
        """ Makes score between 0 and 1000 """
        return min(max(score,0),1000)


       
