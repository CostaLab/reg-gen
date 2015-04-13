from __future__ import print_function
import os
import sys
import ConfigParser
from optparse import OptionParser,BadOptionError,AmbiguousOptionError
import shutil

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
        self.annotation_dump_dir = os.path.join(self.data_dir,self.organism)
        self.gene_alias = os.path.join(self.data_dir,self.organism,self.config.get('GenomeData','gene_alias'))

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

    def get_annotation_dump_dir(self):
        """
        Returns the current path to the gencode annotation gtf file.
        """
        return self.annotation_dump_dir

    def get_gene_alias(self):
        """
        Returns the current path to the gene alias txt file.
        """
        return self.gene_alias


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

    get_mtf_list():
    Returns the list of current paths to the mtf (motif annotation) files
    in the given repositories.

    get_fpr_list():
    Returns the list of current paths to the fpr (motif thresholds) files
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
        self.fpr_list = []
        for current_repository in self.repositories_list:
            self.pwm_list.append(os.path.join(self.data_dir,self.config.get('MotifData','pwm_dataset'),current_repository))
            self.logo_list.append(os.path.join(self.data_dir,self.config.get('MotifData','logo_dataset'),current_repository))
            self.mtf_list.append(os.path.join(self.data_dir,self.config.get('MotifData','pwm_dataset'),current_repository+".mtf"))
            self.fpr_list.append(os.path.join(self.data_dir,self.config.get('MotifData','pwm_dataset'),current_repository+".fpr"))

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

    def get_fpr_list(self):
        """
        Returns the list of current paths to the fpr files.
        """
        return self.fpr_list

class HmmData(ConfigurationFile):
    """
    Represent HMM data.

    Authors: Eduardo G. Gusmao.

    Inherits Co7nfigurationFile.

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

class ImageData(ConfigurationFile):
    """
    Represent image data.

    Authors: Eduardo G. Gusmao.

    Inherits ConfigurationFile.

    Methods:

    get_rgt_logo():
    Returns the rgt logo image file location.

    get_css_file():
    Returns the css file location.

    get_default_motif_logo():
    Returns the default motif logo file location.
    """

    def __init__(self):
        """
        Initializes ImageData.
        """
        ConfigurationFile.__init__(self)
        self.rgt_logo = os.path.join(self.data_dir,"fig","rgt_logo.gif")
        self.css_file = os.path.join(self.data_dir,"fig","style.css")
        self.default_motif_logo = os.path.join(self.data_dir,"fig","default_motif_logo.png")
        self.tablesorter = os.path.join(self.data_dir,"fig","jquery.tablesorter.min.js")
        self.jquery = os.path.join(self.data_dir,"fig","jquery-1.11.1.js")
        self.jquery_metadata = os.path.join(self.data_dir,"fig","jquery.metadata.js")
        self.tdf_logo = os.path.join(self.data_dir,"fig","tdf_logo.png")
        self.viz_logo = os.path.join(self.data_dir,"fig","viz_logo.png")


    def get_rgt_logo(self):
        """
        Returns the rgt logo image file location.
        """
        return self.rgt_logo

    def get_css_file(self):
        """
        Returns the css file location.
        """
        return self.css_file

    def get_default_motif_logo(self):
        """
        Returns the default motif logo file location.
        """
        return self.default_motif_logo
    
    def get_sorttable_file(self):
        """
        Returns the default sorttable code location.
        """
        return self.default_motif_logo
    
    def get_jquery(self):
        """
        Returns the default sorttable code location.
        """
        return self.jquery

    def get_tablesorter(self):
        """
        Returns the default sorttable code location.
        """
        return self.tablesorter
    
    def get_jquery_metadata(self):
        """
        Returns the default sorttable code location.
        """
        return self.jquery_metadata
        
    def get_tdf_logo(self):
        """
        Returns the default TDF logo.
        """
        return self.tdf_logo

    def get_viz_logo(self):
        """
        Returns the default TDF logo.
        """
        return self.viz_logo
    
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

class SequenceType:
    """
    Class of sequence type
    Author: Joseph Kuo
    
    Constants:
    
    DNA, RNA
    """
    DNA = 0
    RNA = 1

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
                pass
                #largs.append(e.opt_str)

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
            "FP_WRONG_EXPMAT": [3,0,"The experimental matrix could not be loaded. Check if it is correctly formatted and that your python version is >= 2.7."],
            "FP_ONE_REGION": [4,0,"You must provide one 'regions' bed file in the experiment matrix."],
            "FP_NO_DNASE": [5,0,"You must provide one 'reads' file termed DNASE in the experiment matrix."],
            "FP_NO_HISTONE": [6,0,"You must provide at least one 'reads' file not termed DNASE (histone modification) in the experiment matrix."],
            "FP_NB_HMMS": [7,0,"You must provide one HMM file or X HMM files where X is the number of histone tracks detected in the experiment matrix."],
            "FP_HMM_FILES": [8,0,"Your HMM file could not be read. If you did not provide any HMM file, you may not have installed scikit correctly."],
            "FP_BB_CREATION": [9,0,"Big Bed file (.bb) could not be created. Check if you have the bedToBigBed script in $PATH."],
            "MM_OUT_FOLDER_CREATION": [10,0,"Could not create output folder."],
            "MM_NO_ARGUMENT": [11,0,"Could not read the arguments. Make sure you provided an experimental matrix."],
            "MM_WRONG_EXPMAT": [12,0,"The experimental matrix could not be loaded. Check if it is correctly formatted and that your python version is >= 2.7."],
            "MM_WRONG_RANDPROP": [13,0,"Proportion of random regions is too low (<= 0)."],
            "MM_LOW_NPROC": [14,0,"Number of processor is too low (<= 0)."],
            "ME_OUT_FOLDER_CREATION": [15,0,"Could not create output folder."],
            "ME_FEW_ARG": [16,0,"There are too few arguments. Please use -h option in order to verify command usage."],
            "ME_WRONG_EXPMAT": [17,0,"The experimental matrix could not be loaded. Check if it is correctly formatted and that your python version is >= 2.7."],
            "ME_MATCH_NOTFOUND": [18,0,"Motif matching file was not found. Are you sure you performed --matching before --enrichment?"],
            "ME_BAD_MATCH": [19,0,"Motif matching file is incorrect. Please perform --matching again."],
            "ME_LOW_NPROC": [20,0,"Number of processor is too low (<= 0)."],
            "ME_RAND_NOTFOUND": [21,0,"Random regions not found. Are you sure --matching was correctly performed?"],
            "ME_BAD_RAND": [22,0,"Could not read random regions."],
            "ME_RAND_NOT_BED_BB": [23,0,"Random regions are not in bed or bigbed format."],
            "MM_PSEUDOCOUNT_0": [24,0,"There were errors involved in the creation of Position Weight Matrices. Some distributions of numpy  and/or scipy does not allow for pseudocounts == 0. Please increase the pseudocount (or use default value of 0.1) and try again."],
            "XXXXXXX": [25,0,"Xxxxxx"]
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
            "FP_ONE_REGION": [1,"There are more than one 'regions' file in the experiment matrix. Only the first will be used."],
            "FP_MANY_DNASE": [2,"There are more than one DNASE 'reads' file. Only the first one will be used."],
            "FP_MANY_HISTONE": [3,"It is recomended that no more than three histone modifications should be used."],
            "FP_DNASE_PROC": [4,"The DNase file could not be processed."],
            "FP_HISTONE_PROC": [5,"The Histone file could not be processed."],
            "FP_SEQ_FORMAT": [6,"The DNase+Histone sequence could not be formatted to be input for scikit."],
            "FP_HMM_APPLIC": [7,"The scikit HMM encountered errors when applied."],
            "MM_MANY_ARG": [8,"There are more than one arguments, only the first experimental matrix will be used."],
            "ME_MANY_ARG": [9,"There are more than two arguments, only the two first arguments will be used."],
            "ME_MANY_GENESETS": [10,"There are more than one geneset associated to a group, only the first one will be used."],
            "ME_FEW_GENESETS": [11,"There seems to be a geneset column, but no geneset was found for some entries."],
            "XXXXXXX": [12,"Xxxxxx"]
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

class Html:
    """
    Represent an HTML file.

    Authors: Eduardo G. Gusmao.
    """

    def __init__(self, name, links_dict, fig_dir=None, fig_rpath="../fig", cluster_path_fix="", 
                 RGT_header=True, other_logo=None, homepage=None):
        """ 
        Initializes Html.
        IMPORTANT = cluster_path_fix is going to be deprecated soon. Do not use it.

        Variables:
        xxxxx -- Position Frequency Matrix.
        relative_dir -- Define the directory to store CSS file and RGT logo so that the html code can read from it. Default is None.
        
        """

        # Variable initializations
        self.name = name
        self.links_dict = links_dict
        self.cluster_path_fix = cluster_path_fix
        self.document = []
        self.image_data = ImageData()
        self.other_logo = other_logo
        self.homepage = homepage
        
        # Initialize document
        if fig_dir:
            self.copy_relevent_files(fig_dir)
            self.create_header(relative_dir=fig_rpath, RGT_name=RGT_header, other_logo=other_logo)
        else:
            self.create_header(relative_dir=fig_rpath, RGT_name=RGT_header, other_logo=other_logo)
        
        self.add_links()
        
        
    def copy_relevent_files(self, target_dir):
        try:
            os.stat(target_dir)
        except:
            os.mkdir(target_dir)
        shutil.copyfile(src=self.cluster_path_fix+self.image_data.get_rgt_logo(), dst=os.path.join(target_dir,"rgt_logo.gif"))
        shutil.copyfile(src=self.cluster_path_fix+self.image_data.get_css_file(), dst=os.path.join(target_dir,"style.css"))
        #shutil.copyfile(src=self.cluster_path_fix+self.image_data.get_jquery(), dst=os.path.join(target_dir,"jquery-1.11.1.js"))
        shutil.copyfile(src=self.cluster_path_fix+self.image_data.get_tablesorter(), dst=os.path.join(target_dir,"jquery.tablesorter.min.js"))
        #shutil.copyfile(src=self.cluster_path_fix+self.image_data.get_jquery_metadata(), dst=os.path.join(target_dir,"jquery.metadata.js"))
        #shutil.copyfile(src=self.cluster_path_fix+self.image_data.get_tablesorter(), dst=os.path.join(target_dir,"jquery.metadata.js"))
        if self.other_logo:
            if self.other_logo == "TDF":
                shutil.copyfile(src=self.cluster_path_fix+self.image_data.get_tdf_logo(), dst=os.path.join(target_dir,"tdf_logo.png"))
            if self.other_logo == "viz":
                shutil.copyfile(src=self.cluster_path_fix+self.image_data.get_viz_logo(), dst=os.path.join(target_dir,"viz_logo.png"))
            
        
        
        
    def create_header(self, relative_dir=None, RGT_name=True, other_logo=None):
        """ 
        Creates default document header.
        
        Return:
        None -- Appends content to the document.
        """
        self.document.append('<script src="http://ajax.googleapis.com/ajax/libs/jquery/1.11.2/jquery.min.js"></script>') 
            
        if relative_dir:
            #self.document.append('<script type="text/javascript" src="'+relative_dir+'/jquery-1.11.1.js"></script>')
            self.document.append('<script type="text/javascript" src="'+relative_dir+'/jquery.tablesorter.min.js"></script>')
            #self.document.append('<script type="text/javascript" src="'+relative_dir+'/jquery.metadata.js"></script>')
        else:
            #self.document.append('<script type="text/javascript" src="'+self.cluster_path_fix+self.image_data.get_jquery()+'"></script>')
            self.document.append('<script type="text/javascript" src="'+self.cluster_path_fix+self.image_data.get_tablesorter()+'"></script>')
            #self.document.append('<script type="text/javascript" src="'+self.cluster_path_fix+self.image_data.get_jquery_metadata()+'"></script>')
        
        
        self.document.append("<html>")
        self.document.append("<head><meta http-equiv=\"Content-Type\" content=\"text/html\"><title>RGT "+self.name+"</title>")

        self.document.append("<style type=\"text/css\">")
        self.document.append("<!--")
        
        if relative_dir:
            self.document.append("@import url(\""+relative_dir+"/style.css\");")
        else:
            self.document.append("@import url(\""+self.cluster_path_fix+self.image_data.get_css_file()+"\");")
        
        self.document.append("-->")
        self.document.append("</style></head>")
        self.document.append("<body topmargin=\"0\" leftmargin=\"0\" rightmargin=\"0\" bottommargin=\"0\" marginheight=\"0\" marginwidth=\"0\" bgcolor=\"#FFFFFF\">")        
        
        self.document.append("<h3 style=\"background-color:white; border-top:3px solid gray; border-bottom:3px solid gray;\">")
        self.document.append("<table border=\"0\" width=\"100%\" cellspacing=\"0\" cellpadding=\"0\">")
        self.document.append("  <tr>")

        # Logo
        
        if relative_dir:            
            self.document.append("    <td width=\"5%\">")
            if self.homepage: self.document.append("<a href=\""+self.homepage+"\">")
            if other_logo == "TDF":
                self.document.append("    <img border=\"0\" src=\""+relative_dir+"/tdf_logo.png"+"\" width=\"130\" height=\"100\">")
            elif other_logo == "viz":
                self.document.append("    <img border=\"0\" src=\""+relative_dir+"/viz_logo.png"+"\" width=\"130\" height=\"100\">")
            else:
                self.document.append("    <img border=\"0\" src=\""+relative_dir+"/rgt_logo.gif\" width=\"130\" height=\"100\">")
            if self.homepage: self.document.append("</a>")
            self.document.append("    </td>")
            
        else:
            self.document.append("    <td width=\"5%\"><img border=\"0\" src=\""+self.cluster_path_fix+self.image_data.get_rgt_logo()+"\" width=\"130\" height=\"100\"></td>")

        # Gap
        self.document.append("     <td width=\"5%\"></td>")
        # Title
        if RGT_name:
            self.document.append("    <td width=\"90%\"><p align=\"left\"><font color=\"black\" size=\"5\">Regulatory Genomics Toolbox - "+self.name+"</font></td>")
        else:
            self.document.append("    <td width=\"90%\"><p align=\"left\"><font color=\"black\" size=\"5\">"+self.name+"</font></td>")
        
        
        self.document.append("  </tr>")
        self.document.append("</table>")
        self.document.append("</h3>")

    def add_links(self):
        """ 
        Adds all the links.
        
        Return:
        None -- Appends links to the document.
        """
        for k in self.links_dict.keys():

            self.document.append("<a class=\"pure-button\" href=\""+os.path.join(self.cluster_path_fix,self.links_dict[k])+"\">"+k+"</a>")

        #self.document.append("<table border=\"0\" width=\"100%\" cellspacing=\"0\" cellpadding=\"0\">")
        #self.document.append("  <tr>")
        #self.document.append("    <td width=\"100%\"><font color=\"black\" face=\"Arial\" size=\"4\"><b>&nbsp;&nbsp;")
        #link_str = "    "+" &nbsp;&nbsp; |&nbsp;&nbsp; ".join(["<a href=\""+os.path.join(self.cluster_path_fix,self.links_dict[k])+"\">"+k+"</a>" for k in self.links_dict.keys()])
        #self.document.append(link_str)
        #self.document.append("    </b></font></td>")
        #self.document.append("  </tr>")
        #self.document.append("</table>")
            
    
    def create_footer(self):
        """ 
        Adds footer.
        
        Return:
        None -- Appends footer to the document.
        """
        self.document.append("<br><br>")
        self.document.append("<p align=\"center\"><font face=\"Arial\" color=\"#000000\" size=\"2\">")
        self.document.append("For more details please visit the <a href=\"http://www.regulatory-genomics.org/\"> RGT Website </a>")
        self.document.append("</font></p>")
        self.document.append("<h3 style=\"background-color:white; border-top:3px solid gray;\"></h3>")
        self.document.append("</body>")
        self.document.append("</html>")

    def add_heading(self, heading, align = 50, color = "black", face = "Arial", size = 5, bold = True, idtag=None):
        """ 
        Creates a heading.
        
        Keyword arguments:
        heading -- The heading title.
        align -- Alignment of the heading. Can be either an integer (interpreted as left margin) 
                 or string (interpreted as HTML positional argument). (default 50)
        color -- Color of the heading. (default "black")
        face -- Font of the heading. (default "Arial")
        size -- Size of the heading (HTML units [1,7]). (default 5)
        bold -- Whether the heading is bold. (default True)
        id -- Add ID tag in the heading element

        Return:
        None -- Appends heading to the document.
        """

        if idtag:
            idstr = ' id="'+idtag+'"'
        else:
            idstr = ""

        # Creating header
        content_str = ""
        if(isinstance(align,int)): content_str += "<p style=\"margin-left: "+str(align)+"\""+idstr+">"
        elif(isinstance(align,str)): content_str += "<p align=\""+align+"\""+idstr+">"
        else: pass # ERROR
        content_str += "<font color=\""+color+"\" face=\""+face+"\" size=\""+str(size)+"\""+idstr+">"
        if(bold): content_str += "<b>"
        self.document.append(content_str)

        # Printing heading name
        self.document.append(heading)

        # Creating footing
        end_str = ""
        if(bold): end_str += "</b>"
        end_str += "</font></p>"
        self.document.append(end_str)

    def add_zebra_table(self, header_list, col_size_list, type_list, data_table, align = 50, 
                        cell_align = 'center', auto_width=False, colorcode=None, header_titles=None,
                        border_list=None, sortable=False):
        """ 
        Creates a zebra table.

        Keyword arguments:
        header_list -- A list with the table headers in correct order.
        col_size_list -- A list with the column sizes (integers).
        type_list -- A string in which each character represents the type of each row.
                     s = string (regular word or number)
                     i = image
                     l = link 
        data_table -- A table containing the data to be input according to each data type defined.
                      s = string
                      i = tuple containing: ("file name", width) width = an integer
                      l = tuple containing: ("Name","Link")
        align -- Alignment of the heading. Can be either an integer (interpreted as left margin) 
                 or string (interpreted as HTML positional argument). (default 50)
        cell_align -- Alignment of each cell in the table (default center)
        auto_width -- Adjust the column width by the content automatically regardless of defined col size
        colorcode --
        header_titles -- Given a list corresponding to the header_list, which defines all the explanation in hint windows
        border_list -- 
        Return:
        None -- Appends table to the document.
        """
        #if header_notes: self.document.append("<style> .ami div {display:none} .ami:hover div {display:block} </style>")
        
        if not border_list:
            border_list = [""] * len(data_table[0])
        if auto_width: auto= " table-layout: auto"
        else: auto=""


        # Starting table
        type_list = type_list.lower()
        if(isinstance(align,int)): self.document.append("<p style=\"margin-left: "+str(align)+"\">")
        elif(isinstance(align,str)): self.document.append("<p align=\""+align+"\">") 
        else: pass # TODO ERROR
        
        # Table header
        #self.document.append("<table id=\"myTable\" class=\"tablesorter\">")
        
        if sortable: 
            sortableclass=" class=\"tablesorter\""
            tableid="sortable"
        else:
            sortableclass=""
            tableid="hor-zebra"

        self.document.append("<table id=\""+tableid+"\""+sortableclass+auto+">")

        if colorcode:
            for line in colorcode:
                self.document.append(line)

        #############################
        ##### Header ################
        self.document.append("  <thead>")
        if (isinstance(header_list[0], list)):
        # For headers more than one row
            for r, row_list in enumerate(header_list):
                self.document.append("    <tr>")
                header_str = []

                merge_num = [1] * len(row_list)
                for i, value in enumerate(row_list):
                    if value == None:
                        merge_num[last_true] += 1
                        merge_num[i] -= 1
                    else:
                        last_true = i

                for i in range(0,len(row_list)):
                    if merge_num[i] > 1:
                        if header_titles:
                            header_str.append("<th scope=\"col\" width=\""+str(col_size_list[i])+
                                              "\" align=\""+'center'+"\""+" colspan=\""+str(merge_num[i])+"\" "+
                                              "title=\""+header_titles[r][i]+"\""+border_list[i+merge_num[i]-1]+auto+" >"+
                                              row_list[i]+"</th>")
                        else:
                            header_str.append("<th scope=\"col\" width=\""+str(col_size_list[i])+
                                              "\" align=\""+'center'+"\""+" colspan=\""+str(merge_num[i])+"\""+
                                              border_list[i+merge_num[i]-1]+auto+">"+row_list[i]+"</th>")
                        
                    elif merge_num[i] == 0:
                        continue
                    else:
                        if header_titles:
                            header_str.append("<th scope=\"col\" width=\""+str(col_size_list[i])+
                                              "\" align=\""+cell_align+"\" "+
                                              "title=\""+header_titles[r][i]+"\""+border_list[i]+auto+">"+
                                              row_list[i]+"</th>")
                        else:
                            header_str.append("<th scope=\"col\" width=\""+str(col_size_list[i])+
                                              "\" align=\""+cell_align+"\""+border_list[i]+auto+">"+
                                              row_list[i]+"</th>")

                header_str = "    "+"\n    ".join(header_str)
                self.document.append(header_str)
                self.document.append("    </tr>")

        else:
            self.document.append("    <tr>")
            header_str = []
            for i in range(0,len(header_list)):
                if header_titles:
                    header_str.append("<th scope=\"col\" width=\""+str(col_size_list[i])+"\" align=\""+cell_align+"\" "+
                                      "title=\""+header_titles[i]+"\" >"+header_list[i]+"</th>")
                else:
                    header_str.append("<th scope=\"col\" width=\""+str(col_size_list[i])+"\" align=\""+cell_align+"\">"+
                                      header_list[i]+"</th>")
                
            header_str = "    "+"\n    ".join(header_str)
            self.document.append(header_str)
            self.document.append("    </tr>")
        self.document.append("  </thead>")

        
        #############################
        ##### Table body ############
        self.document.append("  <tbody>")
        for i in range(0,len(data_table)):

            # Row type
            if(i%2==0) and not sortable: self.document.append("    <tr class=\"odd\">")
            else: self.document.append("    <tr>")

            # Body data
            for j in range(0,len(data_table[i])):
                if(type_list[j] == "s"):
                    self.document.append("      <td align=\""+cell_align+"\" "+border_list[j]+">"+data_table[i][j]+"</td>")
                elif(type_list[j] == "i"): 
                    self.document.append("      <td align=\""+cell_align+"\"><img src=\""+self.cluster_path_fix+
                                         data_table[i][j][0]+"\" width="+str(data_table[i][j][1])+" ></td>")
                elif(type_list[j] == "l"):
                    self.document.append("      <td align=\""+cell_align+"\"><a href=\""+data_table[i][j][1]+"\">"+
                                         data_table[i][j][0]+"</a></td>")
                else: pass # TODO ERROR

            # Row ending
            self.document.append("    </tr>")
        
        # Finishing table
        self.document.append("</tbody></table></p>")

    def add_fixed_rank_sortable(self):
        """Add jquery for fixing the first column of the sortable table"""
        scripts = ["<script>",
                   "// add custom numbering widget",
                   "$.tablesorter.addWidget({",
                   "    id: 'numbering',",
                   "    format: function(table) {",
                   "        var c = table.config;",
                   "        $('tr:visible', table.tBodies[0]).each(function(i) {",
                   "            $(this).find('td').eq(0).text(i + 1);",
                   "        });",
                   "    }",
                   "});",
                   "",
                   "$('.tablesorter').tablesorter({",
                   "    // prevent first column from being sortable",
                   "    headers: {",
                   "        0: { sorter: false }",
                   "    },",
                   "    // apply custom widget",
                   "    widgets: ['numbering']",
                   "});",
                   "</script>"]
        for s in scripts:
            self.document.append(s)

    def add_figure(self, figure_path, notes=None, align=50, color="black", face="Arial", size=3, 
                   bold=False, width="800", more_images=None):
        """ 
        Add a figure with notes underneath.
        
        Keyword arguments:
        figure_path -- The path to the figure.
        notes -- A list of strings for further explanation
        align -- Alignment of the heading. Can be either an integer (interpreted as left margin) 
                 or string (interpreted as HTML positional argument). (default 50)
        
        Return:
        None -- Appends the figure to the document.
        """        
        if(isinstance(align,int)): img_str = "<p style=\"margin-left: "+str(align)+"\">"
        elif(isinstance(align,str)): img_str = "<p align=\""+str(align)+"\">"
        else: pass # TODO ERROR

        img_str += '<img src="'+ figure_path +'" width='+width+'>'
        
        if more_images:
            for im in more_images:
                img_str += '<img src="'+ im +'" width='+width+'>'
        
        img_str += '</p>'

        self.document.append(img_str)
        if notes:
            if(isinstance(align,int)): 
                note_str = "<p style=\"margin-left: "+str(align)+"\"><font color=\""+color+"\" face=\""+face+"\" size=\""+str(size)+"\ align=\""+ str(align) + "\">"
            elif(isinstance(align,str)):
                note_str = "<p align=\""+str(align)+"\"><font color=\""+color+"\" face=\""+face+"\" size=\""+str(size)+"\ align=\""+ str(align) + "\">"            
            else: pass # TODO ERROR
            
            if(bold): note_str += "<b>"
            for line in notes:
                note_str += line + "<br>" 
            if(bold): note_str += "</b>"
            note_str += "</font></p>"
            self.document.append(note_str)

    def add_free_content(self, content_list):
        """ 
        Adds free HTML to the document.

        Keyword arguments:
        content_list -- List of strings. Each string is interpreted as a line in the HTML document.
        
        Return:
        None -- Appends content to the document.
        """
        for e in content_list: self.document.append(e)

    def add_list(self, list_of_items, ordered=False):
        """
        Add a list to the document
        """
        codes = ""

        if ordered: codes += "<ol>"
        else: codes += "<ul>"
        
        for item in list_of_items:
            codes += "<li style=\"margin-left: 50\">"+item+"</li>"
        
        if ordered: codes += "</ol>"
        else: codes += "</ul>"

        self.document.append(codes)
        
    def write(self, file_name):
        """ 
        Write HTML document to file name.

        Keyword arguments:
        file_name -- Complete file name to write this HTML document.
        
        Return:
        None -- Creates file with this HTML document.
        """

        # Add footer - finalize document
        self.create_footer()

        # Writing document to file
        f = open(file_name,"w")
        for e in self.document: f.write(e+"\n")
        f.close()

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

    @staticmethod
    def overlap(t1, t2):
        """Checks if one interval contains any overlap with another interval.

        Keyword arguments:
        t1 -- First tuple.
        t2 -- Second tuple.
  
        Returns:
        Returns -1 if i1 is before i2; 1 if i1 is after i2; and 0 if there is any overlap.
        """
        if(t1[1] <= t2[0]): return -1 # interval1 is before interval2
        if(t2[1] <= t1[0]): return 1 # interval1 is after interval2
        return 0 # interval1 overlaps interval2

