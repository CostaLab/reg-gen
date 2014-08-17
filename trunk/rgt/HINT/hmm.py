
###################################################################################################
# Libraries
###################################################################################################

# Python
import warnings
warnings.filterwarnings("ignore")

# Internal
from .. Util import ErrorHandler

###################################################################################################
# Classes
###################################################################################################

class HMM:
    """
    Represent an HMM

    Authors: Eduardo G. Gusmao.

    Methods:

    load_hmm(input_file_name):
    Loads an HMM based on a .hmm file.
    """

    def __init__(self):
        """ 
        Initializes HMM.

        Variables:
        states -- Number of states of the HMM (integer).
        dim -- Number of dimensions of the HMM (integer).
        pi -- Initial state probabilities vector (list).
        A -- Transition matrix (list of lists).
        means -- Matrix containing the mean vector of each state (list of lists).
        covs -- List of covariance matrices (list of list of lists).
        """
        self.states = 0
        self.dim = 0
        self.pi = []
        self.A = []
        self.means = []
        self.covs = []

    def load_hmm(self, input_file_name):
        """ 
        Loads all objects of this class based on input_file_name.
        These objects are used to create a scikit-based HMM.

        Keyword arguments:
        input_file_name -- File location + name in HMM format. See manual for full description of such format.
        
        Return:
        None -- It loads variables of HMM with information from input_file_name
        """

        # Creating file
        input_file = open(input_file_name,"r")

        # Reading number of states
        hmm_states = int(input_file.readline().strip().split(" ")[1])
        input_file.readline()

        # Reading initial state probabilities
        pi = [float(e) for e in input_file.readline().strip().split(" ")]
        input_file.readline()

        # Reading transition matrix
        A = []
        for i in range(0,hmm_states): A.append([float(e) for e in input_file.readline().strip().split(" ")])
        input_file.readline()
        first_emis_line = input_file.readline()

        # Reading emission probabilities of multivariate HMM
        E = []
        for i in range(0,hmm_states): 
            if(i==0): Elist = first_emis_line.strip().split("#")
            else: Elist = input_file.readline().strip().split("#")
            E.append([[float(e) for e in Elist[0].strip().split(" ")],[float(e) for e in Elist[1].strip().split(" ")]])
        dim_nb = len(E[0][0])

        # Closing file
        input_file.close()

        # Preparing scikit structures
        means_vec = [e[0] for e in E]
        cov_vec = []
        for e in E:
            new_mat = []
            for i in range(0,dim_nb):
                new_vec = []
                for j in range(0,dim_nb):
                    new_vec.append(e[1][(dim_nb*i)+j])
                new_mat.append(new_vec)
            cov_vec.append(new_mat)

        # Create scikit HMM
        self.states = hmm_states
        self.dim = dim_nb
        self.pi = pi
        self.A = A
        self.means = means_vec
        self.covs = cov_vec


