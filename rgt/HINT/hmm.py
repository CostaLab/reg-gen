import numpy as np
from scipy import linalg
from hmmlearn.hmm import GaussianHMM
from hmmlearn.base import ConvergenceMonitor
from sklearn.utils import check_array
from hmmlearn.utils import iter_from_X_lengths

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
        input_file_name -- File location + name in HMM format.
        See manual for full description of such format.
        
        Return:
        None -- It loads variables of HMM with information from input_file_name
        """

        # Creating file
        input_file = open(input_file_name, "r")

        # Reading number of states
        hmm_states = int(input_file.readline().strip().split(" ")[1])
        input_file.readline()

        # Reading initial state probabilities
        pi = [float(e) for e in input_file.readline().strip().split(" ")]
        input_file.readline()

        # Reading transition matrix
        A = []
        for i in range(0, hmm_states): A.append([float(e) for e in input_file.readline().strip().split(" ")])
        input_file.readline()
        first_emis_line = input_file.readline()

        # Reading emission probabilities of multivariate HMM
        E = []
        for i in range(0, hmm_states):
            if i == 0:
                Elist = first_emis_line.strip().split("#")
            else:
                Elist = input_file.readline().strip().split("#")
            E.append([[float(e) for e in Elist[0].strip().split(" ")], [float(e) for e in Elist[1].strip().split(" ")]])
        dim_nb = len(E[0][0])

        # Closing file
        input_file.close()

        # Preparing scikit structures
        means_vec = [e[0] for e in E]
        cov_vec = []
        for e in E:
            new_mat = []
            for i in range(0, dim_nb):
                new_vec = []
                for j in range(0, dim_nb):
                    new_vec.append(e[1][(dim_nb * i) + j])
                new_mat.append(new_vec)
            cov_vec.append(new_mat)

        # Create scikit HMM
        self.states = hmm_states
        self.dim = dim_nb
        self.pi = pi
        self.A = A
        self.means = means_vec
        self.covs = cov_vec

    def save_hmm(self, output_file_name, precision=6):
        """
        Save all objects of this class to the output_file_name.

        Keyword arguments:
        output_file_name -- File location + name in HMM format.

        Return:
        None -- It saves variables of HMM with information to output_file_name
        """
        output_file_name += ".hmm"
        with open(output_file_name, "w") as output_file:
            output_file.write("states " + str(self.states) + "\n")

            output_file.write("initial" + "\n")
            output_file.write(str(self.pi[0]))
            for e in self.pi[1:]:
                output_file.write(" " + str(e))

            output_file.write("\n" + "transitions" + "\n")
            for e in self.A:
                output_file.write(str(round(e[0], precision)))
                for v in e[1:]:
                    output_file.write(" " + str(round(v, precision)))
                output_file.write("\n")

            output_file.write("emissions" + "\n")
            for idx in range(self.states):
                output_file.write(str(round(self.means[idx][0], precision)))
                for e in self.means[idx][1:]:
                    output_file.write(" " + str(round(e, precision)))
                output_file.write(" # ")
                output_file.write(str(round(self.covs[idx][0], precision)))
                for e in self.covs[idx][1:]:
                    output_file.write(" " + str(round(e, precision)))
                output_file.write("\n")


class SemiSupervisedGaussianHMM(GaussianHMM):
    def __init__(self, n_components=1, covariance_type='diag', min_covar=1e-3, startprob_prior=1.0,
                 transmat_prior=1.0, means_prior=0, means_weight=0, covars_prior=1e-2, covars_weight=1,
                 algorithm="viterbi", random_state=None, n_iter=5, tol=1e-2, verbose=False,
                 params="stmc", init_params="stmc", states_prior=None, fp_state=None):
        GaussianHMM.__init__(self, n_components=n_components, covariance_type=covariance_type,
                             min_covar=min_covar, startprob_prior=startprob_prior, transmat_prior=transmat_prior,
                             means_prior=means_prior, means_weight=means_weight,
                             covars_prior=covars_prior, covars_weight=covars_weight,
                             algorithm=algorithm, random_state=random_state,
                             n_iter=n_iter, tol=tol, verbose=verbose,
                             params=params, init_params=init_params)

        self.covariance_type = covariance_type
        self.min_covar = min_covar
        self.means_prior = means_prior
        self.means_weight = means_weight
        self.covars_prior = covars_prior
        self.covars_weight = covars_weight
        self.states_prior = states_prior
        self.fp_state = fp_state

    def fit(self, X, lengths=None):
        """Estimate model parameters.
        An initialization step is performed before entering the
        EM algorithm. If you want to avoid this step for a subset of
        the parameters, pass proper ``init_params`` keyword argument
        to estimator's constructor.
        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Feature matrix of individual samples.
        lengths : array-like of integers, shape (n_sequences, )
            Lengths of the individual sequences in ``X``. The sum of
            these should be ``n_samples``.
        Returns
        -------
        self : object
            Returns self.
        """
        X = check_array(X)
        self._init(X, lengths=lengths)
        self._check()

        self.monitor_ = ConvergenceMonitor(self.tol, self.n_iter, self.verbose)
        for iter in range(self.n_iter):
            stats = self._initialize_sufficient_statistics()
            curr_logprob = 0
            for i, j in iter_from_X_lengths(X, lengths):
                framelogprob = self._compute_log_likelihood(X[i:j])
                logprob, fwdlattice = self._do_forward_pass(framelogprob)
                curr_logprob += logprob
                bwdlattice = self._do_backward_pass(framelogprob)
                posteriors = self._compute_posteriors(fwdlattice, bwdlattice)

                # fix posteriors
                if self.states_prior is not None and self.fp_state is not None:
                    for k in range(len(self.states_prior)):
                        if self.states_prior[k] == 0:
                            # non footprint states
                            posteriors[k][self.fp_state] = 0.0
                            posteriors[k] = posteriors[k] / sum(posteriors[k])

                        elif self.states_prior[k] == 1:
                            # footprint states
                            posteriors[k] = 0.0 / sum(posteriors[k])
                            posteriors[k][self.fp_state] = 1.0

                self._accumulate_sufficient_statistics(stats, X[i:j], framelogprob, posteriors, fwdlattice, bwdlattice)

            self._do_mstep(stats)

            self.monitor_.report(curr_logprob)
            if self.monitor_.converged:
                break

        return self


def _compute_log_likelihood(self, X):
        return log_multivariate_normal_density(
            X, self.means_, self._covars_, self.covariance_type)

def log_multivariate_normal_density(X, means, covars, covariance_type='diag'):
    """Compute the log probability under a multivariate Gaussian distribution.
    Parameters
    ----------
    X : array_like, shape (n_samples, n_features)
        List of n_features-dimensional data points. Each row corresponds to a
        single data point.
    means : array_like, shape (n_components, n_features)
        List of n_features-dimensional mean vectors for n_components Gaussians.
        Each row corresponds to a single mean vector.
    covars : array_like
        List of n_components covariance parameters for each Gaussian. The shape
        depends on `covariance_type`:
            (n_components, n_features)      if 'spherical',
            (n_features, n_features)    if 'tied',
            (n_components, n_features)    if 'diag',
            (n_components, n_features, n_features) if 'full'
    covariance_type : string
        Type of the covariance parameters.  Must be one of
        'spherical', 'tied', 'diag', 'full'.  Defaults to 'diag'.
    Returns
    -------
    lpr : array_like, shape (n_samples, n_components)
        Array containing the log probabilities of each data point in
        X under each of the n_components multivariate Gaussian distributions.
    """
    log_multivariate_normal_density_dict = {
        'spherical': _log_multivariate_normal_density_spherical,
        'tied': _log_multivariate_normal_density_tied,
        'diag': _log_multivariate_normal_density_diag,
        'full': _log_multivariate_normal_density_full}
    return log_multivariate_normal_density_dict[covariance_type](
        X, means, covars)

def _log_multivariate_normal_density_diag(X, means, covars):
    """Compute Gaussian log-density at X for a diagonal model."""
    n_samples, n_dim = X.shape
    lpr = -0.5 * (n_dim * np.log(2 * np.pi) + np.sum(np.log(covars), 1)
                  + np.sum((np.array(means) ** 2) / covars, 1)
                  - 2 * np.dot(X, (means / covars).T)
                  + np.dot(np.array(X) ** 2, (1.0 / covars).T))
    return lpr

def _log_multivariate_normal_density_spherical(X, means, covars):
    """Compute Gaussian log-density at X for a spherical model."""
    cv = covars.copy()
    if covars.ndim == 1:
        cv = cv[:, np.newaxis]
    if cv.shape[1] == 1:
        cv = np.tile(cv, (1, X.shape[-1]))
    return _log_multivariate_normal_density_diag(X, means, cv)


def _log_multivariate_normal_density_tied(X, means, covars):
    """Compute Gaussian log-density at X for a tied model."""
    cv = np.tile(covars, (means.shape[0], 1, 1))
    return _log_multivariate_normal_density_full(X, means, cv)


def _log_multivariate_normal_density_full(X, means, covars, min_covar=1.e-7):
    """Log probability for full covariance matrices."""
    n_samples, n_dim = X.shape
    nmix = len(means)
    log_prob = np.empty((n_samples, nmix))
    for c, (mu, cv) in enumerate(zip(means, covars)):
        try:
            cv_chol = linalg.cholesky(cv, lower=True)
        except linalg.LinAlgError:
            # The model is most probably stuck in a component with too
            # few observations, we need to reinitialize this components
            try:
                cv_chol = linalg.cholesky(cv + min_covar * np.eye(n_dim),
                                          lower=True)
            except linalg.LinAlgError:
                raise ValueError("'covars' must be symmetric, "
                                 "positive-definite")

        cv_log_det = 2 * np.sum(np.log(np.diagonal(cv_chol)))
        cv_sol = linalg.solve_triangular(cv_chol, (X - mu).T, lower=True).T
        log_prob[:, c] = - .5 * (np.sum(np.array(cv_sol) ** 2, axis=1) +
                                 n_dim * np.log(2 * np.pi) + cv_log_det)

    return log_prob
