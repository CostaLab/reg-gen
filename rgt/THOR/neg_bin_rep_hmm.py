#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
THOR detects differential peaks in multiple ChIP-seq profiles associated
with two distinct biological conditions.

Copyright (C) 2014-2016 Manuel Allhoff (allhoff@aices.rwth-aachen.de)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

@author: Manuel Allhoff
"""


import string
from hmmlearn.hmm import _BaseHMM

import sys
from math import fabs
from scipy.special import logsumexp

import numpy as np
from .neg_bin import NegBin

import warnings


def _get_pvalue_distr(mu, alpha, tracker):
    """Derive NB1 parameters for p-value calculation"""
    mu = mu[0,0]
    alpha = alpha[0,0] / 10000.
    tracker.write(text=str(mu), header="Neg. Bin. distribution for p-value estimates (mu)")
    tracker.write(text=str(alpha), header="Neg. Bin. distribution for p-value estimates (alpha)")
    
    nb = NegBin(mu, alpha)
    return {'distr_name': 'nb', 'distr': nb}

def get_init_parameters(s0, s1, s2, **info):
    """For given training set (s0: Background, s1: Gaining, s2: loseing) get inital mu, alpha for NB1."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        mu = np.matrix([np.mean([x[i] for x in s]) for i in range(2) for s in [s0, s1, s2]]).reshape(2, 3)
        var = np.matrix([np.var([x[i] for x in s]) for i in range(2) for s in [s0, s1, s2]]).reshape(2, 3)

    alpha = (var - mu) / np.square(mu)

    # alpha[np.isnan(alpha)] = 0.001

    alpha[alpha < 0] = 0.001

    for el in [mu, alpha]:
        high = min(el[0,1], el[1,2]) + 0.5 * fabs(el[0,1] - el[1,2])
        low = min(el[1,1], el[0,2]) + 0.5 * fabs(el[1,1] - el[0,2])
        med = np.mean([el[0,0], el[1,0]])
        el[0,1] = high
        el[1,2] = high
        el[1,1] = low
        el[0,2] = low
        el[0,0] = med
        el[1,0] = med

    return alpha, mu
    
class NegBinRepHMM(_BaseHMM):
    def __init__(self, alpha, mu, dim_cond_1, dim_cond_2, init_state_seq=None, n_components=3, covariance_type='diag',
                  startprob_prior=1.0, transmat_prior=1.0, func=None,
                 algorithm="viterbi", means_prior=None, means_weight=0,
                 covars_prior=1e-2, covars_weight=1,
                 random_state=None, n_iter=30, thresh=1e-2,
                 params=string.ascii_letters,
                 init_params=string.ascii_letters):
    
        _BaseHMM.__init__(self, n_components,
                          startprob_prior=startprob_prior,
                          transmat_prior=transmat_prior, algorithm=algorithm,
                          random_state=random_state, n_iter=n_iter,
                          tol=thresh, params=params,
                          init_params=init_params)
        
        self.dim = [dim_cond_1, dim_cond_2] #dimension of one emission
        self.n_features = 2 #sum(self.dim) #emission dimension
        self.alpha = alpha
        self.mu = mu
        self._update_distr(self.mu, self.alpha)
        self.func = func
        self.em_prob = 0
    
    
    def fit(self, obs, three_para):
        """Estimate model parameters.

        An initialization step is performed before entering the EM
        algorithm. If you want to avoid this step, pass proper
        ``init_params`` keyword argument to estimator's constructor.

        Parameters
        ----------
        obs : list
            List of array-like observation sequences, each of which
            has shape (n_i, n_features), where n_i is the length of
            the i_th observation.

        Notes
        -----
        In general, `logprob` should be non-decreasing unless
        aggressive pruning is used.  Decreasing `logprob` is generally
        a sign of overfitting (e.g. a covariance parameter getting too
        small).  You can fix this by getting more training data,
        or strengthening the appropriate subclass-specific regularization
        parameter.
        """

        # what does this mean??
        self._init(obs, self.init_params)

        logprob = []
        for i in range(self.n_iter):
            # Expectation step
            stats = self._initialize_sufficient_statistics()
            curr_logprob = 0
            for seq in obs:
                framelogprob = self._compute_log_likelihood(seq)
                lpr, fwdlattice = self._do_forward_pass(framelogprob)
                bwdlattice = self._do_backward_pass(framelogprob)
                gamma = fwdlattice + bwdlattice
                posteriors = np.exp(gamma.T - logsumexp(gamma, axis=1)).T
                curr_logprob += lpr
                self._accumulate_sufficient_statistics(
                    stats, seq, framelogprob, posteriors, fwdlattice,
                    bwdlattice)
            logprob.append(curr_logprob)

            # Check for convergence.
            if i > 0 and logprob[-1] - logprob[-2] < self.tol:
                break

            # Maximization step
            self._do_mstep(stats, three_para)
        #print("Logprob of all M-steps: %s" %logprob, file=sys.stderr)
        self.em_prob = logprob[-1]
        return self
    
    def _update_distr(self, mu, alpha):
        """Update distributions assigned to each state with new mu and alpha"""
        raw1 = [NegBin(mu[0, 0], alpha[0, 0]), NegBin(mu[0, 1], alpha[0, 1]), NegBin(mu[0, 2], alpha[0, 2])]
        raw2 = [NegBin(mu[1, 0], alpha[1, 0]), NegBin(mu[1, 1], alpha[1, 1]), NegBin(mu[1, 2], alpha[1, 2])]
        
        self.neg_distr = np.matrix([raw1, raw2]) #matrix of all Neg. Bin. Distributions, columns=HMM's state (3), row=#samples (2)
        
    def get_alpha(self, m):
        """Return alpha for a given mu based on empirical variance"""
        var = self.func(m)
        try:
            return max((var - m) / m**2, 1e-300)
        except Warning:
            if m**2 > 1e-300:
                return max((var - m) / m**2, 1e-300)
            else:
                return 1e-300
    
    def _compute_log_likelihood(self, X):
        matrix = []
        lookup = {}
        for x in X: #over all observations
            row = []
            for i in range(self.n_components): #over number of HMM's state
                r_sum = 0
                for j in range(self.n_features): #over dim
                    it = list(range(self.dim[0])) if j == 0 else list(range(self.dim[0], self.dim[0] + self.dim[1])) #grab proper ob
                    for k in it:
                        index = (int(x[k]), i, j)
                        if index in lookup:
                            r_sum += lookup[index]
                        else:
                            y = float(self.neg_distr[j,i].logpdf(x[k]))
                            lookup[index] = y
                            r_sum += y
                row.append(r_sum)
        
            matrix.append(row)
        return np.asarray(matrix)
    
    
    def _generate_sample_from_state(self, state, random_state=None):
        output = []
        for i, d in enumerate(self.dim):
            for _ in range(d):
                output.append( self.neg_distr[i,state].rvs() )
        
        return np.array(output)
    
    def _initialize_sufficient_statistics(self):
        stats = super(NegBinRepHMM, self)._initialize_sufficient_statistics()
        stats['post'] = np.zeros([self.n_features, self.n_components])
        stats['post_emission'] = np.zeros([self.n_features, self.n_components]) #dim X states
        
        return stats
    
    def _help_accumulate_sufficient_statistics(self, obs, stats, posteriors):
        for t, symbol in enumerate(obs):
            stats['post'][0] += posteriors[t]
            stats['post'][1] += posteriors[t]
            
            pot_it = [list(range(self.dim[0])), list(range(self.dim[0], self.dim[0] + self.dim[1]))] #consider both classes
            for j, it in enumerate(pot_it):
                for i in it:
                    stats['post_emission'][j] += posteriors[t] * symbol[i]
        
        stats['post'][0] = stats['post'][0] * self.dim[0]
        stats['post'][1] = stats['post'][1] * self.dim[1]
        
        stats['posterior'] = np.copy(posteriors)
    
    def _valid_posteriors(self, posteriors, obs):
    
        warnings.filterwarnings('error')
        
        for i in range(len(obs)):
            state_1 = False
            c1, c2 = np.mean(obs[i][:self.dim[0]]), np.mean(obs[i][self.dim[0]:]) #counts of samples
    
            if posteriors[i][0] > 0.5: 
                state_1 = True
                
            if c1 > c2: #state 1
                if fabs(posteriors[i][2] - 1) < 1e-200:
                    posteriors[i] = np.array([1, 0, 0])
                else:
                    if not state_1 and posteriors[i][2] > posteriors[i][1]:
                        try:                        
                            post_s2 = 0
                            post_s0 = posteriors[i][0] / (posteriors[i][0] + posteriors[i][1])
                            post_s1 = posteriors[i][1] / (posteriors[i][0] + posteriors[i][1])
                            posteriors[i] = np.array([post_s0, post_s1, post_s2])
                        except RuntimeWarning:
                            print(posteriors[i], c1, c2, file=sys.stderr)
            
            if c2 > c1: #state 2
                if fabs(posteriors[i][1] - 1) < 1e-200:
                    posteriors[i] = np.array([1, 0, 0])
                else:
                    if not state_1 and posteriors[i][1] > posteriors[i][2]:
                        try:
                            post_s1 = 0
                            post_s0 = posteriors[i][0] / (posteriors[i][0] + posteriors[i][2])
                            post_s2 = posteriors[i][2] / (posteriors[i][0] + posteriors[i][2])
                            posteriors[i] = np.array([post_s0, post_s1, post_s2])
                        except RuntimeWarning:
                            print(posteriors[i], c1, c2, file=sys.stderr)

        warnings.resetwarnings()

        return posteriors

    def _accumulate_sufficient_statistics(self, stats, obs, framelogprob,
                                      posteriors, fwdlattice, bwdlattice
                                      ):
        super(NegBinRepHMM, self)._accumulate_sufficient_statistics(
            stats, obs, framelogprob, posteriors, fwdlattice, bwdlattice,
            )
        posteriors = self._valid_posteriors(posteriors, obs)
        self._help_accumulate_sufficient_statistics(obs, stats, posteriors)        
    
    def _do_mstep(self, stats, three_para):
        super(NegBinRepHMM, self)._do_mstep(stats)
        
        if three_para:
            self.mu[0,1] = (stats['post_emission'][0][1] + stats['post_emission'][1][2]) / (stats['post'][0][1] + stats['post'][1][2])
            self.mu[1,1] = (stats['post_emission'][1][1] + stats['post_emission'][0][2]) / (stats['post'][1][1] + stats['post'][0][2])
            self.mu[0,0] = (stats['post_emission'][0][0] + stats['post_emission'][1][0]) / (stats['post'][0][0] + stats['post'][1][0])
            
            self.mu[1,2] = self.mu[0,1]
            self.mu[0,2] = self.mu[1,1]
            self.mu[1,0] = self.mu[0,0]
        else:
            self.mu[0,1] = (stats['post_emission'][0][1] + stats['post_emission'][1][2]) / (stats['post'][0][1] + stats['post'][1][2])
            self.mu[1,1] = (stats['post_emission'][1][1] + stats['post_emission'][0][2] + stats['post_emission'][0][0] + stats['post_emission'][0][1]) / (stats['post'][1][1] + stats['post'][0][2] + stats['post'][0][0] + stats['post'][0][1])
            
            self.mu[0,0] = self.mu[1,1]
            self.mu[1,2] = self.mu[0,1]
            self.mu[0,2] = self.mu[1,1]
            self.mu[1,0] = self.mu[0,0]
        
        self.alpha = np.matrix([[self.get_alpha(m) for m in np.asarray(self.mu[i])[0]] for i in range(self.n_features)])
        self._update_distr(self.mu, self.alpha)
       
    def merge_distr(self):
        f = self.count_s2 / float(self.count_s1 + self.count_s2) #TODO exp_data.
        
        for el in [self.mu, self.alpha]:
            high = min(el[0,1], el[1,2]) + f * fabs(el[0,1] - el[1,2])
            low = min(el[1,1], el[0,2]) + f * fabs(el[1,1] - el[0,2])
            med = np.mean([el[0,0], el[1,0]])
            el[0,1] = high
            el[1,2] = high
            el[1,1] = low
            el[0,2] = low
            el[0,0] = low #min(med, low)
            el[1,0] = low #min(med, low)
        
        self._update_distr(self.mu, self.alpha)

    
if __name__ == '__main__':
    alpha = np.matrix([[0.2, 0.2, 0.2], [0.2, 0.2, 0.2]])
    mu = np.matrix([[15.,100.,10.], [10.,10.,100.]])
    f = lambda x: 0.4*x**2 + 1
    
    dim_cond_1 = 2
    dim_cond_2 = 3
    mean = [0, 0]
    cov = [[1, 0], [0, 100]]
    X= np.random.multivariate_normal(mean, cov, size=5)

    m2 = NegBinRepHMM(alpha = alpha, mu = np.matrix([[50.,130.,110.], [60.,100.,120.]]), dim_cond_1 = dim_cond_1, dim_cond_2 = dim_cond_2, func=f)
    m2.fit([X], three_para=False)
    
    posteriors = m2.predict_proba(X)
    e = m2.predict(X)
    for i, el in enumerate(X):
        print(el, e[i], sep='\t', file=sys.stderr)
    #print(np.max(posteriors, axis=1))
    print(m2.mu)
