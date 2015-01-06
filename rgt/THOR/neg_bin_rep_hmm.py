#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on Mar 6, 2013

@author: manuel
'''
from __future__ import print_function
from scipy.stats import binom
from sklearn.hmm import _BaseHMM
import string
import numpy as np
from sklearn import hmm
from sklearn.hmm import _BaseHMM
import sys
from time import time
from math import fabs, log
from sklearn.utils.extmath import logsumexp
# from help_hmm import _init, _add_pseudo_counts, _valid_posteriors

from scipy.stats import nbinom
from random import randint
from scipy.special import gamma
from operator import mul

import scipy.special as special
import scipy.optimize as optimize
import numpy as np
from neg_bin import NegBin

import warnings
warnings.filterwarnings('error')


def get_init_parameters(s0, s1, s2, **info):
    """For given training set (s0: Background, s1: Gaining, s2: loseing) get inital mu, alpha for NB1."""
    mu = np.matrix([np.mean(map(lambda x: x[i], s)) for i in range(2) for s in [s0, s1, s2]]).reshape(2, 3)
    var = np.matrix([np.var(map(lambda x: x[i], s)) for i in range(2) for s in [s0, s1, s2]]).reshape(2, 3)
    
    alpha = (var - mu) / np.square(mu)
    alpha[alpha < 0] = 1e-300
    
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
    def __init__(self, alpha, mu, dim_cond_1, dim_cond_2, init_state_seq=None, n_components=3, covariance_type='diag', startprob=[1, 0, 0],
                 transmat=None, startprob_prior=None, transmat_prior=None, func=None,
                 algorithm="viterbi", means_prior=None, means_weight=0,
                 covars_prior=1e-2, covars_weight=1,
                 random_state=None, n_iter=10, thresh=1e-2,
                 params=string.ascii_letters,
                 init_params=string.ascii_letters):
    
        _BaseHMM.__init__(self, n_components, startprob, transmat,
                          startprob_prior=startprob_prior,
                          transmat_prior=transmat_prior, algorithm=algorithm,
                          random_state=random_state, n_iter=n_iter,
                          thresh=thresh, params=params,
                          init_params=init_params)
        
        self.dim = [dim_cond_1, dim_cond_2] #dimension of one emission
        self.n_features = 2 #sum(self.dim) #emission dimension
        self.alpha = alpha
        self.mu = mu
        self._update_distr(self.mu, self.alpha)
        self.func = func
        
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
                    it = range(self.dim[0]) if j == 0 else range(self.dim[0], self.dim[0] + self.dim[1]) #grab proper ob
                    for k in it:
                        index = (int(x[k]), i, j)
                        if lookup.has_key( index ):
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
        stats['post_emission'] = np.zeros([self.n_features, self.n_components])
        
        return stats
    
    def _help_accumulate_sufficient_statistics(self, obs, stats, posteriors):
        for t, symbol in enumerate(obs):
            stats['post'][0] += posteriors[t]
            stats['post'][1] += posteriors[t]
            
            pot_it = [range(self.dim[0]), range(self.dim[0], self.dim[0] + self.dim[1])] #consider both classes
            for j, it in enumerate(pot_it):
                for i in it:
                    stats['post_emission'][j] += posteriors[t] * symbol[i]
        
        stats['post'][0] = stats['post'][0]*self.dim[0]
        stats['post'][1] = stats['post'][1]*self.dim[1]
        
        stats['posterior'] = np.copy(posteriors)
        

    def _accumulate_sufficient_statistics(self, stats, obs, framelogprob,
                                      posteriors, fwdlattice, bwdlattice,
                                      params):
        super(NegBinRepHMM, self)._accumulate_sufficient_statistics(
            stats, obs, framelogprob, posteriors, fwdlattice, bwdlattice,
            params)
        
        self._help_accumulate_sufficient_statistics(obs, stats, posteriors)        
    
    def _help_do_mstep(self, stats):
        for i in range(self.n_features):
            print('help_m_step', 'i', stats['post_emission'][i], stats['post'][i], file=sys.stderr)
            #print('Perform M-step of EM algorithm (max. 20 times)', file=sys.stderr)
            self.mu[i] = stats['post_emission'][i] / stats['post'][i]
        
        tmp_a = [map(lambda m: self.get_alpha(m), np.asarray(self.mu[i])[0]) for i in range(self.n_features)]
        self.alpha = np.matrix(tmp_a)
        self._update_distr(self.mu, self.alpha)
    
    def _count(self, posts):
        c_1, c_2 = 1, 1
        
        for s0, s1, s2 in posts:        
            if s0 > 0.5:
                c_1 += 0
                c_2 += 0
            elif s1 >= s2:
                c_1 += 1
            elif s2 > s1:
                c_2 += 1
        return c_1, c_2
        
    
    def _do_mstep(self, stats, params):
        super(NegBinRepHMM, self)._do_mstep(stats, params)
        self._help_do_mstep(stats)
        self.count_s1, self.count_s2 = self._count(stats['posterior'])
        self.merge_distr()
       
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
            el[0,0] = med
            el[1,0] = med
        
        self._update_distr(self.mu, self.alpha)

    
if __name__ == '__main__':
    alpha = np.matrix([[0.2, 0.2, 0.2], [0.2, 0.2, 0.2]])
    mu = np.matrix([[10.,100.,10.], [10.,10.,100.]])
    para_func = [[1, 2, 0.3], [1, 2, 0.3]]
    para_func = [np.array([-0.02178527,  0.48686578,  0.21833156]), np.array([ 0.76335214, -0.94956275,  0.70959764])]
    
    #tmp_a = []
    #for i in range(2):
        #tmp_a.append(map(lambda m: get_alpha(para_func, i, m), np.asarray(mu[i])[0]))
    
    #alpha = np.matrix(tmp_a)
    #print(mu)
    #print(alpha)
    
    dim_cond_1 = 5
    dim_cond_2 = 5

    m = NegBinRepHMM(alpha = alpha, mu = mu, dim_cond_1 = dim_cond_1, dim_cond_2 = dim_cond_2)
    
    X, Z = m.sample(20)
#     for i, el in enumerate(X):
#         print(el, Z[i], sep='\t')
    
    
#    X = np.array([[80,14], [34,92], [15,95],[15,5],[44,2]])
#    n_ = [ sum([x[i] for x in X]) for i in range(2) ]
#    X=np.array([[46,5],[41, 3],[43,4],[43,2],[45,4],[39,3],[18,36],[28,28],[43,1],[23,35]])
#     print('obs           ', X)
#     print('hidden states ', Z)
# #    X = np.array([[12,2],[11, 5],[12,4],[10,2],[4,4],[3,3],[2,1],[2,14],[4,11],[2,9]])
    
    m2 = NegBinRepHMM(alpha = alpha, mu = np.matrix([[50.,130.,110.], [60.,100.,120.]]), dim_cond_1 = dim_cond_1, dim_cond_2 = dim_cond_2, para_func = para_func)
    m2.fit([X])
    
    posteriors = m.predict_proba(X)
    e = m2.predict(X)
    for i, el in enumerate(X):
        print(el, Z[i], e[i], Z[i] == e[i], sep='\t', file=sys.stderr)
    
    #print(m2.transmat_)
    #print(m2.alpha)
    #print(m2.mu)
    #print(m2.predict_proba(X))
#     
#     
# #    logprob, posteriors = m2.eval(X)
# #    print('logprob:', logprob)
# #    print('posteriors:', posteriors)
#     
#     print('estim. states ', m2.predict(X))
#     print(m2.predict_proba(X))
#     print(m2.n)
#     print(m2.p)
#     print(m2._get_transmat())
#     init_state = m2.predict(X)
#     m3 = BinomialHMM2d3s(n_components=3, n=n_)
#     m3.fit([X], init_params='advanced')
#     print(m3._get_transmat())
#     print(m3.p)
#     m2.eval(X)
    