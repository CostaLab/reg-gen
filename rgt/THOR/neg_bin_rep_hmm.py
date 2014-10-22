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
import mpmath
from scipy.stats import rv_discrete

def get_init_parameters(s1, s2, **info):
    
    #tmp = sum( [ first_overall_coverage[i] + second_overall_coverage[i] for i in indices_of_interest]) / 2
    n_ = np.array([info['count'], info['count']])
    #print('n_: ', n_, file=sys.stderr)
    
    #_, s1, s2 = _get_training_sets(indices_of_interest, first_overall_coverage, second_overall_coverage, name, verbose, x, threshold, diff_cov)
    
    #get observation that occurs most often:
    m_ =[float(np.argmax(np.bincount(map(lambda x: x[0], s1)))), float(np.argmax(np.bincount(map(lambda x: x[1], s2)))) ]
    #print('m_', m_, file=sys.stderr)
    
    p_ = [[-1,-1,-1],[-1,-1,-1]] #first: 1. or 2. emission, second: state
    
    p_[0][0] = 1. / n_[0]
    p_[1][0] = 1. / n_[1]
       
    p_[0][1] = m_[0] / n_[0]
    p_[1][1] = p_[1][0]
    
    p_[0][2] = p_[0][0]
    p_[1][2] = m_[1] / n_[1]
    
    #print('p_', p_, file=sys.stderr)
    
    return n_, p_

class NegBin():
    def __init__(self, alpha, mu):
        nbin_mpmath = lambda k, p, r: mpmath.gamma(k + r)/(mpmath.gamma(k+1)*mpmath.gamma(r))*np.power(1-p, r)*np.power(p, k)
        self.nbin = np.frompyfunc(nbin_mpmath, 3, 1)
        self.p = (alpha * mu) /  float((1 + alpha * mu))
        self.r = 1./alpha
        c = 5000
        self.dist = rv_discrete(values=([i for i in range(c)], map(lambda x: float(x), [self.pdf(i) for i in range(c)])), name='dist')

    def pdf(self, k):
        return self.nbin(k, self.p, self.r)
    
    def logpdf(self, k):
        return log(self.nbin(k, self.p, self.r))
    
    def rvs(self, size=1):
        return self.dist.rvs(size=size)
    
    
    
class NegBinRepHMM(_BaseHMM):
    def __init__(self, alpha, mu, dim_cond_1, dim_cond_2, init_state_seq=None, n_components=3, covariance_type='diag', startprob=[1, 0, 0],
                 transmat=None, startprob_prior=None, transmat_prior=None,
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
        self.n_features = 2 #emission dimension
        self.alpha = alpha
        self.mu = mu
        raw1 = [NegBin(self.alpha[0, 0], self.mu[0, 0]), NegBin(self.alpha[0, 1], self.mu[0, 1]), NegBin(self.alpha[0, 2], self.mu[0, 2])]
        raw2 = [NegBin(self.alpha[1, 0], self.mu[1, 0]), NegBin(self.alpha[1, 1], self.mu[1, 1]), NegBin(self.alpha[1, 2], self.mu[1, 2])]
        
        self.neg_distr = np.matrix([raw1, raw2]) #matrix of all Neg. Bin. Distributions, columns=HMM's state (3), row=#samples (2)
        
#         self.init_state_seq = init_state_seq
#         self.count_s1, self.count_s2 = 0, 0

    def _get_emissionprob(self):
        return self.p
    
    def _compute_log_likelihood(self, X):
        #t = time()
        matrix = []
        lookup = {}
        for x in X: #over all observations
            row = []
            for i in range(self.n_components): #over number of HMM's state
                r_sum = 0
                for j in range(2): #over all samples
                    it = range(self.dim[0]) if j == 0 else range(self.dim[0], self.dim[0] + self.dim[1]) #grab proper ob
                    for k in it:
                        index = (int(x[k]), i, j)
                        if lookup.has_key( index ):
                            r_sum += lookup[index]
                        else:
                            y = self.neg_distr[j,i].logpdf(x[k])
                            lookup[index] = y
                            r_sum += y
                row.append(r_sum)
        
            matrix.append(row)
#         print('time to compute log-likelihood matrix (compute_log_likelihood): ',\
#               time()-t, file=sys.stderr)
        return np.asarray(matrix)
    
    
    def _generate_sample_from_state(self, state, random_state=None):
        output = []
        for i, d in enumerate(self.dim):
            for _ in range(d):
                output.append( self.neg_distr[i,state].rvs() )
        
        return np.array(output)
    
    def _initialize_sufficient_statistics(self):
        stats = super(NegBinRepHMM, self)._initialize_sufficient_statistics()
        stats['post'] = np.zeros(self.n_components)
        stats['post_emission'] = np.zeros([self.n_features, self.n_components])
        
        return stats
    
    def _help_accumulate_sufficient_statistics(self, obs, stats, posteriors):
        #posteriors = _valid_posteriors(posteriors, obs)
        for t, symbol in enumerate(obs):
            stats['post'] += posteriors[t]
            pot_it = [range(self.dim[0]), range(self.dim[0], self.dim[0] + self.dim[1])]
            for j, it in enumerate(pot_it):
                for i in it:
                    
                    stats['post_emission'][j] += posteriors[t] * symbol[i]
        #stats['obs'] = np.copy(obs)
        #stats['posterior'] = np.copy(posteriors)

    def _accumulate_sufficient_statistics(self, stats, obs, framelogprob,
                                      posteriors, fwdlattice, bwdlattice,
                                      params):
        super(NegBinRepHMM, self)._accumulate_sufficient_statistics(
            stats, obs, framelogprob, posteriors, fwdlattice, bwdlattice,
            params)
        
        self._help_accumulate_sufficient_statistics(obs, stats, posteriors)        
    
    def _help_do_mstep(self, stats):
        #add pseudo counts for nan entries
#        help_denum = _add_pseudo_counts( stats['post'] )
#        self.p[0] = stats['post_emission'][0] / (self.n[0] * help_denum)
#        self.p[1] = stats['post_emission'][1] / (self.n[1] * help_denum)
#        self.p[0] = _add_pseudo_counts(self.p[0])
#        self.p[1] = _add_pseudo_counts(self.p[1])
#        self.merge_distr()
#        print('m-step: ',self.p, file=sys.stderr)
#        tmp=np.array(map(lambda x: x*self.n[0], self.p))
#        print(np.reshape(tmp, (-1,3)), file=sys.stderr)
        for i in range(self.n_features):
            self.mu[i] = stats['post_emission'][i] / stats['post'][i]
        print(self.mu, file=sys.stderr)
    
    def _count(self, posts):
        c_1, c_2 = 0, 0
        
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
#        self.count_s1, self.count_s2 = self._count(stats['posterior'])
#        self._help_do_mstep(stats)
        
        
       
    def merge_distr(self):
        f = self.count_s2 / float(self.count_s1 + self.count_s2)
        p_high = self.p[0][1] + f * fabs(self.p[0][1] - self.p[1][2])
        #print('merge: ', f, self.p, p_high, file=sys.stderr)
        
        self.p[0][1] = p_high
        self.p[1][2] = p_high
        
        p_low = self.p[1][1] + f * fabs(self.p[1][1] - self.p[0][2])
        self.p[1][1] = p_low
        self.p[0][2] = p_low
        #print('merge: ', f, self.p, p_high, p_low, file=sys.stderr)
        tmp=np.array(map(lambda x: x*self.n[0], self.p))
        #print(np.reshape(tmp, (-1,3)), file=sys.stderr)
    
#    def _init(self, obs, params):
#        _init(self, obs, params)



if __name__ == '__main__':
    from numpy.random import negative_binomial
    
    alpha = np.matrix([[10.,10.,10.], [10.,10.,10.]])
    mu = np.matrix([[10.,20.,30.], [40.,50.,60.]])
    dim_cond_1 = 5
    dim_cond_2 = 4
    
#    transmat_ = np.array(([[0.7, 0.2, 0.1], [0.1, 0.8, 0.1], [0.2,0.2,0.6]]))
#     tmp = [[0.000001, 0.0000098, 0.000001], [0.000001, 0.000001, 0.0000098]]
    #[[0.01, 0.98, 0.01], [0.01, 0.01, 0.98]]
    
    r = 1 / alpha[0, 0]
    p = 1 / (1 + alpha[0,0] * mu[0,0])
    
    for _ in range(10):
        print(negative_binomial)
    
    
    #m = NegBinRepHMM(alpha = alpha, mu = mu, dim_cond_1 = dim_cond_1, dim_cond_2 = dim_cond_2)
    
    #X, Z = m.sample(10)
    #print(X)
    
#     X, Z = m.sample(40) #returns (obs, hidden_states)
    
#    X = np.array([[80,14], [34,92], [15,95],[15,5],[44,2]])
#    n_ = [ sum([x[i] for x in X]) for i in range(2) ]
#    X=np.array([[46,5],[41, 3],[43,4],[43,2],[45,4],[39,3],[18,36],[28,28],[43,1],[23,35]])
#     print('obs           ', X)
#     print('hidden states ', Z)
# #    X = np.array([[12,2],[11, 5],[12,4],[10,2],[4,4],[3,3],[2,1],[2,14],[4,11],[2,9]])
    #m2 = NegBinRepHMM(alpha = alpha, mu = np.matrix([[15.,15.,15.], [14.,16.,12.]]), dim_cond_1 = dim_cond_1, dim_cond_2 = dim_cond_2)
    #m2.fit([X])
    
    #print('states        ', Z)
    #print('estim. states ', m2.predict(X))
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
    