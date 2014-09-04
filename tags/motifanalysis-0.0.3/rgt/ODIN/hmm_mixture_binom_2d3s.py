#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on Mar 6, 2013

@author: manuel
'''
from __future__ import print_function
from scipy.stats import binom
from sklearn.hmm import _BaseHMM
import string , numpy as np
from sklearn import hmm
from sklearn.hmm import _BaseHMM
import sys
from time import time
from math import fabs
from sklearn.utils.extmath import logsumexp
from help_hmm import _init, _add_pseudo_counts, _valid_posteriors

class BinomialHMM2d3s(_BaseHMM):
    def __init__(self, n, init_state_seq=None, p = [[[0.1, 0.1], [0.2, 0.9], [0.9, 0.2]], [[0.1, 0.1], [0.2, 0.9], [0.9, 0.2]]], a=[0.5, 0.5], n_components=2, covariance_type='diag', startprob=None,
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
        self.a = a
        self.n = n
        self.p = p
        self.n_features = 2 #emission dimension
        self.init_state_seq = init_state_seq
        self.count_s1, self.count_s2 = 0, 0
        self.distr_magnitude = 2 #

    def _get_emissionprob(self):
        return self.p
    
    def _compute_log_likelihood(self, X):
        t = time()
        matrix = []
        lookup = {}
        k = 0
        for x in X:
            row = []
            for i in range(self.n_components): #state
                sum = 0
                for j in range(self.n_features): #dim
                    for l in range(self.distr_magnitude):
                        index = (x[j], self.n[j], self.p[j][i][l])
                        if lookup.has_key( index ):
                            sum += lookup[index] * self.a[l]
                            k += 1
                        else:
                            y = binom.logpmf(x[j], self.n[j], self.p[j][i][l])
                            lookup[index] = y
                            sum += y * self.a[l]
                row.append(sum)
                
            matrix.append(row)
#         print('time to compute log-likelihood matrix (compute_log_likelihood): ',\
#               time()-t, file=sys.stderr)
        #print("HALLO")
        #print(np.asarray(matrix))
        return np.asarray(matrix)

    def _generate_sample_from_state(self, state, random_state=None):
        return np.array( map(lambda x: round(x), [sum([binom.rvs(self.n[i], self.p[i][state][j]) * self.a[i] for i in range(self.distr_magnitude)]) for j in range(self.n_features)]) )
    
    def _initialize_sufficient_statistics(self):
        stats = super(BinomialHMM2d3s, self)._initialize_sufficient_statistics()
        stats['post'] = np.zeros(self.n_components)
        stats['post_emission'] = np.zeros([self.n_features, self.n_components])
        return stats
    
    def _help_accumulate_sufficient_statistics(self, obs, stats, posteriors):
        posteriors = _valid_posteriors(posteriors, obs)
        for t, symbol in enumerate(obs):
            stats['post'] += posteriors[t]
            stats['post_emission'][0] += (posteriors[t] * symbol[0])
            stats['post_emission'][1] += (posteriors[t] * symbol[1])
        #stats['obs'] = np.copy(obs)
        stats['posterior'] = np.copy(posteriors)

    def _accumulate_sufficient_statistics(self, stats, obs, framelogprob,
                                      posteriors, fwdlattice, bwdlattice,
                                      params):
        super(BinomialHMM2d3s, self)._accumulate_sufficient_statistics(
            stats, obs, framelogprob, posteriors, fwdlattice, bwdlattice,
            params)
        
        self._help_accumulate_sufficient_statistics(obs, stats, posteriors)        
    
    def _help_do_mstep(self, stats):
        #add pseudo counts for nan entries
        help_denum = _add_pseudo_counts( stats['post'] )
        
        self.p[0] = stats['post_emission'][0] / (self.n[0] * help_denum)
        self.p[1] = stats['post_emission'][1] / (self.n[1] * help_denum)
        
        self.p[0] = _add_pseudo_counts(self.p[0])
        self.p[1] = _add_pseudo_counts(self.p[1])
        
        self.merge_distr()
        #print('m-step: ',self.p, file=sys.stderr)
        #tmp=np.array(map(lambda x: x*self.n[0], self.p))
        #print(np.reshape(tmp, (-1,3)), file=sys.stderr)
        
    
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
        super(BinomialHMM2d3s, self)._do_mstep(stats, params)
        self.count_s1, self.count_s2 = self._count(stats['posterior'])
        self._help_do_mstep(stats)
        
        
       
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
#    transmat_ = np.array(([[0.7, 0.2, 0.1], [0.1, 0.8, 0.1], [0.2,0.2,0.6]]))
    tmp = [[0.000001, 0.0000098, 0.000001], [0.000001, 0.000001, 0.0000098]]
    tmp = p = [[[0.3, 0.2], [0.6, 0.8], [0.7, 0.8]], [[0.2, 0.2], [0.6, 0.8], [0.7, 0.8]]]
    #[[0.01, 0.98, 0.01], [0.01, 0.01, 0.98]]
    p_ = np.array(tmp)
    n_ = np.array([2000000, 2000000])
    
    m = BinomialHMM2d3s(n_components=3, startprob=[1,0,0], n = [100, 100])

    X, Z = m.sample(10) #returns (obs, hidden_states)
    
#    X = np.array([[80,14], [34,92], [15,95],[15,5],[44,2]])
#    n_ = [ sum([x[i] for x in X]) for i in range(2) ]
#    X=np.array([[46,5],[41, 3],[43,4],[43,2],[45,4],[39,3],[18,36],[28,28],[43,1],[23,35]])
    print('obs           ', X)
    print('hidden states ', Z)
#    X = np.array([[12,2],[11, 5],[12,4],[10,2],[4,4],[3,3],[2,1],[2,14],[4,11],[2,9]])
    m2 = BinomialHMM2d3s(n_components=3, n=n_, p=tmp)
    m2.fit([X])
    
    
#    logprob, posteriors = m2.eval(X)
#    print('logprob:', logprob)
#    print('posteriors:', posteriors)
    
    print('estim. states ', m2.predict(X))
#     print(m2.n)
#     print(m2.p)
#     print(m2._get_transmat())
#     init_state = m2.predict(X)
#     m3 = BinomialHMM2d3s(n_components=3, n=n_)
#     m3.fit([X], init_params='advanced')
#     print(m3._get_transmat())
#     print(m3.p)
#     m2.eval(X)
    