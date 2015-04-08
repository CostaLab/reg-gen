#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on Mar 6, 2013

@author: manuel
'''
from __future__ import print_function
from scipy.stats import binom
from hmmlearn.hmm import _BaseHMM
import string , numpy as np
from math import fabs
from help_hmm import _init, _add_pseudo_counts, _valid_posteriors
# import cProfile

class BinomialHMM2d3s(_BaseHMM):
    def __init__(self, n, init_state_seq=None, p = [[[0.1, 0.1, 0.2], [0.8, 0.9, 0.7], [0.3, 0.2, 0.25]], [[0.1, 0.1, 0.12], [0.2, 0.2, 0.3], [0.9, 0.7, 0.8]]], \
                 c = [[[0.5, 0.4, 0.1], [0.5, 0.4, 0.1], [0.5, 0.4, 0.1]], [[0.5, 0.4, 0.1], [0.5, 0.4, 0.1], [0.5, 0.4, 0.1]]], n_components=3, covariance_type='diag', startprob=None,
                 transmat=None, startprob_prior=None, transmat_prior=None, distr_magnitude = 3,
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
        self.c = c
        self.n = n
        self.p = p # 1) emission-component 2) state 3) component
        self.n_features = 2 #emission dimension
        self.init_state_seq = init_state_seq
        self.count_s1, self.count_s2 = 0, 0
        self.distr_magnitude = distr_magnitude

    def _get_emissionprob(self):
        return self.p
    
    def _compute_log_likelihood(self, X):
        #t = time()
        matrix = []
        lookup = {}
        k = 0
        for x in X:
            row = []
            for state in range(self.n_components): #state
                res = 0
                for dim in range(self.n_features): #dim
                    for comp in range(self.distr_magnitude):
                        index = (x[dim], self.n[dim], self.p[dim][state][comp])
                        if lookup.has_key( index ):
                            res += lookup[index] * self.c[dim][state][comp]
                            k += 1
                        else:
                            y = binom.logpmf(x[dim], self.n[dim], self.p[dim][state][comp])
                            lookup[index] = y
                            res += y * self.c[dim][state][comp]
                row.append(res)
                
            matrix.append(row)
#         print('time to compute log-likelihood matrix (compute_log_likelihood): ',\
#               time()-t, file=sys.stderr)
        #print("HALLO")
        #print(np.asarray(matrix))
        return np.asarray(matrix)

    def _generate_sample_from_state(self, state, random_state=None):
        res = []
        for dim in range(self.n_features):
            erg = round(sum([binom.rvs(self.n[dim], self.p[dim][state][comp]) * self.c[dim][state][comp] for comp in range(self.distr_magnitude)]))
            res.append(erg)
        
        return np.array(res)
    
    def _initialize_sufficient_statistics(self):
        stats = super(BinomialHMM2d3s, self)._initialize_sufficient_statistics()
        stats['post'] = np.zeros(self.n_components)
        stats['post_emission'] = np.zeros([self.n_features, self.n_components])
        stats['post_sum_l'] = np.zeros([self.n_features, self.n_components, self.distr_magnitude])
        stats['post_sum_l_emisson'] = np.zeros([self.n_features, self.n_components, self.distr_magnitude])
        stats['weights'] = np.empty([self.n_features, self.n_components, self.distr_magnitude])
        stats['weights'].fill(0.5)
        return stats
    
    def _get_value(self, state, symbol, dim, stats):
        erg = 0
        for comp in range(self.distr_magnitude):
            erg += self.c[dim][state][comp] * binom.pmf(symbol[dim], self.n[dim], self.p[dim][state][comp])
        return erg
    
    def _help_accumulate_sufficient_statistics(self, obs, stats, posteriors):
        posteriors = _valid_posteriors(posteriors, obs)
        for t, symbol in enumerate(obs):
            stats['post'] += posteriors[t]
            
            for dim in range(self.n_features):
                stats['post_emission'][dim] += (posteriors[t] * symbol[dim])
            
            for dim in range(self.n_features):
                for state in range(self.n_components):
                    for comp in range(self.distr_magnitude):
                        enum = self.c[dim][state][comp] * binom.pmf(symbol[dim], self.n[dim], self.p[dim][state][comp])
                        denum = self._get_value(state, symbol, dim, stats)
                        stats['post_sum_l'][dim][state][comp] += posteriors[t][state] * enum / denum
                        stats['post_sum_l_emisson'][dim][state][comp] += posteriors[t][state] * enum / denum * symbol[dim]
            
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
        for dim in range(self.n_features):
            for state in range(self.n_components):
                for comp in range(self.distr_magnitude):
                    self.p[dim][state][comp] = stats['post_sum_l_emisson'][dim][state][comp] / (self.n[dim] * _add_pseudo_counts(stats['post_sum_l'][dim][state][comp])) 
                    self.p[dim][state][comp] = _add_pseudo_counts(self.p[dim][state][comp])
                    self.c[dim][state][comp] = stats['post_sum_l'][dim][state][comp] / stats['post'][state]

        #self.merge_distr()
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
        #tmp=np.array(map(lambda x: x*self.n[0], self.p))
        #print(np.reshape(tmp, (-1,3)), file=sys.stderr)

if __name__ == '__main__':
#    transmat_ = np.array(([[0.7, 0.2, 0.1], [0.1, 0.8, 0.1], [0.2,0.2,0.6]]))
    #tmp = [[0.000001, 0.0000098, 0.000001], [0.000001, 0.000001, 0.0000098]]
    tmp = [[[0.3, 0.3, 0.2], [0.6, 0.8, 0.7], [0.7, 0.8, 0.7]], [[0.2, 0.2, 0.1], [0.4, 0.8, 0.8], [0.7, 0.8, 0.9]]]
    #[[0.01, 0.98, 0.01], [0.01, 0.01, 0.98]]
    p_ = np.array(tmp)
    n_ = np.array([2000000, 2000000])
    
    m = BinomialHMM2d3s(startprob=[1,0,0], n = [100, 100])

    X, Z = m.sample(10) #returns (obs, hidden_states)
    #print(X,Z)
    
    m2 = BinomialHMM2d3s(n=[100, 100], p=tmp, startprob=[1,0,0])
#     cProfile.run("m2.fit([X])")
    m2.fit([X])
    print('obs', 'states', 'estimated states', sep='\t')
    e = m2.predict(X)
    for i, el in enumerate(X):
        print(el, Z[i], e[i], Z[i] == e[i], sep='\t')
#     print(m.p)
#     print(m2.p)

    