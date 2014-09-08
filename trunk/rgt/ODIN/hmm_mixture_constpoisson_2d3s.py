#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on Mar 6, 2013

@author: manuel
'''
from __future__ import print_function
from scipy.stats import poisson
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
    def __init__(self, distr_magnitude, factors, init_state_seq=None, p = [[[3, 2, 1], [12, 15, 20], [2, 1, 1]], [[3, 2, 3], [4, 2, 1], [15, 16, 18]]], \
                 c = [[[0.5, 0.4, 0.1], [0.5, 0.4, 0.1], [0.5, 0.4, 0.1]], [[0.5, 0.4, 0.1], [0.5, 0.4, 0.1], [0.5, 0.4, 0.1]]], n_components=3, covariance_type='diag', startprob=[1,0,0],
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
        self.c = c # 1) dim 2) state 3) component
        self.p = p # 1) dim 2) state 3) component, parameter of Poisson distribution
        self.n_features = 2 #emission dimension
        self.init_state_seq = init_state_seq
        self.distr_magnitude = distr_magnitude
        self.factors = factors #weight of the posteriors
        assert len(self.factors) == self.distr_magnitude
        self.weights = 0 #weiths for merging distributions
        
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
                        index = (x[dim], self.p[dim][state][comp])
                        if lookup.has_key( index ):
                            res += lookup[index] * self.c[dim][state][comp]
                            k += 1
                        else:
                            y = poisson.logpmf(x[dim], self.p[dim][state][comp])
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
            erg = round(sum([poisson.rvs(self.p[dim][state][comp]) * self.c[dim][state][comp] for comp in range(self.distr_magnitude)]))
            res.append(erg)
        
        return np.array(res)
    
    def _initialize_sufficient_statistics(self):
        stats = super(BinomialHMM2d3s, self)._initialize_sufficient_statistics()
        stats['post'] = np.zeros(self.n_components)
#         stats['post_emission'] = np.zeros([self.n_features, self.n_components])
        stats['post_sum_l'] = np.zeros([self.n_features, self.n_components, self.distr_magnitude])
        stats['post_sum_l_emisson'] = np.zeros([self.n_features, self.n_components, self.distr_magnitude])
        stats['post_sum_l_factor'] = np.zeros([self.n_features, self.n_components, self.distr_magnitude])
        stats['post_l'] = [[[0]*self.distr_magnitude for _ in range(self.n_components)] for _ in range(self.n_features)]
        
        return stats
    
    def _get_value(self, state, symbol, dim, stats):
        erg = 0
        for comp in range(self.distr_magnitude):
            erg += self.c[dim][state][comp] * poisson.pmf(symbol[dim], self.p[dim][state][comp])
        return erg
    
    def _help_accumulate_sufficient_statistics(self, obs, stats, posteriors):
        posteriors = _valid_posteriors(posteriors, obs)
        
        for dim in range(self.n_features):
            for state in range(self.n_components):
                for comp in range(self.distr_magnitude):
                    stats['post_l'][dim][state][comp] = [0] * len(obs)

        for t, symbol in enumerate(obs):
            stats['post'] += posteriors[t]
            for dim in range(self.n_features):
                for state in range(self.n_components):
                    for comp in range(self.distr_magnitude):
                        enum = self.c[dim][state][comp] * poisson.pmf(symbol[dim], self.p[dim][state][comp])
                        denum = self._get_value(state, symbol, dim, stats)
                        help = posteriors[t][state] * enum / denum
                        stats['post_sum_l'][dim][state][comp] += help
                        stats['post_sum_l_emisson'][dim][state][comp] += help * symbol[dim]
                        stats['post_sum_l_factor'][dim][state][comp] += help * self.factors[comp]
                        stats['post_l'][dim][state][comp][t] = help
        
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
                    self.c[dim][state][comp] = stats['post_sum_l'][dim][state][comp] / stats['post'][state]
                    if comp == 0:
                        self.p[dim][state][comp] = stats['post_sum_l_emisson'][dim][state][comp] / (_add_pseudo_counts(stats['post_sum_l_factor'][dim][state][comp])) 
                        self.p[dim][state][comp] = _add_pseudo_counts(self.p[dim][state][comp])
                    else:
                        self.p[dim][state][comp] = self.factors[comp] * self.p[dim][state][0]
        self.merge_distr()
        #print('m-step: ',self.p, file=sys.stderr)
        #tmp=np.array(map(lambda x: x*self.n[0], self.p))
        #print(np.reshape(tmp, (-1,3)), file=sys.stderr)
        
    
    def _count(self, posts_l):
        res = [[[[1, 1] for _ in range(self.distr_magnitude)] for _ in range(self.n_components)] for _ in range(self.n_features)]
        for dim in range(self.n_features):
            for comp in range(self.distr_magnitude):
                for i in range(len(posts_l[dim][1][comp])):
                    if posts_l[dim][1][comp][i] >= posts_l[dim][2][comp][i]:
                        res[dim][1][comp][0] += 1
                    if posts_l[dim][2][comp][i] > posts_l[dim][1][comp][i]:
                        res[dim][2][comp][1] += 1

        return res
        
    
    def _do_mstep(self, stats, params):
        super(BinomialHMM2d3s, self)._do_mstep(stats, params)
        self.weights = self._count(stats['post_l'])
        self._help_do_mstep(stats)
       
    def merge_distr(self):
        for dim in range(self.n_features):
            for comp in range(self.distr_magnitude):
                for state1, state2 in [[1, 2], [2, 1], [0, 0]]:
                    dim2 = 1 if dim == 0 else 0
                    c1, c2 = self.weights[dim][state1][comp][0], self.weights[dim][state1][comp][1]
                    f = c2 / float(c1 + c2)
                    p_norm = self.p[dim][state1][comp] + f * fabs(self.p[dim][state1][comp] - self.p[dim2][state2][comp])
                    #print('merge: ', f, self.p, p_norm, file=sys.stderr)
                    
                    self.p[dim][state1][comp] = p_norm
                    self.p[dim2][state2][comp] = p_norm
                
                #print('merge: ', f, self.p, p_high, p_low, file=sys.stderr)
                #tmp=np.array(map(lambda x: x*self.n[0], self.p))
                #print(np.reshape(tmp, (-1,3)), file=sys.stderr)

if __name__ == '__main__':
    #3 components
    distr_magnitude = 3
    factors = [1,2,3]
     
    tmp1 = [[[3, 2, 1], [12, 15, 20], [2, 1, 1]], [[3, 2, 3], [4, 2, 1], [15, 16, 18]]]
    c1 = [[[0.5, 0.4, 0.1], [0.5, 0.4, 0.1], [0.5, 0.4, 0.1]], [[0.5, 0.4, 0.1], [0.5, 0.4, 0.1], [0.5, 0.4, 0.1]]]
     
    tmp2 = [[[4, 4, 1], [15, 16, 20], [2, 2, 3]], [[3, 1, 1], [1, 2, 1], [12, 13, 14]]]
    c2 = [[[0.5, 0.4, 0.1], [0.5, 0.4, 0.1], [0.5, 0.4, 0.1]], [[0.5, 0.4, 0.1], [0.5, 0.4, 0.1], [0.5, 0.4, 0.1]]]
    
#     #2 components
#     distr_magnitude = 2
#     factors = [1,2]
#     
#     tmp1 = [[[3, 2], [12, 15], [2, 1]], [[3, 2], [4, 2], [15, 16]]]
#     c1 = [[[0.5, 0.4], [0.5, 0.4], [0.5, 0.4]], [[0.5, 0.4], [0.5, 0.4], [0.5, 0.4]]]
#     
#     tmp2 = [[[4, 4], [15, 16], [2, 2]], [[3, 1], [1, 2], [12, 13]]]
#     c2 = [[[0.5, 0.4], [0.5, 0.4], [0.5, 0.4]], [[0.5, 0.4], [0.5, 0.4], [0.5, 0.4]]]
    
    
    
    m = BinomialHMM2d3s(p=tmp1, c=c1, distr_magnitude=distr_magnitude, factors = factors)

    X, Z = m.sample(100) #returns (obs, hidden_states)
    
    m2 = BinomialHMM2d3s(p=tmp2, c=c2, distr_magnitude=distr_magnitude, factors = factors)
    m2.fit([X])
      
    print('obs', 'states', 'estimated states', sep='\t')
    e = m2.predict(X)
    for i, el in enumerate(X):
        print(el, Z[i], e[i], Z[i] == e[i], sep='\t')
    print(m.p)
    print(m2.p)
    print(m2.c)
    