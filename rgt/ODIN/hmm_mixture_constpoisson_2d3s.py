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
import cProfile
import trace
from memory_profiler import profile
from array import array

lookup_poisson = {}
lookup_state = {}
lookup_poisson_state = {}
lookup_denum = {}

def get_init_parameters(s1, s2, **info):
    #get observation that occurs most often:
    distr_magnitude = int(info['distr_magnitude'])
    n_components = int(info['n_components'])
    n_features = int(info['n_features'])
    
    #emp_mean =[float(np.argmax(np.bincount(map(lambda x: x[0], s1)))), float(np.argmax(np.bincount(map(lambda x: x[1], s2)))) ]
    
    emp_mean = [np.mean(map(lambda x: x[0], s1)), np.mean(map(lambda x: x[1], s2))]
                                                 
    initial_c = [[[1/float(info['distr_magnitude']) for _ in range(n_components)] for _ in range(distr_magnitude)] for _ in range(n_features)]
    initial_p = [[[0 for _ in range(n_components)] for _ in range(distr_magnitude)] for _ in range(n_features)]
    
    for dim in range(n_features):
        for comp in range(distr_magnitude):
            for state in range(n_components):
                if state == 0:
                    background_value = max(1, emp_mean[dim] / 100.)
                    if comp == 0:
                        initial_p[dim][comp][state] = background_value
                    elif comp > 0:
                        initial_p[dim][comp][state] = (comp+1) * initial_p[dim][0][state]
                
                indices = [1,2] if dim == 0 else [2,1]
                if state > 0:
                    if comp == 0:
                        initial_p[dim][comp][indices[0]] = distr_magnitude * emp_mean[dim] / float(sum(range(1, distr_magnitude+1)))
                        initial_p[dim][comp][indices[1]] = initial_p[dim][comp][0]
                    elif comp > 0:
                        initial_p[dim][comp][state] = (comp+1) * initial_p[dim][0][state]
    
    return np.array(initial_c, np.float64), np.array(initial_p, np.float64)

class PoissonHMM2d3s(_BaseHMM):
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
        self.c = c # 1) dim 2) component 3) state
        self.p = p # 1) dim 2) component, parameter of Poisson distribution 3) state 
        self.n_features = 2 #emission dimension
        self.init_state_seq = init_state_seq
        self.distr_magnitude = distr_magnitude
        self.factors = factors #weight of the posteriors
        assert len(self.factors) == self.distr_magnitude
        self.weights = 0 #weigths for merging distributions
        
    def _compute_log_likelihood(self, X):
        matrix = []
        lookup = {}
        k = 0
        for x in X:
            row = []
            for state in range(self.n_components): #state
                res = 0
                for dim in range(self.n_features): #dim
                    for comp in range(self.distr_magnitude):
                        index = (x[dim], self.p[dim][comp][state])
                        if lookup.has_key( index ):
                            res += lookup[index] * self.c[dim][comp][state]
                            k += 1
                        else:
                            y = poisson.logpmf(x[dim], self.p[dim][comp][state])
                            lookup[index] = y
                            res += y * self.c[dim][comp][state]
                row.append(res)
                
            matrix.append(row)
        return np.asarray(matrix)

    def _generate_sample_from_state(self, state, random_state=None):
        res = []
        for dim in range(self.n_features):
            erg = round(sum([poisson.rvs(self.p[dim][comp][state]) * self.c[dim][comp][state] for comp in range(self.distr_magnitude)]))
            res.append(erg)
        
        return np.array(res)
    
    def _initialize_sufficient_statistics(self):
        stats = super(PoissonHMM2d3s, self)._initialize_sufficient_statistics()
        stats['post'] = np.zeros(self.n_components)
        stats['post_sum_l'] = np.zeros([self.n_features, self.distr_magnitude, self.n_components])
        stats['post_sum_l_emisson'] = np.zeros([self.n_features, self.distr_magnitude, self.n_components])
        stats['post_sum_l_factor'] = np.zeros([self.n_features, self.distr_magnitude, self.n_components])
        stats['weights'] = [[[[1, 1] for _ in range(self.n_components)] for _ in range(self.distr_magnitude)] for _ in range(self.n_features)]
        
        return stats
    
    def _get_poisson(self, x, p):
        if lookup_poisson.has_key((x, p)):
                value_poisson = lookup_poisson[(x, p)]
        else:
            value_poisson = poisson.pmf(x, p)
            lookup_poisson[(x, p)] = value_poisson
        return value_poisson    
    
    def _get_value(self, state, symbol, dim):
        help_i = [self.c[dim][i][state] for i in range(self.distr_magnitude)]
        index = (state, symbol[dim], tuple(help_i))
        
        if index not in lookup_state:
            res = 0
            for comp in range(self.distr_magnitude):
                res += self.c[dim][comp][state] * self._get_poisson(symbol[dim], self.p[dim][comp][state])
            lookup_state[index] = res
        
        return lookup_state[index]
    
    
    def _help_accumulate_sufficient_statistics(self, obs, stats, posteriors):
        posteriors = _valid_posteriors(posteriors, obs)
        i = 0
        print("run...!", file=sys.stderr)
        for t, symbol in enumerate(obs):
            stats['post'] += posteriors[t]
            for dim in range(self.n_features):
                for comp in range(self.distr_magnitude):
                    #lookup
                    index = (symbol[dim], tuple([self.p[dim][comp][state] for state in range(self.n_components)]))
                    if index not in lookup_poisson_state: 
                        tmp = np.array([self._get_poisson(symbol[dim], self.p[dim][comp][state]) for state in range(self.n_components)])
                        lookup_poisson_state[index] = tmp
                    h = lookup_poisson_state[index]
                    
                    enum = self.c[dim][comp] * h
                    denum = np.array([self._get_value(state, symbol, dim) for state in range(self.n_components)])
                    
                    i += 1
                    try:
                        help = (posteriors[t] * enum / _add_pseudo_counts(denum))
                    except:
                        print("%s \n" %i, file=sys.stderr)
                        print("%s %s %s \n" %(denum, symbol, dim), file=sys.stderr)
                        print("%s \n" %(self.c), file=sys.stderr)
                        print("%s \n" %(self.p), file=sys.stderr)
                        print("%s \n" %(posteriors[t]), file=sys.stderr)
                        print("%s \n" %(enum), file=sys.stderr)
                        help = np.array([1.0/self.distr_magnitude, 1.0/self.distr_magnitude, 1.0/self.distr_magnitude])
                    stats['post_sum_l'][dim][comp] += help
                    stats['post_sum_l_emisson'][dim][comp] += help * symbol[dim]
                    stats['post_sum_l_factor'][dim][comp] += help * self.factors[comp]
                    
                    
                    if help[1] >= help[2]:
                        stats['weights'][dim][comp][1][0] += 1
                    if help[2] > help[1]:
                        stats['weights'][dim][comp][2][1] += 1
                    
        stats['posterior'] = np.copy(posteriors)
    
    def _accumulate_sufficient_statistics(self, stats, obs, framelogprob,
                                      posteriors, fwdlattice, bwdlattice,
                                      params):
        super(PoissonHMM2d3s, self)._accumulate_sufficient_statistics(
            stats, obs, framelogprob, posteriors, fwdlattice, bwdlattice,
            params)
        
        self._help_accumulate_sufficient_statistics(obs, stats, posteriors)        
    
    def _help_do_mstep(self, stats):
        for dim in range(self.n_features):
            for comp in range(self.distr_magnitude):
                for state in range(self.n_components):
                    self.c[dim][comp][state] = stats['post_sum_l'][dim][comp][state] / stats['post'][state]
                    if comp == 0:
                        self.p[dim][comp][state] = stats['post_sum_l_emisson'][dim][comp][state] / (_add_pseudo_counts(stats['post_sum_l_factor'][dim][comp][state])) 
                        self.p[dim][comp][state] = _add_pseudo_counts(self.p[dim][comp][state])
                    else:
                        self.p[dim][comp][state] = self.factors[comp] * self.p[dim][0][state]
        self.merge_distr(stats['weights'])
    
    def _do_mstep(self, stats, params):
        super(PoissonHMM2d3s, self)._do_mstep(stats, params)
        self._help_do_mstep(stats)
    
    def get_mean(self, state, dim):
        erg = 0
        for comp in range(self.distr_magnitude):
            erg += self.c[dim][comp][state] * self.p[dim][comp][state]
        return erg
    
    def merge_distr(self, weights):
        for dim in range(self.n_features):
            for comp in range(self.distr_magnitude):
                for state1, state2 in [[1, 2], [2, 1], [0, 0]]:
                    dim2 = 1 if dim == 0 else 0
                    c1, c2 = weights[dim][comp][state1][0], weights[dim][comp][state1][1]
                    f = c2 / float(c1 + c2)
                    p_norm = self.p[dim][comp][state1] + f * fabs(self.p[dim][comp][state1] - self.p[dim2][comp][state2])
                    
                    self.p[dim][comp][state1] = p_norm
                    self.p[dim2][comp][state2] = p_norm
                
if __name__ == '__main__':
    #3 components
    distr_magnitude = 3
    factors = [1,2,3]
    tmp1 = np.array([[[3, 12, 2], [2, 15, 1], [1, 20, 1]], [[3, 4, 15], [2, 2, 16], [3, 1, 18]]], np.float64)
    c1 = np.array([[[0.2, 0.3, 0.4], [0.3, 0.4, 0.3], [0.5, 0.3, 0.3]], [[0.5, 0.4, 0.6], [0.4, 0.4, 0.3], [0.1, 0.2, 0.1]]], np.float64)
      
    tmp2 = np.array([[[2, 10, 4], [2, 11, 3], [3, 14, 1]], [[1, 4, 14], [3, 3, 15], [2, 12, 20]]], np.float64)
    c2 = np.array([[[0.1, 0.5, 0.3], [0.4, 0.3, 0.4], [0.5, 0.2, 0.3]], [[0.4, 0.3, 0.6], [0.4, 0.5, 0.3], [0.2, 0.2, 0.1]]], np.float64)
    
    #2 components
#     distr_magnitude = 2
#     factors = [1,2]
#       
#     tmp1 = np.array([[[3, 12, 2], [2, 15, 1]], [[3, 4, 15], [2, 2, 16]]], np.float64)
#     c1 = np.array([[[0.6, 0.7, 0.5], [0.4, 0.3, 0.5]], [[0.6, 0.8, 0.5], [0.4, 0.2, 0.5]]], np.float64)
#       
#     tmp2 = np.array([[[2, 9, 1], [3, 10, 2]], [[4, 2, 12], [2, 2, 12]]], np.float64)
#     c2 = np.array([[[0.5, 0.5, 0.5], [0.5, 0.5, 0.5]], [[0.5, 0.5, 0.5], [0.5, 0.5, 0.5]]], np.float64)
    
    #a, b = get_init_parameters([(10,1), (50,3), (50,2), (40,2)], [(0,10),(2,12),(10,16)], distr_magnitude=2, n_components=3, n_features=2)
    #print(b)
    
    
    m = PoissonHMM2d3s(p=tmp1, c=c1, distr_magnitude=distr_magnitude, factors = factors)

    X, Z = m.sample(3000) #returns (obs, hidden_states)
    m2 = PoissonHMM2d3s(p=tmp2, c=c2, distr_magnitude=distr_magnitude, factors = factors)
    
    #tracer = trace.Trace(ignoredirs = [sys.prefix, sys.exec_prefix], trace = 0)
    #tracer.run("m2.fit([X])")
    #r = tracer.results()
    #r.write_results(show_missing=True, coverdir="ergebnis") 

    cProfile.run("m2.fit([X])")
    #m2.fit([X])
      
    #print('obs', 'states', 'estimated states', sep='\t')
    #e = m2.predict(X)
    #for i, el in enumerate(X):
    #    print(el, Z[i], e[i], Z[i] == e[i], sep='\t')
    #print(m.p)
    #print(m2.p)
    #print(m2.c)
    