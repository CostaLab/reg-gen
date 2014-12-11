#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
%prog

Implementation ofNegative Binomial distribution.

@author: Manuel Allhoff

"""
from __future__ import print_function
from mpmath import gamma, rf, loggamma
import numpy as np
from math import log
from numpy.random import random_sample
import sys

class NegBin():
    """Negative Binomial distribution (NB1) with continuous parameter r,
    for NB1 see 
    Noriszura Ismail and Abdul Aziz Jemain. 
    Handling Overdispersion with Negative Binomial and Generalized Poisson Regression Models. 
    
    Implementation based on
    
    - http://stackoverflow.com/questions/11373192/generating-discrete-random-variables-with-specified-weights-using-scipy-or-numpy
    - http://www.nehalemlabs.net/prototype/blog/2013/11/11/negative-binomial-with-continuous-parameters-in-python/
    """
    
    def _get_value_log(self, x, mu, v):
        """log basic 2"""
        try:
            return loggamma(x+v) - loggamma(x+1) - loggamma(v) + v*log(v) - v*log(v+mu) + x*log(mu) - x*log(v+mu)
        except ValueError:
            print('_get_value_log ValueError', x, mu, v, file=sys.stderr)
            return 1
    
    def _get_value(self, x, mu, v):
        try:
            return rf(v, x) / gamma(x+1) * ( v/float(v+mu) ) ** v * ( mu/float(v+mu) ) ** x
        except ValueError:
            print(x, mu, v, file=sys.stderr)
            return 1
            
    def __init__(self, mu, alpha, max_range=200000):
        self.map_pdf = {}
        self.map_logpdf = {}
        self.bins = []
        mu = float(mu)
        self.max_range = max_range
        
        self.alpha = alpha
        self.mu = mu
        
        self.nbin = np.frompyfunc(self._get_value, 3, 1)
        
        self.nbin_log = np.frompyfunc(self._get_value_log, 3, 1)
        
        probs = []
        probs_log = []
        for i in range(self.max_range):
            #print(i, self.mu, 1./self.alpha, file=sys.stderr)
            v = self.nbin(i, self.mu, 1./self.alpha)
            #print(i, self.mu, self.alpha)
            v_log = self.nbin_log(i, self.mu, 1./self.alpha)
            
            probs.append(v)
            probs_log.append(v_log)
            
            self.map_pdf[i] = v
            self.map_logpdf[i] = v_log
        
        self.bins = map(lambda x: float(x), np.add.accumulate(probs))
        
        #c = 5000
        #self.dist = rv_discrete(values=([i for i in range(c)], map(lambda x: float(x), [self.pdf(i) for i in range(c)])), name='dist')

    def pdf(self, x):
        #global map_pdf
        if x >= self.max_range:
            print('overflow pdf ', x, self.max_range, file=sys.stderr)
            return self.map_pdf[self.max_range-1]
        else:
            return self.map_pdf[x]
    
    def logpdf(self, x):
        #global map_logpdf
        if x >= self.max_range:
            print('overflow logpdf ', x, self.max_range, file=sys.stderr)
            return self.map_logpdf[self.max_range-1]
        else:
            return self.map_logpdf[x]
    
    def rvs(self):
        return np.digitize(random_sample(1), self.bins)[0]
#     def rvs(self, size=1):
#         return self.dist.rvs(size=size)
    

if __name__ == '__main__':
    neg_bin = NegBin(15, 0.2, max_range=50)
    s=0
    ew = 0
    for i in range(100):
        v= neg_bin.pdf(i)
        v_log = neg_bin.logpdf(i)
        s += v
        ew += i * v
        print(i, v, log(v), v_log, sep='\t')
    print(s, ew)

#     s = 0
#     for i in range(1000):
#         s += neg_bin.rvs()
#     print(s)
        
    
    
    
    
    
    
    
    
    
