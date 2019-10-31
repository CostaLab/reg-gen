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


from mpmath import gamma, rf, loggamma
import numpy as np
from math import log
from numpy.random import random_sample
import sys
from math import fabs
from scipy.stats import nbinom

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
            #print('_get_value_log ValueError', x, mu, v, file=sys.stderr)
            return 1
    
    def _get_value(self, x, mu, v):
        try:
            return rf(v, x) / gamma(x+1) * ( v/float(v+mu) ) ** v * ( mu/float(v+mu) ) ** x
        except ValueError:
            print(x, mu, v, file=sys.stderr)
            return 1
            
    def __init__(self, mu, alpha):
        self.map_pdf = {}
        self.map_logpdf = {}
        self.bins = []
        mu = float(mu)
        
        self.alpha = alpha
        self.mu = mu
        
        self.nbin = np.frompyfunc(self._get_value, 3, 1)
        self.nbin_log = np.frompyfunc(self._get_value_log, 3, 1)
        
    def pdf(self, x):
        if x not in self.map_pdf:
            v_log = self.nbin(x, self.mu, 1./self.alpha)
            self.map_pdf[x] = v_log
        
        return self.map_pdf[x]
    
    def logpdf(self, x):
        if x not in self.map_logpdf:
            v_log = self.nbin_log(x, self.mu, 1./self.alpha)
            self.map_logpdf[x] = v_log
        
        return self.map_logpdf[x]
        
    def rvs(self):
        if not self.bins:
            probs = []
            i = 0
            v_old = -1
            while True:
                v = self.nbin(i, self.mu, 1./self.alpha)
                probs.append(v)
                
                i += 1
                if fabs(v_old - v) < 10**-10:
                    break
                v_old = v
            self.bins = [float(x) for x in np.add.accumulate(probs)]
        return np.digitize(random_sample(1), self.bins)[0]

if __name__ == '__main__':
    neg_bin = NegBin(0.1, 0.00000000001)
    s=0
    ew = 0
    distr = {'n': 10, 'p': 0.1}
    
    for i in range(10):
        #v = neg_bin.pdf(i)
        v = nbinom.logpmf(i, distr['n'], distr['p'])
        #v_log = neg_bin.logpdf(i)
        #s += v
        #ew += i * v
        #print(i, v, log(v), v_log, sep='\t')
        print(i, v, sep='\t')
    #print(s, ew)
