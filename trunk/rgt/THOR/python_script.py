#!/usr/bin/env python
#  -*- coding: utf-8 -*-
'''
Created on Mar 6, 2013

@author: manuel

Consider the 2x1 contingency table with x and y 
and fixed marginal column sum n=x+y. Value x follows h1, 
value y follows h2.
    
'''

from __future__ import print_function
from time import time
from scipy.stats import nbinom, binom
from math import log, fabs, exp, log10
from sklearn.utils.extmath import logsumexp
import numpy as np
import sys
from rgt.THOR.neg_bin import NegBin

lookup_pmf = {}
lookup_pvalue = {}

def get_log_value(x, distr):
    if not lookup_pmf.has_key(x):
        v = float(distr['distr'].logpdf(x))
        lookup_pmf[x] = v
    
    return lookup_pmf[x]

def compute_pvalue(distr, N, x, current_p):
    """Compute log2 pvalue"""
    sum_num = []
    sum_denum = []
    
    for i in range(N+1):
        p1 = get_log_value(i, distr)
        p2 = get_log_value(N - i, distr)
        p = p1 + p2
        
        if current_p >= p:
            sum_num.append(p)
        
        sum_denum.append(p)

    return logsumexp(np.array(sum_num)) - logsumexp(np.array(sum_denum))
    
def get_log_pvalue(x, y, side, distr):
    """compute log10 p-value"""
    N = x + y
    
    if side == 'l':
        x, y = y, x
        side = 'r'
    
    if not lookup_pvalue.has_key((x, y)):
        current_p = get_log_value(x, distr) + get_log_value(y, distr)
        pvalue = compute_pvalue(distr, N, x, current_p) / log(10)
        lookup_pvalue[(x, y)] = pvalue
    return lookup_pvalue[(x, y)]

def get_d(lower, upper, s_l, s_r, m):
    return sum([m.pdf(s_l + i) * m.pdf(s_r - i) for i in range(lower, upper + 1)])

if __name__ == '__main__':
    mu, alpha = 10, 0.1
    m = NegBin(mu, alpha)
    distr = {'distr_name': 'nb', 'distr': m}
     
    x, y, side = 800, 600, 'l'
#     x, y, side = 30, 20, 'l'
    
    print(x, y, -get_log_pvalue(x, y, side, distr), sep='\t')
     
    x, y, side = 30, 20, 'l'
    print(x, y, -get_log_pvalue(x, y, side, distr), sep='\t')
    
    #x, y, side = 30, 10, 'l'
    #print(x, y, -get_log_pvalue(x, y, side, distr), sep='\t')
    
    x, y, side = 40, 10, 'l'
    print(x, y, -get_log_pvalue(x, y, side, distr), sep='\t')
    
    print()
    S = 10
    for x in range(S+1):
        y = S - x
        side = 'l' if x > y else 'r'
               
        pvalue = 10**get_log_pvalue(x, y, side, distr)
     
        #print(x, y, pvalue, sep='\t')

    y_1 = 800
    y_2 = 600
      
    y__1 = 40
    y__2 = 10
      
    a = get_d(1, y_2, y_1, y_2, m)
    b = get_d(0, y_1 + y_2, 0, y_1+y_2, m)
      
    c = get_d(1, y__2, y__1, y__2, m)
    d = get_d(0, y__1 + y__2, 0, y__1 + y__2, m)
      
    print(a, b)
    print(c, d)
      
    print(-log10(a/b), -log10(c/d))



