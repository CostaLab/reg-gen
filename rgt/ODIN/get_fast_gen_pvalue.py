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

def get_value(x, distr):
    if distr['distr_name'] == 'binomial':
        if lookup_pmf.has_key(x):
            return lookup_pmf[x]
        else:
            v = binom.pmf(x, distr['n'], distr['p'])
            lookup_pmf[x] = v
            return v

def get_log_value(x, distr):
    if distr['distr_name'] == 'binomial':
        if lookup_pmf.has_key(x):
            return lookup_pmf[x]
        else:
            v = binom.logpmf(x, distr['n'], distr['p'])
            lookup_pmf[x] = v
            return v
    if distr['distr_name'] == 'nb':
        return distr['distr'].logpdf(x)
        #return nbinom.logpmf(x, distr['n'], distr['p'])

def _comp(i, x, side, current_p, p):
    if side == 'r':
        return i <= x
    elif side == 'l':
        return i >= x
    else:
        return current_p >= p


def compute_pvalue(distr, N, side, current_p, x):
    """Compute log2 pvalue"""
    sum_num = []
    sum_denum = []
    it = range(N/2+1) if side == 'r' else range(N+1, -1, -1)
    
    for i in it:
        p1 = get_log_value(i, distr)
        p2 = get_log_value(N - i, distr)
        p = p1 + p2
        
        if _comp(i, x, side, current_p, p):
        #if p > current_p:
            sum_num.append(p)
        sum_denum.append(p)
        
    if distr['distr_name'] == 'nb':
        sum_num = map(lambda x: float(x), sum_num)
        sum_denum = map(lambda x: float(x), sum_denum)
    
    return logsumexp(np.array(sum_num)) - (log(2) + logsumexp(np.array(sum_denum)))
    
def get_log_pvalue_new(x, y, side, distr):
    """compute log10 p-value"""
    N = x + y
    if side == 'l':
        x, y = y, x
        side = 'r'
    
    if lookup_pvalue.has_key((x, y, 'r')):
        return lookup_pvalue[(x, y, 'r')]
    else:
        current_p = get_log_value(x, distr) + get_log_value(y, distr)
        pvalue = compute_pvalue(distr, N, side, current_p, x) / log(10)
        lookup_pvalue[(x, y, side)] = pvalue
        return pvalue

def change_nb_WP2NB1(n, p):
    alpha = 1./n
    mu = (1./p - 1) / alpha
    
    return mu, alpha

def change_nb_NB12WP(mu, alpha):
    alpha = float(alpha)
    mu = float(mu)
    p = 1./(1 + mu * alpha)
    n = 1./alpha
    
    return n, p

if __name__ == '__main__':
    mu = 1.03882161264 
    alpha = 0.1
    
    m = NegBin(mu, alpha)
    distr = {'distr_name': 'nb', 'distr': m}
    
    #,0.510793370086
    for x, y in [([800, 900],[600, 500]), ([200, 190], [40,50])]:
        side = 'l' if x > y else 'r'
        var =  np.var( x + y )
        mu = np.mean( x + y )
        alpha = max((var - mu) / np.square(mu), 0.00000000001)
        m = NegBin(mu, alpha)
        distr = {'distr_name': 'nb', 'distr': m}
        print(x, y, -get_log_pvalue_new(int(sum(x)), int(sum(y)), side, distr), sep='\t')
    
    
    #n = 90
    #p = 0.01
    #distr={'distr_name': 'binomial', 'p':p, 'n':n}
    
    #n, p = change_nb_NB12WP(mu, alpha)
    
    #n, p = 10, 0.1 #working
    #print(n, p)
    #distr={'distr_name': 'nb', 'p': p, 'n': n}
    
    #for i in range(10):
        #print(nbinom.logpmf(i, distr['n'], distr['p']), m.logpdf(i), fabs(nbinom.logpmf(i, distr['n'], distr['p']) - round(m.logpdf(i), 11)) < 10**-10, sep='\t')
    
    #,,0.335898768556
    x, y, side = 800, 600, 'l'
    #print(x, y, -get_log_pvalue_new(x, y, side, distr), sep='\t')
    
    
    #for x,y in [(800, 600), (12, 5)]:
        #side = 'l' if x > y else 'r'
    #x, y, side = 12, 5, 'l'
        #print(x, y, -get_log_pvalue_new(x, y, side, distr), sep='\t')
    
    #print()
    
    S = 30
    for x in range(S+1):
        y = S - x
        side = 'l' if x > y else 'r'
        cur = 10**(get_log_value(x, distr))
        cur2 = 10**(get_log_value(y, distr))
        #print(x, y, cur, cur2, 10**get_log_pvalue_new(x, y, side, distr), sep='\t')
        #print(cur, cur2, sep='\t')
    