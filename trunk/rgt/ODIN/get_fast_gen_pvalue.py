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
from math import log, fabs, exp
from sklearn.utils.extmath import logsumexp
import numpy as np
import sys

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

def _comp(i, x, side, current_p, p):
    if side == 'r':
        return i <= x
    elif side == 'l':
        return i >= x
    else:
        return current_p >= p


def _check_side(side, i, S):
    if side == 'r' and i > S/2:
        return True
    elif side == 'l' and i <= S/2:
        return True
    else:
        return False
    
def compute_pvalue(distr, N, side, current_p, x):
    sum_num = []
    sum_denum = []
    it = range(N/2+1) if side == 'r' else range(N+1, -1, -1)
    
    for i in it:
        p1 = get_log_value(i, distr)
        p2 = get_log_value(N - i, distr)
        p = p1 + p2
        
        if _comp(i, x, side, current_p, p):
            sum_num.append(p)
        sum_denum.append(p)
        
        #with t is precision. does not help so much, so ignore it
        #if side == 'r' and i > S/2 and p < t:
        #    print(i ,len(it), p, t)
        #    sum_denum += (S - i) * [p]
        #    if _comp(i, x, side, current_p, p):
        #        sum_num += (S - i) * [p]
        #    break
        
        #if side == 'l' and i < S/2 and p < t:
        #    sum_denum += i * p
        #    if _comp(i, x, side, current_p, p):
        #        sum_num += i * p            
        #    break
    #sum_num += sum_num
    #sum_denum += sum_denum
    return logsumexp(np.array(sum_num)) - (log(2) + logsumexp(np.array(sum_denum)))

def get_log_pvalue_new(x, y, side, distr):
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

if __name__ == '__main__':
    n = 30000000
    p = 0.01
    distr={'distr_name': 'binomial', 'p':p, 'n':n}
    #side = 'l' #r bei y > x, sonst l
    #x=870
    #y=1123
    #side = 'r' if y > x else 'l'
    #rep = 1
    
    values = []
    with open('/home/manuel/pvalues_cut', 'r') as f:
        for line in f:
            line=line.strip()
            line=line.split(',')
            values.append((int(line[0]), int(line[1]), float(line[2])))

    i = 0
    t = time()
    k = 0
    for x, y, pvalue in values:
        k += 1
        if k == 100:
            break
        side = 'r' if y > x else 'l'
        p = -get_log_pvalue_new(x, y, side, distr)
        print(p)
        if fabs(pvalue - p) > 1:
            i += 1
            #print('r',x, y, pvalue, p, sep= '\t')
    t2 = time()-t
    print('new approach:', t2, i)
     
#    t = time()
#    i = 0
#    for x, y, pvalue in values:
#        i += 1
#        if i % 10 == 0:
#            print(i, file=sys.stderr)
#        
#        side = 'r' if y > x else 'l'
#        p = get_log_pvalue(x, y, side, distr_name='binomial', p=p, n=n)
#        if fabs(pvalue - p) > 10:
#            print(x, y, pvalue, p)
#    t2_o = time()-t
#    print('old approach:', t2_o)
#    print('The new method is %s time(s) faster... awesome!' %((t2_o/t2)))
    
    
#     x,y = 10,12
#     side = 'r' if y > x else 'l'
#     print(get_pvalue_new(x,y,distr,side))
#     print(get_log_pvalue(x, y, side, distr_name='binomial', p=p, n=n))
    
    
    #new
    #t = time()
    #for i in range(rep):
    #    pvalue = get_pvalue_new(x, y, distr, side)
    #told = time()-t
    #print('time for new: ', told, pvalue)
    
    #old
    #t2 = time()
    #for i in range(rep):
    #    pvalue = get_log_pvalue(x, y, side, distr_name='binomial', p=p, n=n)
    #t2old = time()-t2
    #print('time for old: ', t2old, pvalue)
    #print('The new method is %s time(s) faster... awesome!' %((t2old/told)))
    
    #print(get_log_pvalue(12, 3, 'r', distr_name='binomial', p=p, n=n))
    #print(get_log_pvalue(3, 12, 'l', distr_name='binomial', p=p, n=n))
    
    
    #S = 4
    #for x in range(S+1):
    #    print('new', x, S-x, get_pvalue_new(x, S-x, distr, side))
    #    #a = get_pvalue(x, S-x, side, distr_name='binomial', p=p, n=n)
    #    #print('old', x, S-x, a, log(a, 10))
    #    print('old', x, S-x, get_log_pvalue(x, S-x, side, distr_name='binomial', p=p, n=n))
    #    #print(get_pvalue_new(x, S-x, distr, side) - get_pvalue(x, S-x, side, distr_name='binomial', p=p, n=n) < 10**-10)
