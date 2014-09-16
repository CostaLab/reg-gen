from __future__ import print_function
import sys
import numpy as np
from collections import Counter
from scipy import ndimage
import os.path
from math import log, fabs
from random import randint
import warnings

EPSILON=1e-320


    
def _valid_posteriors(posteriors, obs):
    
    warnings.filterwarnings('error')
    
    for i in range(len(obs)):
        state_1 = False
        c1, c2 = obs[i][0], obs[i][1] #counts of samples

        if posteriors[i][0] > 0.5: 
            state_1 = True
            
        if c1 > c2: #state 1
            if fabs(posteriors[i][2] - 1) < 1e-200:
                posteriors[i] = np.array([1, 0, 0])
            else:
                if not state_1 and posteriors[i][2] > posteriors[i][1]:
                    try:                        
                        post_s2 = 0
                        post_s0 = posteriors[i][0] / (posteriors[i][0] + posteriors[i][1])
                        post_s1 = posteriors[i][1] / (posteriors[i][0] + posteriors[i][1])
                        posteriors[i] = np.array([post_s0, post_s1, post_s2])
                    except RuntimeWarning:
                        print(posteriors[i], c1, c2, file=sys.stderr)
        
        if c2 > c1: #state 2
            if fabs(posteriors[i][1] - 1) < 1e-200:
                posteriors[i] = np.array([1, 0, 0])
            else:
                if not state_1 and posteriors[i][1] > posteriors[i][2]:
                    try:
                        post_s1 = 0
                        post_s0 = posteriors[i][0] / (posteriors[i][0] + posteriors[i][2])
                        post_s2 = posteriors[i][2] / (posteriors[i][0] + posteriors[i][2])
                        posteriors[i] = np.array([post_s0, post_s1, post_s2])
                    except RuntimeWarning:
                        print(posteriors[i], c1, c2, file=sys.stderr)
    
    return posteriors

def _init(self, obs, params):
    """init HMM for EM algorithm"""
    if params == 'advanced':
        if self.init_state_seq is None:
            raise ValueError("Initial state sequence must be defined when initialzing HMM.")
        state_seq = np.array(map(lambda x: x[0], self.init_state_seq)) 
        #init start probabilities
        startprob = [0.0] * self.n_components
        startprob[self.init_state_seq[0][0]] = 1
        self.startprob_ = startprob
        
        #init transmission_probabilities
        #look at state 0 and consider the following state, count these states with Counter 
        #and transform it to a dict, keys of dict give the following states (0,1,2), 
        #the corresponding values the number of times
        #Example: consider state 1, result; {0: 2, 1: 1, 2: 1} : state 0 follows state 1 2 times ... 
        matrix = []
        for j in range(3):
            counts = {0:0, 1:0, 2:0}
            counts.update(dict( Counter( [state_seq[i] for i in filter(lambda x: x < len(state_seq),\
                        map(lambda x: x+1, np.where(state_seq == j)[0])) ] ) ))
            if sum(counts.values()):  
                v = map(lambda x: float(x) / sum(counts.values()), counts.values() )
            else:
                v = [1.0 / self.n_components] * self.n_components 
            matrix.append(v)
        self._set_transmat(np.array(matrix))
        
        
        #init emission probabilities
        matrix = []
        for state in state_seq:
            row = [0] * self.n_components
            row[state] = 1
            matrix.append(row)
        posteriors = np.array(matrix)
        stats = self._initialize_sufficient_statistics()
        obs = map(lambda x: (x[1], x[2]), self.init_state_seq) 
        self._help_accumulate_sufficient_statistics(obs, stats, posteriors)
        self._help_do_mstep(stats)
        
    else:
        super(BinomialHMM2d3s, self)._init(obs, params)
    
    print("initital HMM:", file=sys.stderr)
    print("start probabilities", file=sys.stderr)
    print(self.startprob_, file=sys.stderr)
    print("transistion matrix", file=sys.stderr)
    print(self.transmat_, file=sys.stderr)
    print("emission distributin parameters", file=sys.stderr)
    print(self._get_emissionprob(), file=sys.stderr)



def _add_pseudo_counts(arr):
    if type(arr) is np.ndarray:
        tmp = np.array([1e-323 if x < 1e-323 else x for x in arr], np.float64)
        #tmp2 = np.array([1.0 - 1.0e-5 if x == 1.0 else x for x in tmp], np.float64)
        return tmp
    else:
        tmp = 1e-323 if arr < 1e-323 else arr
        #tmp2 = 1.0 - 1.0e-10 if tmp == 1.0 else tmp
        return tmp

def sub_neginf(m):
    return 1e-323 if m == -np.inf else m

#    def get_posteriors(self, obs):
#        """other version with more output"""
#        print('start calculation postertiors', file=sys.stderr)
#        
#        obs = np.asarray(obs)
#        framelogprob = self._compute_log_likelihood(obs)
#        
#        
#        t=time()
#        logprob, fwdlattice = self._do_forward_pass(framelogprob)
#        print('forward step: ', time()-t, file=sys.stderr)
#        
#        t=time()
#        bwdlattice = self._do_backward_pass(framelogprob)
#        print('backward-step: ', time()-t, file=sys.stderr)
#        
#        t=time()
#        gamma = fwdlattice + bwdlattice
#        # gamma is guaranteed to be correctly normalized by logprob at
#        # all frames, unless we do approximate inference using pruning.
#        # So, we will normalize each frame explicitly in case we
#        # pruned too aggressively.
#        posteriors = np.exp(gamma.T - logsumexp(gamma, axis=1)).T
#        posteriors += np.finfo(np.float32).eps
#        posteriors /= np.sum(posteriors, axis=1).reshape((-1, 1))
#        print('calculation: ', time()-t, file=sys.stderr)
#        return posteriors



class histogram():
    hist = {}
    density = {}
    
    def load(self, path):
        with open(path) as f:
            for line in f:
                line.strip()
                line = line.split('\t')
                self.set_value(int(line[0]), float(line[1]))
    
    def combine(self, h1, h2):
        keys = h1.density.keys() + h2.density.keys()
        for k in keys:
            v1 = h1.density[k] if h1.density.has_key(k) else 0
            v2 = h2.density[k] if h2.density.has_key(k) else 0
            self.set_value(k, (v1+v2)/2.0)
    
        self._normalize()
    
    def get_random_data(self):
        r = randint(0, 1000)
        s = 0
        for i in range(0, self.get_max_key()+1):
            s += self.get_prob_density(i) * 1000
            if s > r:
                return i
    
        while i is None:
            i = self.get_random_data()
        return i
    
    def _normalize(self):
        t = float(sum(self.density.values()))
        for i in self.density.keys():
            self.density[i] = self.density[i] / t

    def merge(self, h):
        key_set = self.density.keys() + h.density.keys()
        key_set = set(key_set)
        for k in key_set:
            if self.density.has_key(k) and h.density.has_key(k):
                self.density[k] = (self.density[k] + h.density[k]) / 2
            if not self.density.has_key(k) and h.density.has_key(k):
                self.density[k] = h.density.has_key(k) /2
            if self.density.has_key(k) and not h.density.has_key(k):
                self.density[k] = self.density[k] / 2
        
        self._normalize()
#        key_set = self.get_keys() + h.get_keys()
#        key_set = set(key_set)
#        for k in key_set:
#            if self.hist.has_key(k) and h.hist.has_key(k):
#                self.hist[k] = (self.hist[k] + h.hist[k]) / 2
#            if not self.hist.has_key(k) and h.hist.has_key(k):
#                self.hist[k] = h.hist.has_key(k) /2
#            if self.hist.has_key(k) and not h.hist.has_key(k):
#                self.hist[k] = self.hist[k] / 2
    
#        self.get_density()
    
    def __init__(self):
        self.hist = {}
        self.density = {}
    
    def __str__(self):
        k = self.hist.keys()
        k.sort()
        for i in k:
            print(i, self.hist[i], sep='\t')
        return ''
    
    def add(self, x):
        """Increment value for x by 1"""
        if self.hist.has_key(x):
            self.hist[x] += 1.
        else:
            self.hist[x] = 1.
        self._get_density()
        self._normalize()

    def set_value(self, x, y):
        """set density value x to y"""
        self.density[x] = y
    
    def get_max_key(self):
        """get highest key in histogram, that is the rightmost x-axis value higher than the pseudo count"""
        return max(self.density.keys())
    
    def get_integral(self):
        """return integral of histogram, do not consider pseudo counts"""
        if sum(self.hist.values()) < EPSILON:
            pass
        return sum(self.hist.values())
    
#    def remove(self, x):
#        """decrement value x by 1"""
#        if self.hist.has_key(x):
#            self.hist[x] -= 1
#            if self.hist[x] <= 0:
#                del self.hist[x]
    
    def _get_density(self):
        for k in self.hist.keys():
            self.density[k] = self._get_prob(k)
    
    
    def print_density(self, *filename):
        if len(filename) == 1:
            h = open(filename[0], 'w')
#        self.get_density()
        for k in self.density.keys():
            if len(filename) == 1:
                print(k, self.density[k], sep='\t', file=h)
            else:
                print(k, self.density[k], sep='\t')
    
    
    def smooth(self, max):
         r = 2
         #mirror at the egdes
         f = [self.density[i] if self.density.has_key(i) else 0. for i in range(max)]
         for i in range(len(f)):
             if f[i] < EPSILON:
                 if i < r: #at a edge
                     a = f[i:r+1] * 2
                 else:
                     a = f[i-r:i] + f[i+1:i+r+1]
                 f[i] = sum(a) / float(len(a))
         
         for i in range(len(f)):
             self.density[i] = f[i]
        
         self._normalize() 
#         t = float(sum(self.density.values()))
#         for i in self.density.keys():
#            self.density[i] = self.density[i] / t
         
#         self.get_density()
#         g = map(lambda x: x<EPSILON, f)
#         if g.count(True):
#             self.smooth(max)
    def delete_histogram(self):
        self.hist = {}
    
    def get_prob_density(self, x):
        """return probability of x or pseudo count"""
#        self.get_density()
        if self.density.has_key(x):
            return self.density[x]
        else:
            return EPSILON
        
    def _get_prob(self, x):
        """return probability of x or pseudo count"""
        if self.hist.has_key(x):
            return float(self.get(x))/self.get_integral()
        else:
            return EPSILON
    
    def get(self, x):
        """return value of x or pseudo count"""
        if self.hist.has_key(x):
            return self.hist[x]
        else:
            return EPSILON
    
#    def get_keys(self):
#        """get keys, that is values of the x-axis without pseudo counts"""
#        return self.hist.keys()
    
#    def save(self, filename):
#        """Save in filename"""
#        h = open(filename, 'w')
#        for k in self.hist.keys():
#            print(k, self.hist[k], sep='\t', file=h)
#        h.close
    