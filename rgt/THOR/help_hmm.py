#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
THOR detects differential peaks in multiple ChIP-seq profiles between
two distinct biological conditions.

Copyright (C) 2014-2016 Manuel Allhof (allhoff@aices.rwth-aachen.de)

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


import sys
import numpy as np
from math import fabs
import warnings


def _count(posts):
    c_1, c_2 = 1, 1
    
    for s0, s1, s2 in posts:        
        if s0 > 0.5:
            c_1 += 0
            c_2 += 0
        elif s1 >= s2:
            c_1 += 1
        elif s2 > s1:
            c_2 += 1
    return c_1, c_2

def _valid_posteriors(posteriors, obs, dim):
    
    warnings.filterwarnings('error')
    
    for i in range(len(obs)):
        state_1 = False
        c1, c2 = np.mean(obs[i][:dim[0]]), np.mean(obs[i][dim[0]:]) #counts of samples

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

    warnings.resetwarnings()

    return posteriors

