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


## actually here, we also need to get mean_var distribution
"""

from __future__ import print_function
import sys
import numpy as np
from random import sample
from math import fabs
from scipy import sparse, optimize

import configuration

import warnings
warnings.filterwarnings('error')

import matplotlib as mpl
#see http://stackoverflow.com/questions/4931376/generating-matplotlib-graphs-without-a-running-x-server
mpl.use('Agg')
import matplotlib.pyplot as plt


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
    
    return posteriors


def _func_quad_2p(x, a, c):
    """Return y-value of y=max(|a|*x^2 + x + |c|, 0),
    x may be an array or a single float"""
    res = []
    if type(x) is np.ndarray:
        for el in x:
            res.append(max(el, fabs(a) * el ** 2 + el + fabs(c)))

        return np.asarray(res)
    else:
        return max(x, fabs(a) * x ** 2 + x + fabs(c))


def _write_emp_func_data(data, name):
    """Write mean and variance data"""
    assert len(data[0]) == len(data[1])
    f = open(configuration.FOLDER_REPORT_DATA + name + '.data', 'w')
    for i in range(len(data[0])):
        print(data[0][i], data[1][i], sep='\t', file=f)
    f.close()


def _plot_func(plot_data, outputdir):
    """Plot estimated and empirical function"""

    maxs = []  # max for x (mean), max for y (var)
    for i in range(2):
        tmp = np.concatenate((plot_data[0][i], plot_data[1][i]))  # plot_data [(m, v, p)], 2 elements
        maxs.append(max(tmp[tmp < np.percentile(tmp, 90)]))

    for i in range(2):
        x = np.linspace(0, max(plot_data[i][0]), int(np.ceil(max(plot_data[i][0]))))
        y = _func_quad_2p(x, plot_data[i][2][0], plot_data[i][2][1])

        for j in range(2):
            # use matplotlib to plot function and datapoints
            # and save datapoints to files
            ext = 'original'
            if j == 1:
                plt.xlim([0, maxs[0]])
                plt.ylim([0, maxs[1]])
                ext = 'norm'
            ax = plt.subplot(111)
            plt.plot(x, y, 'r', label='fitted polynomial')  # plot polynom
            plt.scatter(plot_data[i][0], plot_data[i][1], label='empirical datapoints')  # plot datapoints
            ax.legend()
            plt.xlabel('mean')
            plt.ylabel('variance')
            plt.title('Estimated Mean-Variance Function')
            name = "_".join(['mean', 'variance', 'func', 'cond', str(i), ext])
            _write_emp_func_data(plot_data[i], name)
            plt.savefig(configuration.FOLDER_REPORT_PICS + name + '.png')
            plt.close()


def sm_var(sm, axis=None):
    """get variance of one sparse matrix in axis
    if axis=None, then for whole data
    if axis=0, for column
    if axis=1 for row
    """
    # get shape of sm_var and axis
    if axis is None:
        tmp = sm.data **2
        return np.mean(tmp) - (sm.mean()) ** 2
    elif axis == 0:
        # in column order, but now it's 1_D array
        tmp = sm[:,:] # create a new sparse-matrix
        tmp.data **= 2
        e_x2 = tmp.mean(axis=0)
        return e_x2 - np.square(sm.mean(axis=0))


def _get_sm_data_rep(one_sample_cov, sample_size):
    """Return list of (mean, var) points for samples 0 and 1
    overall_coverage is a list of sparse matrix
    """
    # data_rep = []
    # firstly to get union of non-zeros columns
    tmp_cov = sum(one_sample_cov)
    # sample indices and then sample it
    idx = sample(tmp_cov.indices, min(sample_size,len(tmp_cov.indices)))
    idx.sort()
    cov = sparse.vstack([one_sample_cov[j][:,idx] for j in range(len(one_sample_cov))])
    # count the mean and var of it
    ms = cov.mean(axis=0)
    vars = sm_var(cov,axis=0)
    ms = np.insert(np.squeeze(np.asarray(ms)), len(idx),0)
    vars = np.insert(np.squeeze(np.asarray(vars)),len(idx),0)
    # we add 0 to ms and vars, but consider computation time, could we find a better way??
    idx = np.logical_and(ms < np.percentile(ms, 99.75) * (ms > np.percentile(ms, 1.25)),
                      vars < np.percentile(vars, 99.75) * (vars > np.percentile(vars, 1.25)))
    final_ms = ms[idx]
    final_vars = vars[idx]

    # we need to return mean and var list, maybe one by one, we could filter data
    return final_ms, final_vars


def fit_sm_mean_var_distr(overall_coverage, name, debug, verbose, outputdir, report, poisson, sample_size=5000):
    """Estimate empirical distribution (quadr.) based on empirical distribution
    On paper, it says for smaller samples, we use this, but for big samples, should we change methods ???
    change this method and combine it with get_data_rep;
    main thing is repeat to get data
    """
    plot_data = []  # means, vars, paras
    dim = overall_coverage['dim'] # now it's a list, so we can't do it , but we could use dictionary and add dimension into it
    done = False
    # for this function, if we get bad result, we do it again, to get better samples again
    # and fit it again...until a good fit
    while not done:
        res = []
        for i in range(dim[0]):
            try:
                m, v = _get_sm_data_rep(overall_coverage['data'][i], sample_size)

                if debug:
                    np.save(str(name) + "-emp-data" + str(i) + ".npy", zip(m,v))

                if len(m) > 0 and len(v) > 0:
                    try:
                        p, _ = optimize.curve_fit(_func_quad_2p, m, v)  # fit quad. function to empirical data
                    except:
                        print("Optimal parameters for mu-var-function not found, get new datapoints", file=sys.stderr)
                        break  # restart for loop
                else: # one length is zero, it means??
                    p = np.array([0, 1])

                res.append(p)
                plot_data.append((m, v, p))

            except RuntimeError:
                print("Optimal parameters for mu-var-function not found, get new datapoints", file=sys.stderr)
                break  # restart for loop
        # after loop successfully ends, we are done
        done = True
    if report:
        _plot_func(plot_data, outputdir)

    if poisson:
        print("Use Poisson distribution as emission", file=sys.stderr)
        p[0] = 0
        p[1] = 0 # here it should be 0, we assume var = mean
        res = [np.array([0, 0]), np.array([0, 0])]

    return lambda x: _func_quad_2p(x, p[0], p[1]), res

