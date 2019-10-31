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


import re
import numpy as np
from ..Util import Html
from collections import OrderedDict
from os import path
# import sys
from datetime import datetime


class Tracker:
    data = []

    def __init__(self, p, bamfiles, genome, chrom_sizes, dims, inputs, options, version):
        self.file = open(p, 'w')
        self.bamfiles = bamfiles
        if genome is None:
            self.genome = "None"
        else:
            self.genome = genome
        self.chrom_sizes = chrom_sizes
        self.dims = dims
        self.inputs = inputs
        self.samples = [path.splitext(path.basename(x))[0] for x in bamfiles]
        self.options = options
        self.version = version
    
    def write(self, text, header):
        if header:
            self.file.write('#' + header + '\n')
        if isinstance(text, np.ndarray):
            text = re.sub(' +',' ', str(text).replace('\n ', '\n').replace('[','').replace(']','').strip().\
                          replace('\n ', '\n').replace('\n ', '\n').replace('\n ', '\n'))
        else:
            if isinstance(text, list):
                text = " ".join(text)
            if isinstance(text, int) or isinstance(text, float):
                text = str(text)
            
        self.file.write(text + '\n')
        
        tmp = text.replace('\n', ' ').split(' ')
        if len(tmp) == 1:
            text = tmp[0]
        else:
            new = []
            for i, el in enumerate(tmp):
                i += 1
                new.append(el)
                if i > 0 and i % 3 == 0 and i < len(tmp) and len(tmp) % 3 == 0:
                    new.append('<br>')
                
            text = " ".join(new)
        self.data.append((header, text))

    def _read_hk(self, p):
        d = []
        if path.isfile(p):
            with open(p) as f:
                for line in f:
                    line = line.split(" ")
                    d.append([line[0], str(float(line[1]))])
        return d
    
    def make_hmm(self, html):
        d_all = []
        for el in self.data:
            if not el[0].startswith('Neg') and not el[0].startswith('Parameters') and\
            not el[0].startswith('Ext') and not el[0].startswith('Scaling'):
                d_all.append((el[0], el[1]))
        
        html.add_list([el[0] + '<br>' + el[1] for el in d_all])
        
        info = "THOR uses an HMM to segment the ChIP-seq signals. The HMM's emission is a Negative Binomial distribution. Inital\
        values (before HMM training) and final values (after HMM training) of the parameters mu and alpha are given. Furthermore,\
        the transition matrix is given. The rows and columns of the matrix describe: background state, DP found in first \
        condition and DP found in second condition. See [1] for further details about the HMM training."
        self._write_text(html, info)
        
    def make_pre_post(self, html):
        d_all = []
        for el in self.data:
            if el[0].startswith('Neg') or el[0].startswith('Parameters'):
                d_all.append((el[0], el[1]))
        
        html.add_list([el[0] + '<br>' + el[1] for el in d_all])
        
        info = "THOR's p-value calculation is based on a Negative Binomial distribution with \
        parameter mu which gives the location and paramter alpha which gives the dispersion. <br>Moreover, THOR uses a\
        polynomial function of degree 2 with parameters <i>a</i> and <i>c</i> to empirically describe the mean-variance relationship in\
        the data. <br><b>Note:</b> Parameter <i>a</i> describes the function's quadratic characteristic and can be used as indicator \
        for overdispersion. In our experience, for experiments with high <i>a</i> (higher than 1e-2) it is better\
        to normalize with housekeeping genes (see [1] for more details)."
        self._write_text(html, info)
    
    def make_ext_scaling_table(self, html):
        """make table: sample, ext, scaling"""
        d_all = []
        for el in self.data:
            if el[0].startswith('Ext'):
                d_all.append(el[1])
            if el[0].startswith('Scaling'):
                d_all.append(el[1])
        
        exts = d_all[0].split(' ')
        factors = d_all[1].split(' ')
        
        l = []
                
        for i in range(len(exts)):
            l.append([self.samples[i], exts[i], factors[i]])
        
        html.add_zebra_table(header_list=['Sample', 'Extension Size', 'Scaling Factor'], col_size_list=[1,150,150], type_list='sss', data_table=l, auto_width = True)
        
        info = "Since the Chip-seq protocol makes it necessary to extend the experiment's reads, we follow [2] and define a strand cross-correlation function to derive the\
        extension size from it.<br> To bring all ChIP-seq profiles to the same scale, THOR normalizes each sample with a scaling factor. \
        The scaling factor is either based on TMM [3] or the house-keeping gene approach [1]."
        self._write_text(html, info)
  
    def make_ext_config(self, html):
        """make table about configuration: feature, path"""
        b = []
        for el in self.bamfiles:
            b += [el, '<br>']
        
        b = b[:len(b) - 1]
        b = " ".join(b)
        
        if self.options.housekeeping_genes:
            norm = 'Housekeeping Genes approach'
        elif self.options.scaling_factors_ip:
            norm = 'Predefined Values'
        else:
            norm = 'TMM'
        
        d = [['Name', path.splitext(path.basename(self.options.name))[0]], ['time', str(datetime.now()).split(".")[0]], ['BAM files', b], ['genome', self.genome], ['chromosome sizes', self.chrom_sizes],\
             ['Normalization Strategy', norm], ['version', self.version], ['merge DPs', str(self.options.merge)], ['p-value cutoff', str(self.options.pcutoff)],\
             ['deadzones', str(self.options.deadzones)]]
        
        if self.inputs:
            a = []
            for el in self.inputs:
                a += [el, '<br>']
            
            a = a[:len(a) - 1]
            a = " ".join(a)
        
            d.append(['BAM input-DNA files', a])
        
        html.add_zebra_table(header_list=['Feature', ''], col_size_list=[190, 1000], type_list='ss', data_table=d, cell_align='left')
    
    def _write_text(self, html, t):
        r = "<table style=\"width:90%\"  align=\"center\"  cellpadding=\"0\" cellspacing=0>"
        r += "<tr><td>"
        #t = t.replace('<br>', '<br>' + '&nbsp;'*9)
        #t = t.replace('<p>', '<p>' + '&nbsp;'*9)
        #t = '<p>' + '&nbsp;'*10 + t + '<br>'
        r += t
        r += "</td></tr>"
        r += "</table>"
        html.add_free_content([r])
    
    def make_html(self):
        html_header = "THOR"
        from rgt.THOR.dpc_help import FOLDER_REPORT

        #Links
        links_dict = OrderedDict()
        links_dict['Experimental Configuration'] = 'index.html#extinfo'
        links_dict['Sample Information'] = 'index.html#sampleinfo'
        links_dict['HMM Information'] = 'index.html#hmminfo'
        links_dict['Mean Variance Function Estimate'] = 'index.html#mvfunction'
        
        p = path.join(FOLDER_REPORT, 'pics/fragment_size_estimate.png')
        if path.isfile(p):
            links_dict['Fragment Size Estimate'] = 'index.html#fsestimate'

        p = path.join(FOLDER_REPORT, 'pics/data/sample.data')
        if path.isfile(p):
            links_dict['Housekeeping Gene Normalization'] = 'index.html#norm'
        
        links_dict['References'] = 'index.html#ref'
        links_dict['Contact'] = 'index.html#contact'

        # copy basic rgt logo, style etc to local directory inside report
        fig_path = path.join(FOLDER_REPORT, "fig")
        html = Html(name=html_header, links_dict=links_dict, fig_dir=fig_path, fig_rpath="fig")
        
        try:
            html.add_heading("Experimental Configuration", idtag='extinfo')
            self.make_ext_config(html)
        except:
            pass

        html.add_heading("Pre- and post-processing Features", idtag='prepostinfo')
        self.make_pre_post(html)
        
        try:
            html.add_heading("Sample Information", idtag='sampleinfo')
            self.make_ext_scaling_table(html)
        except:
            pass

        #Run Info
        try:
            html.add_heading("HMM Information", idtag='hmminfo')
            self.make_hmm(html)
        except:
            pass

        #Mean Variance Function
        try:
            p = path.join(FOLDER_REPORT, 'pics/mean_variance_func_cond_0_original.png')
            if path.isfile(p):
                html.add_heading("Mean Variance Function", idtag='mvfunction')
                html.add_figure(path.relpath(p, FOLDER_REPORT), align="left", width="45%",
                                more_images=['pics/mean_variance_func_cond_1_original.png'])
                info = "THOR uses a polynomial function to empirically describe the relationship between mean and variance in the data.\
                The data the plot is based on can be found at report/pics/data for further downstream analysis."
                self._write_text(html, info)
        except:
            pass

        #Fragment Size Estimate
        try:
            p = path.join(FOLDER_REPORT, 'pics/fragment_size_estimate.png')
            if path.isfile(p):
                html.add_heading("Fragment Size Estimate", idtag='fsestimate')
                html.add_figure(path.relpath(p, FOLDER_REPORT), align="left", width="45%")
                info = "THOR estimates the fragmentation sizes of each sample's reads. Here, the cross-correlation function [1] is shown. Their maxima give the\
                fragmentation extension sizes.<br> The data the plot is based on can be found at report/pics/data for further downstream analysis."
                self._write_text(html, info)
        except:
            pass

        #HK normalization
        try:
            p = path.join(FOLDER_REPORT, 'pics/data/gene.data')
            if path.isfile(p):
                d = self._read_hk(p)
                html.add_heading("Housekeeping Gene Normalization", idtag='norm')
                html.add_zebra_table(header_list=['gene', 'quality q'], col_size_list=[1,150], type_list='s'*len(d), data_table=d)
                info = "For active histone marks, housekeeping genes given by [4] can be used for normalization [1]. Here, the genes for the experiments are\
                evaluated. For each gene i, we estimate the normalization factors with gene i and without gene i and compute the sums of squared deviations q.\
                High values (higher than 2) indicate striking genes which should be considered to be left our for normalization.,<br> One can also \
                use other genes or regions for normalization.<br> The data the plot is based on can be found at report/pics/data for further downstream analysis."
                self._write_text(html, info)
                
            p = path.join(FOLDER_REPORT, 'pics/data/sample.data')
            if path.isfile(p):
                d = self._read_hk(p)
                html.add_zebra_table(header_list=['sample', 'quality p'], col_size_list=[1,150], type_list='s'*len(d), data_table=d)
                info = "We evaluate the effect of samples to the normalization factors. For sample j, we estimate the normalization factors with sample j\
                and without sample j and compute the sums of squared deviations p. High values (higher than 2) indicate striking samples which should be\
                considered to be left out for the analysis.<br> The data the plot is based on can be found at report/pics/data for further downstream analysis."
                self._write_text(html, info)
        except:
            pass

        html.add_heading("References", idtag='ref')
        info = "[1] M. Allhoff, J. F. Pires, K. Ser&eacute;, M. Zenke, and I. G. Costa. Differential Peak Calling of ChIP-Seq \
        Signals with Replicates with THOR. <i>submitted.</i> <br>\
        [2] A. Mammana, M. Vingron, and H.-R. Chung. Inferring nucleosome positions with their histone mark annotation from chip data. \
        Bioinformatics, 29(20):2547-2554, 2013. <br>\
        [3] M. D. Robinson and A. Oshlack. A scaling normalization method for differential expression analysis of RNA-seq data. \
        Genome Biology, 11(3):R25, 2010. <br>\
        [4] E. Eisenberg and E. Y. Levanon. Human housekeeping genes, revisited. Trends in genetics: TIG, 29(10):569-574, 2013."
        self._write_text(html, info)
        
        html.add_heading("Contact", idtag='contact')
        info = "If you have any questions, please don't hesitate to contact us: allhoff@aices.rwth-aachen.de"
        self._write_text(html, info)
        
        html.write(path.join(FOLDER_REPORT, "index.html"))
