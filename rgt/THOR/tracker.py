#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import re
import numpy as np
from rgt.Util import Html, ConfigurationFile
from collections import OrderedDict
from os import path
import sys
from datetime import datetime

class Tracker:
    data = []
    def __init__(self, p, bamfiles, genome, chrom_sizes, dims, inputs, options):
        self.file = open(p, 'w')
        self.bamfiles = bamfiles
        self.genome = genome
        self.chrom_sizes = chrom_sizes
        self.dims = dims
        self.inputs = inputs
        self.samples = map(lambda x: path.splitext(path.basename(x))[0], bamfiles)
        self.options = options
    
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
        
        info = "this is a info text.<br>this should be a new line test."
        self._write_text(html, info)
        
    def make_pre_post(self, html):
        d_all = []
        for el in self.data:
            if el[0].startswith('Neg') or el[0].startswith('Parameters'):
                d_all.append((el[0], el[1]))
        
        html.add_list([el[0] + '<br>' + el[1] for el in d_all])
        
        info = "this is a info text.<br>this should be a new line test."
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
        
        html.add_zebra_table(header_list=['Sample', 'Extension Size', 'Scaling Factor'], col_size_list=[1,1,1], type_list='sss', data_table=l, auto_width = True)
        
        info = "this is a info text.<br>this should be a new line test."
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
            norm = TMM
        
        d = [['Name', options.name], ['time', str(datetime.now())], ['BAM files', b], ['genome', self.genome], ['chromosome sizes', self.chrom_sizes],\
             ['Normalization Strategy', norm], ['version', self.options.version], ['merge DPs', self.options.merge], ['p-value cutoff', self.options.pcutoff],\
             ['deadzones', self.options.deadzones]]
        
        if self.inputs:
            a = []
            for el in self.inputs:
                a += [el, '<br>']
            
            a = a[:len(a) - 1]
            a = " ".join(a)
        
            d.append(['BAM input-DNA files', a])
        
        html.add_zebra_table(header_list=['Feature', 'Value'], col_size_list=[190, 1000], type_list='ss', data_table=d, cell_align='left')
    
    def _write_text(self, html, t):
        t = t.replace('<br>', '<br>' + '&nbsp;'*10)
        t = t.replace('<p>', '<p>' + '&nbsp;'*10)
        t = '<p>' + '&nbsp;'*10 + t + '<br>'
        html.add_free_content([t])
    
    def make_html(self):
        html_header = "THOR run"
        from rgt.THOR.dpc_help import FOLDER_REPORT
        #Links
        links_dict = OrderedDict()
        links_dict['Experimental Configuration'] = 'index.html#extinfo'
        links_dict['Sample Information'] = 'index.html#sampleinfo'
        links_dict['HMM Information'] = 'index.html#runhmminfo'
        links_dict['Mean Variance Function Estimate'] = 'index.html#mvfunction'
        
        
        
        p = path.join(FOLDER_REPORT, 'pics/fragment_size_estimate.png')
        if path.isfile(p):
            links_dict['Fragment Size Estimate'] = 'index.html#fsestimate'
        
        p = path.join(FOLDER_REPORT, 'pics/data/sample.data')
        if path.isfile(p):
            links_dict['Housekeeping Gene Normalization'] = 'index.html#norm'
        
        config_class = ConfigurationFile()
        html = Html(name=html_header, links_dict=links_dict, fig_rpath= config_class.data_dir + '/fig/')
        
        #try:
        html.add_heading("Experimental Configuration", idtag = 'extinfo')
        self.make_ext_config(html)
        #except:
        #    pass
        
        html.add_heading("Pre/Postprocessing Features", idtag = 'prepostinfo')
        self.make_pre_post(html)
        
        try:
            html.add_heading("Sample Information", idtag = 'sampleinfo')
            self.make_ext_scaling_table(html)
        except:
            pass
        
        #Run Info
        try:
            html.add_heading("HMM Information", idtag = 'hmminfo')
            self.make_hmm(html)
        except:
            pass
        
        #Mean Variance Function
        try:
            p = path.join(FOLDER_REPORT, "pics/mean_variance_func_cond_0_original.png")
            if path.isfile(p):
                html.add_heading("Mean Variance Function", idtag='mvfunction')
                html.add_figure(p, align="left", width="45%", more_images=[path.join(FOLDER_REPORT, 'pics/mean_variance_func_cond_1_original.png')])
                info = "this is a info text.<br>this should be a new line test."
                self._write_text(html, info)
        except:
            pass
        
        #Fragment Size Estimate
        try:
            p = path.join(FOLDER_REPORT, 'pics/fragment_size_estimate.png')
            if path.isfile(p):
                html.add_heading("Fragment Size Estimate", idtag = 'fsestimate')
                html.add_figure(p, align="left", width="45%")
        except:
            pass
        
        #HK normalization
        try:
            p = path.join(FOLDER_REPORT, 'pics/data/gene.data')
            if path.isfile(p):
                d = self._read_hk(p)
                html.add_heading("Housekeeping Gene Normalization", idtag = 'norm')
                html.add_zebra_table(header_list=['gene', 'quality q'], col_size_list=[1,1], type_list='s'*len(d), data_table=d, auto_width = True)
            
            p = path.join(FOLDER_REPORT, 'pics/data/sample.data')
            if path.isfile(p):
                d = self._read_hk(p)
                html.add_zebra_table(header_list=['sample', 'quality p'], col_size_list=[1,1], type_list='s'*len(d), data_table=d, auto_width = True)
        except:
            pass
        
        html.write(path.join(FOLDER_REPORT, "index.html"))
        
        
        