#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import re
import numpy as np
from rgt.Util import Html, ConfigurationFile
from collections import OrderedDict
from os import path
import sys

class Tracker:
    data = []
    def __init__(self, path):
        self.file = open(path, 'w')
    
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
        
    def make_html(self):
        html_header = "THOR run"
        from rgt.THOR.dpc_help import FOLDER_REPORT
        #Links
        links_dict = OrderedDict()
        # XXX only if path exists
        links_dict['Run Information'] = 'index.html#runinfo'
        links_dict['Mean Variance Function Estimate'] = 'index.html#mvfunction'
        links_dict['Fragment Size Estimate'] = 'index.html#fsestimate'
        
        config_class = ConfigurationFile()
        html = Html(name=html_header, links_dict=links_dict, fig_rpath= config_class.data_dir + '/fig/')
        
        #Run Info
        html.add_heading("Run Information", idtag = 'runinfo')
        html.add_list([el[0] + '<br>' + el[1] for el in self.data])
        
        #Mean Variance Function
        p = path.join(FOLDER_REPORT, "pics/mean_variance_func_cond_0_original.png")
        if path.isfile(p):
            html.add_heading("Mean Variance Function", idtag='mvfunction')
            html.add_figure(p, align="left", width="45%", more_images=[path.join(FOLDER_REPORT, 'pics/mean_variance_func_cond_1_original.png')])
        
        #Fragment Size Estimate
        p = path.join(FOLDER_REPORT, 'pics/fragment_size_estimate.png')
        if path.isfile(p):
            html.add_heading("Fragment Size Estimate", idtag = 'fsestimate')
            html.add_figure(p, align="left", width="45%")
        
        
        
        
        html.write(path.join(FOLDER_REPORT, "index.html"))
        
        
        