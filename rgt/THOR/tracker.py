#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import numpy as np

class Tracker:
    def __init__(self, path):
        self.file = open(path, 'w')
    
    def write(self, text, header=None):
        if header:
            self.file.write('#' + header + '\n')
        if isinstance(text, np.ndarray):
            self.file.write(re.sub(' +',' ', str(text).replace('\n ', '\n').replace('[','').replace(']','').strip().replace('\n ', '\n').replace('\n ', '\n').replace('\n ', '\n')) + '\n')
        else:
            if isinstance(text, list):
                text = ",".join(text)
            if isinstance(text, int) or isinstance(text, float):
                text = str(text)
            self.file.write(text + '\n')