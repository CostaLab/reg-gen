#!/usr/bin/env python
# -*- coding: utf-8 -*-

class Tracker:
    def __init__(self, path):
        self.file = open(path, 'w')
    
    def write(self, text, header=None):
        if header:
            self.file.write('#' + header + '\n')
        
        self.file.write(text + '\n')