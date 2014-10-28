#!/usr/bin/env python
# -*- coding: utf-8 -*-

class Tracker:
    def __init__(self, path):
        self.file = open(path, 'w')
    
    def write(self, text=t, header=None):
        if header is not None:
            self.file.write('#' + header + '\n')
        else:
            self.file.write(t + '\n')