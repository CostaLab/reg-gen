#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ODIN detects differential peaks in multiple ChIP-seq profiles associated
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
