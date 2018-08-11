#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 11 22:36:50 2018

@author: tobias
"""

import os

print os.walk('all exports for calculation').next()

'''
for path, directories, files in os.walk('./'):
     if file in files:
          print('found %s' % os.path.join(path, file))
'''