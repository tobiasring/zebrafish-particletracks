# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 21:08:53 2018

@author: ring
"""

from matplotlib import use
use('Qt4Agg')
import numpy as np
import matplotlib.pyplot as plt

def intersect(p1, p2, p3, p4):
    A = np.array([[p2[0], -p4[0]], [p2[1], -p4[1]]])
    b = np.array([p3[0] - p1[0], p3[1] - p1[1]])
    [l,g] = np.linalg.solve(A, b)
    intersec = [p1[0] + l*p2[0], p1[1] + l*p2[1]]
    # print 'Schnittpunkt: x = %.1f, y = %.1f' %(intersec[0], intersec[1])
    return intersec
    
if __name__ == '__main__':
    a = [0,0]
    b = [1,1]
    c = [-1,1]
    d = [2,-2]
    intersect(a, b, c, d)