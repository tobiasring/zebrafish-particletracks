# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 22:08:02 2018

@author: ring
"""
import matplotlib.pyplot as plt

import numpy as np
import glob
import os
import csv
import sys

def getParticleTracks(path2file, minTrackLength = 1):
    all_tracks = []
    # read reference line information
    refline = []
    if os.path.exists(path2file + 'axis points.csv'):
        with open(path2file + 'axis points.csv', 'rb') as axisfile:
            reader = csv.reader(axisfile)
            for row in reader:
                if row[5] != 'X':
                    x = float(row[5])
                    y = float(row[6])
                    refline.append([x, y])
    elif os.path.exists(path2file + 'axis points.txt'):
        with open(path2file + 'axis points.txt', 'r') as axisfile:
            for line in axisfile:
                if 'Max' not in line:
                    x = float(line.split()[5])
                    y = float(line.split()[6])
                    refline.append([x, y])
                
    # read particle track data
    data = np.genfromtxt(path2file + 'tracks.txt', delimiter = '\t', skip_header = 1, missing_values = 'NA')
    nObservs = int(np.max(data[:,0]))
    nTracks = int(np.max(data[:,1]))
    for trackID in range(1, nTracks + 1):
        thistrack = particle_track(path2file)
        thistrack.path = path2file
        thistrack.trackID = trackID
        thistrack.ref = refline
        thistrack.minLength = minTrackLength
        for j in range(nObservs):
            if data[j, 1] == trackID:
                thistrack.xpos.append(data[j, 3])
                thistrack.ypos.append(data[j, 4])
        if len(thistrack.xpos) >= minTrackLength:
            thistrack.compDataOfInterest()
            all_tracks.append(thistrack)
    return all_tracks

def intersect(p1, p2, p3, p4):
    # p1 and p2 define the first line
    # p3 and p4 define the second line
    A = np.array([[p2[0], -p4[0]], [p2[1], -p4[1]]])
    b = np.array([p3[0] - p1[0], p3[1] - p1[1]])
    [l,g] = np.linalg.solve(A, b)
    intersec = [p1[0] + l*p2[0], p1[1] + l*p2[1]]
    # print 'Intersection: x = %.1f, y = %.1f' %(intersec[0], intersec[1])
    return intersec
    
def pointLineDist(point, line):
    p = np.array([point[0], point[1], 0])
    ga = np.array([line[0][0], line[0][1], 0])
    gb = np.array([line[1][0] - line[0][0], line[1][1] - line[0][1], 0])
    return np.linalg.norm(np.cross(p - ga, gb))/np.linalg.norm(gb)
    
def saveTracks2File(path2file, listOfTracks):
    print 'Writing data to file in ASCII format...'
    fname = path2file.split(os.sep)[-2]
    with open(path2file + fname + '_trackDATA.txt', 'w') as datafile:
        datafile.write('# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n')
        datafile.write('# This is the particle track data of file: %s\n' %(fname + '/tracks.txt'))
        datafile.write('# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n\n')
        for track in listOfTracks:
            datafile.write('# --------------------------------------\n')
            datafile.write('# Data of Track with ID %d\n' %track.trackID)
            datafile.write('# --------------------------------------\n')
            if track.minLength == 1:
                datafile.write('#   * Number of track points: %d\n' %track.npos)
            else:
                datafile.write('#   * Number of track points: %d\n' %track.npos)
                datafile.write('#    -> tracks < %d points discarded!\n' %track.minLength)
            datafile.write('#   * Time interval in this track:  %.3f secs\n' %track.dt)
            datafile.write('#   * Total duration of this track: %.3f secs\n' %track.ges_t)
            datafile.write('#   * Reference line coordinates:\n')
            datafile.write('#         x = %5.2f, y = %5.2f\n' %(track.ref[0][0], track.ref[0][1]))
            datafile.write('#         x = %5.2f, y = %5.2f\n' %(track.ref[1][0], track.ref[1][1]))
            datafile.write('#   * Time running  [secs.]:  %10.6f\n' %track.t_running)
            datafile.write('#   * Time paused   [secs.]:  %10.6f\n' %track.t_paused)
            datafile.write('#   * Relative time paused:   %10.6f\n' %track.fraction_paused)
            datafile.write('#   * Mean velocity [my / s]: %10.6f\n' %track.vmean)
            datafile.write('#   * Max velocity  [my / s]: %10.6f\n' %track.vmax)
            datafile.write('#   * Max vel. (pos. dir.):   %10.6f\n' %track.vmax_pos)
            datafile.write('#   * Max vel. (neg. dir.):   %10.6f\n' %track.vmax_neg)
            datafile.write('#   * Net velocity  [my / s]: %10.6f\n' %track.vnet)
            datafile.write('#   * Runlengths [my]\n')
            datafile.write('#     - positive:\n')
            if len(track.pos_legs) != 0:
                for leg in track.pos_legs:
                    datafile.write('#            * %.6f \n' %leg)
            else:
                datafile.write('#            * none detected\n')
            datafile.write('#     - negative:\n')
            if len(track.neg_legs) != 0:
                for leg in track.neg_legs:
                    datafile.write('#            * %.6f \n' %leg)
            else:
                datafile.write('#            * none detected\n')
            datafile.write('# --------------------------------------\n\n')
    print '... finished.'
    return 0
    
class particle_track:
    def __init__(self, path):
        self.path = path
        self.dt = 0     # time between two points in seconds
        self.xpos = []      # x position in 1e-6 m
        self.ypos = []      # y position in 1e-6 m
        if os.path.exists(self.path + 'frame interval.txt'):
            with open(self.path + 'frame interval.txt', 'r') as dtfile:
                self.dt = float(dtfile.read())
        else:
            print 'No frame interval file found! Assuming 1.299 s'
            self.dt = 1.299
        return None
        
    def compDataOfInterest(self):
        self.npos = len(self.xpos)
        self.ges_t = self.npos * self.dt
        # compute velocities per path segment
        self.v = []
        self.ds = []
        self.dist2refline = []
        for i in range(1, self.npos):
            dx = np.abs(self.xpos[i] - self.xpos[i-1])
            dy = np.abs(self.ypos[i] - self.ypos[i-1])
            self.ds.append(np.sqrt(dx**2 + dy**2))
            v = self.ds[-1] / self.dt
            self.v.append(v)
        self.v = np.array(self.v)
        self.vmean = np.mean(self.v)
        n_paused = np.sum(np.where(self.v <= .19, 1, 0))
        self.t_paused = n_paused * self.dt
        self.t_running = (self.npos - n_paused) * self.dt
        self.fraction_paused = self.t_paused / self.ges_t

        '''
        alle V < 0.19 mym/s sind Pausen! aungeben, wieviel Zeit Pause, wieviel Zeit Bewegung!
        '''
        self.vmax = np.max(self.v)
        # compute net velocity
        dxnet = np.abs(self.xpos[-1] - self.xpos[0])
        dynet = np.abs(self.ypos[-1] - self.ypos[0])
        self.dsnet = np.sqrt(dxnet**2 + dynet**2)
        self.vnet = self.dsnet / (self.dt * self.npos)
        # compute Runlength
        self.xintersec = []
        self.yintersec = []
        Lambda = []
        for i in range(self.npos):
            # the following refers to the definition of the ref-line:
            # the second point is the start, the first point is the target
            refline_direction = [self.ref[0][0] - self.ref[1][0], self.ref[0][1] - self.ref[1][1]]
            refline_normal = [refline_direction[1], refline_direction[0]]
            # print refline_direction
            # print refline_normal
            intersection = intersect([self.xpos[i], self.ypos[i]], [self.xpos[i] + refline_normal[0], self.ypos[i] + refline_normal[1]], [self.ref[0][0], self.ref[0][1]], [self.ref[1][0], self.ref[1][1]])
            # print 'intersects at: x = %.3f and y = %.3f' %(intersection[0], intersection[1])
            if refline_direction[0] != 0:
                Lambda.append((intersection[0] - self.ref[0][0])/refline_direction[0])
            elif refline_direction[1] != 0:
                Lambda.append((intersection[1] - self.ref[0][1])/refline_direction[1])
            else:
                Lambda.append(0)
        indicator = np.sign([Lambda[i+1] - Lambda[i] for i in range(len(Lambda)-1)]) # np.diff(Lambda))
        self.pos_legs = []
        self.neg_legs = []
        counter = self.ds[0]
        for i in range(1, len(indicator)):
            if indicator[i] == indicator[i-1]:
                counter += self.ds[i]
            else:
                if indicator[i-1] == 1:
                    self.pos_legs.append(counter)
                elif indicator[i-1] == -1:
                    self.neg_legs.append(counter)
                else:
                    pass
                counter = self.ds[i]
        if indicator[-1] == -1:
            self.neg_legs.append(counter)
        elif indicator[-1] == 1:
            self.pos_legs.append(counter)
        else:
            pass
        vmaxIDX = np.where(self.v == np.max(self.v))[0][0]
        self.vmax = indicator[vmaxIDX] * self.vmax
        if self.pos_legs != []:
            self.vmax_pos = np.max(self.pos_legs) / self.dt
        else:
            self.vmax_pos = np.nan
        if self.neg_legs != []:
            self.vmax_neg = np.max(self.neg_legs) / self.dt
        else:
            self.vmax_neg = np.nan
        return 0

if __name__ == '__main__':
    '''
    The data to be computed has to be stored in the same level as this file.
    Put the name of the data directory in the DirName variable. 
    Everything else should work. I hope so at least. 
    '''

    minLength = 4
    DirName = 'testtracks' # 'all Cadherin Tracking data 091219'

    for (root,dirs,files) in os.walk(DirName, topdown=True):
        if dirs == []:
            path2file = root + os.sep
            data = getParticleTracks(path2file, minLength)
            saveTracks2File(path2file, data)
