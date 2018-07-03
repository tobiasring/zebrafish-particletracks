# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 22:08:02 2018

@author: ring
"""

import numpy as np
import glob
import os
import csv
import sys

# start script - change into working dir
os.chdir('Partikel Tracking 040816')

dirlist = glob.glob('./*')

def getParticleTracks(minTrackLength = 1):
    all_tracks = []
    # read reference line information
    refline = []
    with open('axis points.csv', 'rb') as axisfile:
        reader = csv.reader(axisfile)
        for row in reader:
            if row[5] != 'X':
                x = float(row[5])
                y = float(row[6])
                refline.append([x, y])
                
    # read particle track data
    data = np.genfromtxt('tracks.txt', delimiter = '\t', skip_header = 1, missing_values = 'NA')  
    nObservs = int(np.max(data[:,0]))
    nTracks = int(np.max(data[:,1]))
    for trackID in range(1, nTracks + 1):
        thistrack = particle_track()
        thistrack.trackID = trackID
        thistrack.ref = refline
        thistrack.minLength = minTrackLength
        for j in range(nObservs):
            if data[j, 1] == trackID:
                thistrack.xpos.append(data[j, 3])
                thistrack.ypos.append(data[j, 4])
        thistrack.compDataOfInterest()
        if len(thistrack.xpos) >= minTrackLength:
            all_tracks.append(thistrack)
    return all_tracks

def pointLineDist(point, line):
    p = np.array([point[0], point[1], 0])
    ga = np.array([line[0][0], line[0][1], 0])
    gb = np.array([line[1][0] - line[0][0], line[1][1] - line[0][1], 0])
    return np.linalg.norm(np.cross(p - ga, gb))/np.linalg.norm(gb)
    
def saveTracks2File(listOfTracks):
    print 'Writing data to file in ASCII format...'
    fname = os.getcwd().split('/')[-1].split('.')[0]
    with open(fname + '_trackDATA.txt', 'w') as datafile:
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
            datafile.write('#   * Reference line coordinates:\n')
            datafile.write('#         x = %5.2f, y = %5.2f\n' %(track.ref[0][0], track.ref[0][1]))
            datafile.write('#         x = %5.2f, y = %5.2f\n' %(track.ref[1][0], track.ref[1][1]))
            datafile.write('#   * Mean velocity [my / s]: %10.6f\n' %track.vmean)
            datafile.write('#   * Max velocity [my / s]:  %10.6f\n' %track.vmax)
            datafile.write('#   * Net velocity [my / s]:  %10.6f\n' %track.vnet)
            datafile.write('#   * Overall RunLength [my]: %10.6f\n' %track.overallRunLength)
            datafile.write('#   * Pos. Runlength [my]:    %10.6f\n' %track.posRunLength)
            datafile.write('#   * Neg. Runlength [my]:    %10.6f\n' %track.negRunLength)
            datafile.write('# --------------------------------------\n\n')
    print '... finished.'
    return 0
    
class particle_track:
    def __init__(self):
        self.dt = 1.299     # time between two points in seconds
        self.xpos = []      # x position in 1e-6 m
        self.ypos = []      # y position in 1e-6 m
    
    def compDataOfInterest(self):
        self.npos = len(self.xpos)
        # compute velocities per path segment
        self.v = []
        self.ds = []
        self.dist2refline = []
        for i in range(1, self.npos):
            dx = np.abs(self.xpos[i] - self.xpos[i-1])
            dy = np.abs(self.ypos[i] - self.ypos[i-1])
            self.ds.append(np.sqrt(dx** 2 + dy**2))
            self.v.append(self.ds[-1] / self.dt)
        self.vmean = np.mean(self.v)
        self.vmax = np.max(self.v)
        # compute net velocity
        dxnet = np.abs(self.xpos[-1] - self.xpos[0])
        dynet = np.abs(self.ypos[-1] - self.ypos[0])
        self.dsnet = np.sqrt(dxnet**2 + dynet**2)
        self.vnet = self.dsnet / (self.dt * self.npos)
        # compute Runlength
        self.posRunLength = []
        self.negRunLength = []
        self.overallRunLength = np.sum(self.ds)
        for i in range(self.npos):
            self.dist2refline.append(pointLineDist([self.xpos[i], self.ypos[i]], self.ref))
        posLegs = []
        negLegs = []
        for i in range(self.npos-1):
            if self.dist2refline[i+1] >= self.dist2refline[i]:
                posLegs.append(self.ds[i])
            else:
                negLegs.append(self.ds[i])
        self.posRunLength = np.sum(posLegs)
        self.negRunLength = np.sum(negLegs)
        return 0

if __name__ == '__main__':
    for directory in dirlist:
        os.chdir('./' + directory)
        print 'I am in %s' %os.getcwd()
        minLength = 4
        saveTracks2File(getParticleTracks(minLength))
        os.chdir('../')