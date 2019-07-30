import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
nrods = 51
nsel = 104

speed = np.zeros(nsel)
for ind in range(nsel):
    ##################     CALCULATE   ################
    ###################### CURVATURE  ################
    body = np.loadtxt('selected/%s/body.dat'%ind)## 100 seconds of simulation
    x = body[:, range(1, 154, 3)]
    xc = np.sum(x, axis = 1)/nrods
    y = body[:, range(2, 154, 3)]
    yc = np.sum(y, axis = 1)/nrods
    v = np.sqrt((np.diff(xc))**2 + (np.diff(yc))**2)/0.05 #0.05 correspond to dt
    
    speed[ind] = np.mean(v)

np.savetxt('speed', speed)
