import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
nrods = 51

sel = np.loadtxt('../list_selected')

from mpi4py import MPI
comm=MPI.COMM_WORLD
numproc=comm.size
rank = comm.Get_rank()

##### Curvature loocking for the radio of cincurscribed circle ############
def get_curvature_angle(p1, p2, p3):
    p1, p2, p3 = np.array(p1), np.array(p2), np.array(p3)
    a = np.linalg.norm(p1 - p2)
    b = np.linalg.norm(p2 - p3)
    c = np.linalg.norm(p3 - p1)
    avec, bvec = p1 - p2, p3 - p2
    cosangle = np.dot(avec, bvec)/(a * b)
    sign = np.sign(avec[0]*bvec[1] - avec[1]*bvec[0])
    if (-0.999 < cosangle < 0.999):
        return sign *((a+b+c) * (-a+b+c) * (a-b+c) * (a+b-c))**(0.5) / (1000*a*b*c)
    else:
        return 0
########################################################################
for k, ind in enumerate(sel):#[23]:#
    if k%numproc == comm.rank:
        curv = np.zeros((12000, 49))
        body = np.loadtxt('./filtered_agents/%d/body.dat'%ind) ## first 100 seconds of simulation
        x = body[:, range(1, 154, 3)]
        xc = np.sum(x, axis = 1)/nrods
        y = body[:, range(2, 154, 3)]
        yc = np.sum(y, axis = 1)/nrods
        v = np.sqrt((np.diff(xc))**2 + (np.diff(yc))**2)/0.2/0.0022 #0.05 correspond to dt, 0.0022 is the expected velocity
        
        np.savetxt('./filtered_agents/%d/speed.dat'%ind, v)
        for t in range(12000):
            x = body[t, range(1, 154, 3)]
            y = body[t, range(2, 154, 3)]
            for s in range(49):
                p1, p2, p3 = (x[s], y[s]), (x[s+1], y[s+1]),(x[s+2], y[s+2])
                curv[t, s] = get_curvature_angle(p1, p2, p3)
        np.savetxt('./filtered_agents/%d/curv.dat'%ind, curv)

