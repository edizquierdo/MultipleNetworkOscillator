import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
nrods = 51
nsel = 104

##### Curvature loocking for the radio of cincurscribed circle
def get_curvature_angle(p1, p2, p3):
    p1, p2, p3 = np.array(p1), np.array(p2), np.array(p3)
    a = np.linalg.norm(p1 - p2)
    b = np.linalg.norm(p2 - p3)
    c = np.linalg.norm(p3 - p1)
    avec, bvec = p1 - p2, p3 - p2
    cosangle = np.dot(avec, bvec)/(a * b)
    sign = np.sign(avec[0]*bvec[1] - avec[1]*bvec[0])
    if (-1 < cosangle < 1):
        return [sign *((a+b+c) * (-a+b+c) * (a-b+c) * (a+b-c))**(0.5) / (a*b*c), sign * (np.pi - np.arccos(cosangle))]
    else:
        return [0, 0]


import time
ti = time.time()

for ind in range(104):
    curv_rad, curv_ang = np.zeros((1000, 49)), np.zeros((1000, 49))
    body = np.loadtxt('../../behavior/selected/%s/body.dat'%ind)[200:, :] ## first 100 seconds of simulation
    for t in range(1000):
        x = body[t, range(1, 154, 3)]
        y = body[t, range(2, 154, 3)]
        for s in range(49):
            p1, p2, p3 = (x[s], y[s]), (x[s+1], y[s+1]),(x[s+2], y[s+2])
            [curv_rad[t, s], curv_ang[t, s]] = get_curvature_angle(p1, p2, p3)
    np.savetxt('curvature_data/mean_curv_rad%d'%ind, np.mean(np.abs(curv_rad), axis = 0), fmt = '%.5e')
    np.savetxt('curvature_data/mean_curv_ang%d'%ind, np.mean(np.abs(curv_ang), axis = 0), fmt = '%.5e')
    print(ind, time.time() - ti)


