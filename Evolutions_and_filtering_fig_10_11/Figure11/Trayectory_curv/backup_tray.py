import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
nrods = 51


def get_curvature_angle(p1, p2, p3):
    p1, p2, p3 = np.array(p1), np.array(p2), np.array(p3)
    a = np.linalg.norm(p1 - p2)
    b = np.linalg.norm(p2 - p3)
    c = np.linalg.norm(p3 - p1)
    
    a1, b1 = p1 - p2, p3 - p2
    cosangle = np.dot(a1, b1)/(a * b)
    sign = np.sign(a1[0]*b1[1] - a1[1]*b1[0])
    if (-1 < cosangle < 1):
        return [((a+b+c) * (-a+b+c) * (a-b+c) * (a+b-c))**(0.5) / (a*b*c), sign * (np.pi - np.arccos(cosangle))]
    else:
        return [0, 0]

dt = 40
r = 2* 40 *10**-3 * np.abs(np.sin(np.arccos((np.arange(51) - 25)/(25+0.2))))

plt.close('all')
fig = plt.figure(figsize = [9,8])
gm = gridspec.GridSpec(100, 190)
ax0 = plt.subplot(gm[:, :])


def plot_trajectory(ind):
    body = np.loadtxt('../behavior/selected/%s/body.dat'%ind)[:, :] ## first 100 seconds of simulation
    print(body.shape)
    x = body[:, range(1, 154, 3)]*1000
    xc = np.sum(x, axis = 1)/nrods
    y = body[:, range(2, 154, 3)]*1000
    yc = np.sum(y, axis = 1)/nrods
    a = body[:, range(3, 154, 3)]
    
    xcc, ycc = np.mean(xc), np.mean(yc)
    radii = 0
    for ev in range(12):
        p1, p2, p3 = (xc[ev * dt], yc[ev * dt]), (xc[(ev+1) * dt], yc[(ev+1) * dt]), (xc[(ev+2) * dt], yc[(ev+2) * dt])
        [R, ang] = get_curvature_angle(p1, p2, p3)
        radii += 1/(R*12) # 12 evaluation of radii
    print(radii)

    ax0.plot(xc, yc, color = color, label = label+'\nradius = %.2f mm'%radii)
    ax0.plot([xcc-radii, xcc + radii], [ycc, ycc], '--', color = color)
    ax0.plot([xcc, xcc], [ycc-radii, ycc + radii], '--', color = color)

    ax0.plot(x[wt], y[wt], 'k', lw = 1)
    ax0.plot(x[wt] - r/2.0*np.cos(a[wt]), y[wt] - r/2.0*np.sin(a[wt]), 'k', lw = 1)
    ax0.plot(x[wt] + r/2.0*np.cos(a[wt]), y[wt] + r/2.0*np.sin(a[wt]), 'k', lw = 1)
#    
    polyx = np.r_[x[wt] - r/2.0*np.cos(a[wt]), (x[wt] + r/2.0*np.cos(a[wt]))[::-1]]
    polyy = np.r_[y[wt] - r/2.0*np.sin(a[wt]), (y[wt] + r/2.0*np.sin(a[wt]))[::-1]]
    
    ax0.fill(polyx, polyy, color = color, alpha = 0.6)
    return radii

color = 'green'
label = 'Accepted'
wt = 300 #time of simulation to plot the worm
radii = plot_trajectory(12)

wt = 1280
color = 'red'
label = 'Rejected'
radii = plot_trajectory(10)

ax0.set_aspect('equal')
ax0.legend(fontsize = 16)
ax0.set_xlabel('X coordenate (mm)', fontsize = 16)
ax0.set_ylabel('Y coordenate (mm)', fontsize = 16)

fig.subplots_adjust(left=0.09, bottom=0.06, right=0.99, top=1.00)
plt.savefig('./plot_2.png')
