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

dt = 80
curvature = np.zeros(104)
worm_axis = (np.arange(51)-25)/25.0
worm_shape = 11 * (1 - worm_axis**2)

## Plot worm on times 1, 5, 9 seconds

'''
for ind in range(104):#[0]:#
    plt.close('all')

    fig = plt.figure(figsize = [12,12])

    gm = gridspec.GridSpec(100, 190)
    ax0 = plt.subplot(gm[:, :])
    
    body = np.loadtxt('../behavior/selected/%s/body.dat'%ind)[:, :] ## first 100 seconds of simulation
    x = body[:, range(1, 154, 3)]
    xc = np.sum(x, axis = 1)/nrods
    y = body[:, range(2, 154, 3)]
    yc = np.sum(y, axis = 1)/nrods
    
    xcc, ycc = np.mean(xc), np.mean(yc)
    radii = 0
    
    for ev in range(12):
        p1, p2, p3 = (xc[ev * dt], yc[ev * dt]), (xc[(ev+1) * dt], yc[(ev+1) * dt]), (xc[(ev+2) * dt], yc[(ev+2) * dt])
        [R, ang] = get_curvature_angle(p1, p2, p3)
        radii += 1/(R*12) # 12 evaluation of radii

    curvature[ind] = radii
    print(radii)

    ax0.plot(xc, yc)
    ax0.plot([xcc-radii, xcc + radii], [ycc, ycc])
    ax0.plot([xcc, xcc], [ycc-radii, ycc + radii])

    for i in range(51):
        ax0.plot(x[-1][i], y[-1][i], '.k', ms = 6 + worm_shape[i])

    ax0.set_aspect('equal')
    plt.figtext(0.5, 0.9, 'Curvature = %.2e, Radius = %.2e'%(1./radii, radii), ha = 'center', fontsize = 22)

    fig.subplots_adjust(left=0.05, bottom=0.04, right=0.98, top=0.96)
    plt.savefig('trayectory_figures/plot_sel_%d.png'%ind)
#
#np.savetxt('trayectory_curvature', curvature)
'''
curvature = np.loadtxt('trayectory_curvature')
plt.close('all')

fig = plt.figure(figsize = [12,6])
    
gm = gridspec.GridSpec(100, 190)
ax2 = plt.subplot(gm[:, :])

h = ax2.hist(1./curvature, bins=20, range = [0, 2100], orientation="vertical", color = 'b')
ax2.set_xlim(0,2100)
#ax2.yaxis.tick_right()
#ax2.yaxis.set_label_position("right")
ax2.set_ylabel('Counts', fontsize = 18, rotation = 270, labelpad = 13)
ax2.set_ylim(0, 15)
ax2.set_xlabel('Curvature (1/m)', fontsize = 18, labelpad = 2)
#ax0.set_xticks([])
#ax2.set_yticks([25, 50, 75, 100])
#ax2.get_children()[19].set_facecolor('r')

plt.figtext(0.5, 0.94, 'Trayectory curvature distribution', fontsize = 20, va = 'center', ha = 'center')

plt.savefig('trayectory_curvature_distribution.png')
