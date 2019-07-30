import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
nrods = 51
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 24 
mpl.rcParams['ytick.labelsize'] = 24 
sel = [4, 5, 9, 21, 23, 30, 37, 42, 45, 53, 63, 66, 68, 101, 103]
#################     PLOT   ################
############## neural activity  ################
plt.close('all')
fig = plt.figure(figsize = [34,12])
gm = gridspec.GridSpec(65, 150)
ax0 = [plt.subplot(gm[5 + (i%3)*20 : 20 + (i%3)*20, 65 + (i%2)*43: 105 + (i%2)*43]) for i in [0,4, 2, 3, 1, 5]] # body posture
ax1 = plt.subplot(gm[4:30, :58]) # curvature
ax2 = plt.subplot(gm[32:58, :58]) #velocity
################################################
################################################
###################### CURVATURE  ################
body = np.loadtxt('body_23')## first 50 seconds of simulation
curv = np.loadtxt('curv_23')
v = np.loadtxt('velocity_filtered')
x = body[:, range(1, 154, 3)]
y = body[:, range(2, 154, 3)]
xc = np.sum(x, axis = 1)/nrods
yc = np.sum(y, axis = 1)/nrods
a = body[:, range(3, 154, 3)]
################################################
############ Trayectory #######
###############################
def plot_worm(ax, k):
    dx, dy = x[k]-xc[k], y[k]-yc[k]
#    b = np.arctan((dx[0]-dx[-1])/(dy[0]-dy[-1]))
    b = 3*np.pi/2 + np.arctan((dx[0]-dx[-1])/(dy[0]-dy[-1]))
    tx = np.cos(b)*dx - np.sin(b)*dy
    ty = np.cos(b)*dy + np.sin(b)*dx
    ta = a[k]+b
    ax.plot(tx, ty, 'k', lw = 1)
    ax.plot(tx[0], ty[0], '.g', ms = 18)
    ax.plot(tx - r/2.0*np.cos(ta), ty - r/2.0*np.sin(ta), 'k', lw = 1)
    ax.plot(tx + r/2.0*np.cos(ta), ty + r/2.0*np.sin(ta), 'k', lw = 1)
    ax.set_xlim(-0.0005, 0.0005)
    ax.set_ylim(-0.00015, 0.00015)
    ax.text(-0.0004, 0.0001, 't = %.2f'%(k/20.), ha = 'center', va = 'center', fontsize = 26)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal')

r = 2* 40 *10**-6 * np.abs(np.sin(np.arccos((np.arange(51) - 25)/(25+0.2)))) #body radius, from BBC model

for k, t in enumerate(np.linspace(0, 40, 6)):
    print(t)
    plot_worm(ax0[k], int(t))

################################################
fzl = 30
############ Curvature  #######
###############################
imcurv = ax1.imshow(curv, cmap=plt.get_cmap('seismic'), aspect='auto', vmin = -10, vmax = 10, origin = 'lower')
ax1.set_xlim(0, 120)
ax1.set_xticks([0, 40, 80, 120])
ax1.set_xticklabels([])
ax1.set_yticks([1, 21])
ax1.set_yticklabels(['',''])
ax1.text(-1, 1, 'Head', va = 'center', ha = 'right', fontsize = fzl)
ax1.text(-1, 21, 'Tail', va = 'center', ha = 'right', fontsize = fzl)

############ Velocity #######
###############################
s = np.where(np.array(sel)==23)[0]
ax2.plot(np.linspace(0, 10, len(v[s][0])), 1000*v[s][0], 'k', linewidth = 3)
ax2.axhline(y = 0.22, linestyle =  '--', color = 'r')
ax2.set_ylim(0.05, 0.35)
ax2.set_xlim(0, 6)
ax2.set_yticks([0.1, 0.2, 0.3])
ax2.set_ylabel('Velocity (mm/s)', fontsize = fzl, labelpad = 24)
ax2.set_xlabel('Time (s)', fontsize = fzl, labelpad = 22)
###############################
fz = 70
plt.figtext(0.018, 0.95, 'A', fontsize = fz, ha = 'center', va = 'center')
plt.figtext(0.018, 0.46, 'B', fontsize = fz, ha = 'center', va = 'center')
plt.figtext(0.460, 0.95, 'C', fontsize = fz, ha = 'center', va = 'center')
###############################
fig.subplots_adjust(left=0.08, bottom=0.02, right=1.0, top=0.98)
plt.savefig('behavior.png', dpi = 100)
plt.savefig('behavior.pdf')

