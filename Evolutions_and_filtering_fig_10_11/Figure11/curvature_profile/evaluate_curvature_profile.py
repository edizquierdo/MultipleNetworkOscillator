import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
nrods = 51
nsel = 104

### For each of the 104 selected worms we have a 49 poins vector that have the mean curvature calculated
### from 1000 points in time, that correspond to 50 seconds of simulation (22 cycles aprox at a freq = 0.44Hz)

### To smooth the curvature vector (length 49) we are going to convolve it by a vector of ones of lenght 10 (this give us a vector of 40 points)
### Then we can normalize the curvature vector by it's first second third etc point

pos_norm = 7
n_conv = 9
n_seg = 49 + 1 - n_conv

def normalize(a):
    return a/a[pos_norm]
    
def smooth(x):
    return 1.0* np.convolve(x, np.ones(n_conv), 'valid')/n_conv
    
def get_slope(x):
    x = x[pos_norm:]
    A = np.vstack([np.arange(len(x)), np.ones(len(x))]).T
    m, c = np.linalg.lstsq(A, x)[0]
    return [m, c]

#'''
curv_profile_rad = np.zeros((nsel, n_seg))
curv_profile_ang = np.zeros((nsel, n_seg))

for ind in range(nsel):
    curv_profile_rad[ind] = normalize(smooth(np.loadtxt('curvature_data/mean_curv_rad%d'%ind)))
    curv_profile_ang[ind]  = normalize(smooth(np.loadtxt('curvature_data/mean_curv_ang%d'%ind)))

np.savetxt('curv_profile_rad', curv_profile_rad)
np.savetxt('curv_profile_ang', curv_profile_ang)

#'''
##################################################
############### Figure ###########################
curv_profile_rad = np.loadtxt('curv_profile_rad')
curv_profile_ang = np.loadtxt('curv_profile_ang')
gain = np.loadtxt('gain')

plt.close('all')

fig = plt.figure(figsize = [12,8])

colors = ['#2689ad','#d95f02','#7570b3']

gm = gridspec.GridSpec(40, 98)
ax0 = plt.subplot(gm[2:19, 1:40])
ax1 = plt.subplot(gm[2:19, 50:90])

ax10 = plt.subplot(gm[22:-1, 1:40])
ax11 = plt.subplot(gm[22:-1, 50:90])

axbar = plt.subplot(gm[2:-1, 96:98])


numb_growing_in_rad, numb_growing_in_ang = 0,0

for ind in range(nsel):
    param = 0.5 * (1 + gain[ind])
    color = [2*np.abs(param-0.5), 1 - param, param]
    m , c = get_slope(curv_profile_rad[ind])

    if m > 0:
        ax0.plot(curv_profile_rad[ind], color = color)
        numb_growing_in_rad +=1
    else:
        ax10.plot(curv_profile_rad[ind], color = color)
    
    m , c = get_slope(curv_profile_ang[ind])
    if m > 0:
        ax1.plot(curv_profile_ang[ind], color = color)
        numb_growing_in_ang +=1
    else:
        ax11.plot(curv_profile_ang[ind], color = color)
    
#    color = [1.0*ind/nsel, 1 - 1.0*ind/nsel, 1.0*ind/nsel]
    color = [ 2*np.abs(1.0*ind/nsel-0.5), 1 - 1.0*ind/nsel, 1.0*ind/nsel]
    axbar.axhline(y = ind, linewidth = 4.7, color = color)
    axbar.plot([0], [nsel], '.', color = color)

for ax in [ax0, ax1]:
    ax.set_xlim(pos_norm, n_seg-1)
    ax.set_xticks([])
    ax.set_ylim(0.5, 2.5)

for ax in [ax10, ax11]:
    ax.set_xlim(pos_norm, n_seg-1)
    ax.set_xticks([1+pos_norm, n_seg-2])
    ax.set_ylim(0.0, 1.5)
    ax.set_xticklabels(['Head', 'Tail'], fontsize = 15)

ax = axbar
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xticks([])
ax.set_ylim(0, 103)
ax.set_yticks([0, 52, 103])
ax.set_yticklabels(['0.2', '0.6', '1.0'])


plt.figtext(0.25, 0.88, 'Radious Curvature\n\nPositive slope N=%d'%numb_growing_in_rad, color = 'k', fontsize = 20, ha = 'center')
plt.figtext(0.7, 0.88, 'Angle Curvature\n\nPositive slope N=%d'%numb_growing_in_ang, color = 'k', fontsize = 20, ha = 'center')

plt.figtext(0.25, 0.44, 'Negative slope N=%d'%(nsel - numb_growing_in_rad), color = 'k', fontsize = 20, ha = 'center')
plt.figtext(0.7, 0.44, 'Negative slope N=%d'%(nsel - numb_growing_in_ang), color = 'k', fontsize = 20, ha = 'center')

fig.subplots_adjust(left = 0.05, bottom = 0.03, right = 0.96, top = 0.91)

plt.savefig('curvature_profile_norm.png', dpi = 200)
