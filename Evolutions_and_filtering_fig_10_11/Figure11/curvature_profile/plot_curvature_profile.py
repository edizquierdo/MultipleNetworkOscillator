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

'''
curv_profile_rad = np.zeros((nsel, n_seg))
curv_profile_ang = np.zeros((nsel, n_seg))

for ind in range(nsel):
    curv_profile_rad[ind] = normalize(smooth(np.loadtxt('curvature_data/mean_curv_rad%d'%ind)))
    curv_profile_ang[ind]  = normalize(smooth(np.loadtxt('curvature_data/mean_curv_ang%d'%ind)))

np.savetxt('curv_profile_rad', curv_profile_rad)
np.savetxt('curv_profile_ang', curv_profile_ang)
'''
##################################################
############### Figure ###########################
curv_profile_rad = np.loadtxt('curv_profile_rad')

plt.close('all')

fig = plt.figure(figsize = [6,8])

colors = ['#2689ad','#d95f02','#7570b3']

gm = gridspec.GridSpec(40, 42)
ax0 = plt.subplot(gm[:19, 1:40])
ax10 = plt.subplot(gm[21:, 1:40])



numb_growing_in_rad, numb_growing_in_ang = 0,0

for ind in range(nsel):
    m , c = get_slope(curv_profile_rad[ind])
    if m > 0:
        ax0.plot(curv_profile_rad[ind], color = 'red')
        numb_growing_in_rad +=1
    else:
        ax10.plot(curv_profile_rad[ind], color = 'green')



for ax in [ax0, ax10]:
    ax.set_xlim(pos_norm, n_seg-1)
    ax.set_xticks([])
    ax.set_yticks([0.5, 1.0, 1.5, 2.0])
    ax.set_ylim(0.4, 2.1)
    ax.set_ylabel('Normalized curvature', fontsize = 15)


ax10.set_xticks([1+pos_norm, n_seg-2])
ax10.set_xticklabels(['Head', 'Tail'], fontsize = 15)




plt.figtext(0.5, 0.96, 'Negative slope N=%d'%(numb_growing_in_rad), color = 'k', fontsize = 20, ha = 'center')
plt.figtext(0.5, 0.485, 'Positive slope N=%d'%(nsel - numb_growing_in_rad), color = 'k', fontsize = 20, ha = 'center')

fig.subplots_adjust(left = 0.09, bottom = 0.05, right = 0.99, top = 0.95)

plt.savefig('plot_2.png', dpi = 200)
