import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
nrods = 51

sel = np.loadtxt('../list_selected')


#### Plot
plt.close('all')

fig = plt.figure(figsize = [7,7])
gm = gridspec.GridSpec(80, 105)
ax1 = plt.subplot(gm[2:58, :82])
axc = plt.subplot(gm[2:58, 86:90])
ax2 = plt.subplot(gm[60:77, :82])


curv_prof = np.loadtxt('curv_prof')
speeds_mean = np.loadtxt('speed')

m = ax1.imshow(curv_prof, aspect= 'auto', origin = 'lower', vmin = 0, vmax = 1)
#, interpolation = 'none'

barim = fig.colorbar(m, ax = ax1, cax = axc, ticks=[0, 0.5, 1])

axc.set_ylabel('Bending', labelpad = 18, rotation = -90, fontsize = 15)

ax2.errorbar(np.arange(8), np.mean(speeds_mean, axis = 0), np.std(speeds_mean, axis = 0), marker = '^')

for ax in [ax1, ax2]:
    ax.set_xlim(-0.5, 7.5)
    ax.set_xticks([0, 2, 4, 6])
    ax.set_xticklabels([])
ax2.set_xticklabels(['0', '0.5', '2', '8'])
ax2.set_xlabel('Gap junction overexpression in B-class motorneurons', fontsize = 15)

ax1.set_ylabel('Body\ncoordenates', labelpad = -10, fontsize = 15)
ax2.set_ylabel('Normalized\nSpeed', fontsize = 15)

ax1.set_yticks([4, 45])
ax1.set_yticklabels(['Head', 'Tail'], fontsize = 15)

ax2.set_ylim(-0.02, 1.05)
ax2.set_yticks([0, 0.5, 1.0])

fig.subplots_adjust(left=0.13, bottom=0.06, right=1.01, top=0.99)
plt.savefig('overexpression.png', dpi = 300)
plt.savefig('../../Figures/overexpression.png', dpi = 300)
#'''
