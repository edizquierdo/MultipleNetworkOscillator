import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
from aux import *
################################################################## 
pal=["r", "c", "y", 'b', 'm']
## colors according to figure necc and suff
new_order = [2, 1, 0, 3, 4]
pal = [pal[i] for i in new_order]
##############################################################
plt.close('all')
fig = plt.figure(figsize = [40,22])
#axes
gm = gridspec.GridSpec(54, 290)
axnet1 = plt.subplot(gm[2:23, :90])
ax01 = plt.subplot(gm[25:36, 5:90])
ax11 = plt.subplot(gm[40:52, 5:90])

axnet2 = plt.subplot(gm[2:23, 100:190])
ax02 = plt.subplot(gm[25:36, 105:190])
ax12 = plt.subplot(gm[40:52, 105:190])

axnet3 = plt.subplot(gm[2:23, 200:290])
ax03 = plt.subplot(gm[25:36, 205:290])
ax13 = plt.subplot(gm[40:52, 205:290])
################################################################## 
wrap = lambda x: ((x + np.pi)%(2*np.pi ) - np.pi)
unwrap = lambda x: x if (x<180) else (x-360)
#########################   VB_DB  #############################################
ind = 101
ext_input = 0.0
############### Network ###########################################
phen = np.loadtxt('../entrainment/ind/%d/phenotype.dat'%ind)
phen = phen[[21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 39, 40, 41, 42, 43]]
neurons = [0, 1, 2, 4, 5, 6, 7, 8, 9, 11, 12, 13]
synapses = [0, 1, 2, 3, 4, 5, 14]
plot_segments(axnet1, pal[0], neurons, synapses, phen)
############### Activity minimal network###########################
body = np.loadtxt('../entrainment/ind/%d/body.dat'%ind)[:301, :] 
plot_worm(ax11, ax01, body)

################################################################## 
#########################   DA_AS  #############################################
ind = 103
ext_input = 0.54
############### Network ###########################################
phen = np.loadtxt('../entrainment/ind/%d/phenotype.dat'%ind)
phen = phen[[21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 39, 40, 41, 42, 43]]
neurons = [0, 1, 2, 4, 5, 6, 7, 8, 9, 11, 12, 13]
synapses = [0, 1, 2, 3, 4, 5, 13]
plot_segments(axnet2, pal[1], neurons, synapses, phen)
############### Activity minimal network###########################
body = np.loadtxt('../entrainment/ind/%d/body.dat'%ind)[:301, :] 
plot_worm(ax12, ax02, body)
################################################################## 
#########################   AS_VA  #############################################
ind = 30
ext_input = -0.35
############### Network ###########################################
phen = np.loadtxt('../entrainment/ind/%d/phenotype.dat'%ind)
phen = phen[[21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 39, 40, 41, 42, 43]]
neurons = [0, 1, 2, 4, 5, 6, 7, 8, 9, 11, 12, 13]
synapses = [0, 1, 2, 3, 4, 5, 12]
plot_segments(axnet3, pal[2], neurons, synapses, phen)
############### Activity minimal network###########################
body = np.loadtxt('../entrainment/ind/%d/body.dat'%ind)[:301, :] 
plot_worm(ax13, ax03, body)

################################################################## 
colors = ['#d95f02','#7570b3', '#1b9e77']
lstyle = 'simple, tail_width = %s, head_width = %s, head_length = %s'%(1.8, 4.1, 4.1)
yloc = -0.04
axnet1.annotate('', xy=(0.75, yloc), xycoords='axes fraction', xytext=(0.25, yloc),arrowprops=dict(arrowstyle=lstyle, color=colors[0]))
axnet2.annotate('', xy=(0.75, yloc), xycoords='axes fraction', xytext=(0.25, yloc),arrowprops=dict(arrowstyle=lstyle, color=colors[0]))
axnet3.annotate('', xy=(0.25, yloc), xycoords='axes fraction', xytext=(0.75, yloc),arrowprops=dict(arrowstyle=lstyle, color=colors[1]))

################################################################## 
################################################################## 
################     DATA     ##################
anterio30 = np.loadtxt('../entrainment/data/data_entrain_head_30')
anterio101 = np.loadtxt('../entrainment/data/data_entrain_head_101')
anterio103 = np.loadtxt('../entrainment/data/data_entrain_head_103')

posterio30 = np.loadtxt('../entrainment/data/data_entrain_tail_30')
posterio101 = np.loadtxt('../entrainment/data/data_entrain_tail_101')
posterio103 = np.loadtxt('../entrainment/data/data_entrain_tail_103')
a = [anterio30, anterio101, anterio103, posterio30, posterio101, posterio103]
for k, j in enumerate(a):
    for i in range(6): j[i] = unwrap(j[i])
    a[k] = a[k][[4, 5, 0, 1, 2, 3]]
#    a[k] = np.append(j, j[0])
[anterio30, anterio101, anterio103, posterio30, posterio101, posterio103] = a

    
################################################################## 
lw, ms = 6, 20
colors = ['#d95f02','#7570b3', '#1b9e77']
ax11.plot(anterio101, '--d', color = colors[0], lw = lw, ms = ms, markeredgecolor = 'w')
ax11.plot(posterio101, '--o', color = colors[1], lw = lw, ms = ms)
ax12.plot(anterio103, '--d', color = colors[0], lw = lw, ms = ms)
ax12.plot(posterio103, '--o', color = colors[1], lw = lw, ms = ms)
ax13.plot(anterio30, '--d', color = colors[0], lw = lw, ms = ms, label = 'Head to tail entrainment')
ax13.plot(posterio30, '--o', color = colors[1], lw = lw, ms = ms, label = 'Tail to head entrainment')

ind = [101, 103, 30]
fz = 25
for k, ax in enumerate([ax11, ax12, ax13]):
    ax.set_ylim(-182, 182)
    ax.set_ylabel('Head-Tail phase shift', fontsize = fz)
    ax.set_xticks([0, 1, 2, 3, 4, 5])
    ax.set_xticklabels(['-120', '-60', '0', '60', '120', '180'], fontsize = fz)
    ax.set_yticks([-180, -90, 0, 90, 180])
    ax.set_yticklabels(['-180', '-90', '0', '90', '180'], fontsize = fz)
    ax.set_xlabel('Degrees entrainment displacement', fontsize = fz)

#ax13.legend(fontsize = 24)

legend_neck = ax13.legend(bbox_to_anchor=(-0.15, 0.85, 0.6, .102), 
ncol=1, borderaxespad=0., fontsize = 24, handletextpad=0.0, markerscale = 1.1, borderpad = 0.0, labelspacing = 1.0, frameon=False)

################################################################## 
xpos = 0.015, 0.340, 0.675
ypos = [0.97, 0.6, 0.3]
fzl = 60

[plt.figtext(xpos[i], ypos[0], 'A%d'%(i+1), fontsize = fzl, ha = 'center', va = 'center') for i in range(3)]
[plt.figtext(xpos[i], ypos[1], 'B%d'%(i+1), fontsize = fzl, ha = 'center', va = 'center') for i in range(3)]
[plt.figtext(xpos[i], ypos[2], 'C%d'%(i+1), fontsize = fzl, ha = 'center', va = 'center') for i in range(3)]

fig.subplots_adjust(left=0.018, bottom= 0.00, right=0.99, top=1.06)
plt.savefig('./kinogram_strategies.png')
plt.savefig('./kinogram_strategies.pdf')

