import numpy as np
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import os
#import pandas as pd
#import seaborn as sns

selected = [4, 5, 9, 21, 23, 30, 37, 42, 45, 53, 63, 66, 68, 101, 103] ## correlative numbering
#ordered = [5, 9, 30, 42, 45, 53, 4, 21, 23, 37, 63, 66, 68, 101, 103] ## group by direcctionality 
ordered = [23, 66, 68, 101, 21, 42, 103, 9, 30, 4, 53, 63, 5, 37, 45] ## group by necc & suff 
ordered = [5, 12, 13, 14, 4, 8, 15, 3, 6, 1, 10, 11, 2, 7, 9] ## group by necc & suff 
index = np.argsort(ordered)
index = index*1.0
##################################################################
##################################################################
labels_syn = ['AS'+ r'$\leftrightsquigarrow$' +'VA' + r'$^{+1}$'\
, 'DA'+ r'$\leftrightsquigarrow$' +'AS' + r'$^{+1}$'\
, 'VB'+ r'$\leftrightsquigarrow$' +'DB' + r'$^{+1}$'\
, 'DB'+ r'$\to$' +'DD' + r'$^{+1}$'\
, 'VA' + r'$^{+1}$'+ r'$\to$' +'DD']
pal=["r", "c", "y", 'b', 'm']
##################################################################
necc = np.loadtxt('neccesity/data_ablation.dat')
suff = np.loadtxt('sufficienty/data_suff.dat')
necc[np.where(necc<0)] = 0 
necc[np.where(necc>1)] = 1 
suff[np.where(suff<0)] = 0 
suff[np.where(suff>1)] = 1 
################################################################## 
## re-sorting connections as to make them appear like in the text
new_order = [2, 1, 0, 3, 4]
necc = necc[:,new_order]
suff = suff[:,new_order]
labels_syn = [labels_syn[i] for i in new_order]
pal = [pal[i] for i in new_order]
##################################################################
plt.close('all')
fig = plt.figure(figsize = [16, 5])
gm = gridspec.GridSpec(150, 145)
ax = plt.subplot(gm[7:-19, 0:140])
##########################
s = 0.14
for i in range(15):
    if (index[i]>3):
        index[i] += 0.1
    if (index[i]>7):
        index[i] += 0.1
    if (index[i]>9):
        index[i] += 0.25
    if (index[i]>12):
        index[i] += 0.1
for i in range(15):
    ax.axvspan(index[i],index[i]+0.94, color = '0.91', lw = 0)
##########################
for i in range(5):
    ax.plot(index+s*i+0.1, 10-necc[:,i], 's', color = pal[i], ms = 8, label = labels_syn[i]) ## only for labels purposes
    ax.plot(index+s*i+0.2, 1-necc[:,i], 'o', markerfacecolor = pal[i], markeredgewidth = 2, color = pal[i], ms = 10)
    ax.plot(index+s*i+0.2, suff[:,i], 'o', markerfacecolor = 'w', markeredgewidth = 2,color = pal[i], ms = 10)
ax.set_xticks(index+0.5)

############## Connect important points for explanation ############
lw = 11.5
for i in [4, 11, 12, 13]: ## Connecting points in solution (a)
    j = 0
    plt.plot([(index[i]+s*j+0.2), (index[i]+s*j+0.2)], [1-necc[i, j], suff[i, j]], lw = lw, color = pal[int(j)], alpha = 0.25, zorder = 1)
for i in [3, 7, 14]: ## Connecting points in solution (b)
    j = 1
    plt.plot([(index[i]+s*j+0.2), (index[i]+s*j+0.2)], [1-necc[i, j], suff[i, j]], '-', lw = lw, color = pal[int(j)], alpha = 0.25, zorder = 1)
for i in [2, 5]: ## Connecting points in solution (c)
    j = 2
    plt.plot([(index[i]+s*j+0.2), (index[i]+s*j+0.2)], [1-necc[i, j], suff[i, j]], '-', lw = lw, color = pal[int(j)], alpha = 0.25, zorder = 1)


##########################
legend_connections = plt.legend(bbox_to_anchor=(1.001, 0.5, 0.1, .302), 
ncol=1, borderaxespad=0., fontsize = 13, handletextpad=0.0, markerscale = 1.5, borderpad = 0.0, labelspacing = 1.3, frameon=False)
plt.gca().add_artist(legend_connections)

j = [100, 1001]
l1, = ax.plot(j, j, 'ok', markerfacecolor = 'k', markeredgewidth = 2, ms = 7, label = 'Necessary')
l2, = ax.plot(j, j, 'ok', markerfacecolor = 'w', markeredgewidth = 2, ms = 7, label = 'Sufficient')
plt.legend(handles = [l1, l2], bbox_to_anchor=(1.001, 0.8, 0.1, .20), 
ncol=1, borderaxespad=0., fontsize = 13, handletextpad=0.05, markerscale = 1.4, borderpad = 0.05, labelspacing = 1.3, frameon=False)
##########################
xlabels = ['M%s'%str(a) for a in np.sort(ordered)]
for s in [13, 14, 5]: xlabels[s] = xlabels[s]+'*'
ax.set_xticklabels(xlabels, fontsize = 13)
ax.set_xlim(0, 15.49)
ax.set_ylim(-0.03, 1.03)
ax.set_ylabel('Locomotion performance', fontsize = 17)
ax.set_xlabel('Individual', fontsize = 17)
#####################################################################
compl = ['Simple', 'Redundant', 'Complex']
positions = [4.8, 10.95, 14.05]
[ax.text(positions[i], 1.16, compl[i], fontsize = 17, ha = 'center') for i in range(3)]

pos = np.array([0,9.3, 12.45, 15.5])/15.45
ypos = 1.11
arrowprops=dict(arrowstyle='|-|, widthA=0.2,widthB=0.2', linewidth = 1.25, edgecolor = 'k')
[ax.annotate("",  xy= [pos[i]+0.002, ypos], xycoords='axes fraction', xytext = [pos[i+1]-0.002,ypos], textcoords = 'axes fraction', arrowprops=arrowprops, zorder=0) for i in range(3)]


positions = [2, 5.6, 8.22]
[ax.text(positions[i], 1.07, labels_syn[i], fontsize = 14, ha = 'center') for i in range(3)]

pos = np.array([0.01, 4.00, 7.10, 9.29])/15.45
ypos = 1.015
arrowprops=dict(arrowstyle='|-|, widthA=0.15,widthB=0.15', linewidth = 1.2, edgecolor = 'k')
[ax.annotate("",  xy= [pos[i]+0.002, ypos], xycoords='axes fraction', xytext = [pos[i+1]-0.002,ypos], textcoords = 'axes fraction', arrowprops=arrowprops, zorder=0) for i in range(3)]

#####################################################################

fig.subplots_adjust(left=0.06, bottom=0.00, right=0.93, top=0.90)
plt.savefig('interunit.png', dpi = 100)
plt.savefig('interunit.pdf')


