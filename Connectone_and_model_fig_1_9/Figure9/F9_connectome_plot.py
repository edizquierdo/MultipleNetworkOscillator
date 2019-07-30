import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
from draw_H_segment import *

object_fontsize = 45

######################################################################
#### Figure and axes definition
######################################################################

fig = plt.figure(figsize = [60, 40])
R, Rr, Rgj = 1.15, 1.20, 1.29
gmain = gridspec.GridSpec(200, 20)
ax = plt.subplot(gmain[0:95, :])
axmale = plt.subplot(gmain[105:200, :])
######################################################################
pal=["y", "r", "c"]

'''######################################################################
## HERMAPHRODYTE
######################################################################'''
#VB-DB
nvb = [3, 4, 6, 4, 1, 2, 4, 4]
ndb = [1, 2, 3, 1, 2, 4, 4, 7]
label = ['J,V,C', 'J,C', 'J,C', 'V', 'V,C', 'V', 'V', 'C']
displ = [0.8, 0.2, 0.15, 0.7, 0.25, 0.8, 0.5, 0.8]

for i in range(len(label)):
    n1, n2 = labels.index('VB%.2d'%nvb[i]), labels.index('DB%.2d'%ndb[i])
    Gapjunction(ax, n1, n2, Rgj, pal[0], displ[i], label[i])

#AS-VA
nas = [1, 1, 2, 2, 10, 11, 3, 5, 6, 10]
nva = [3, 6, 3, 4, 12, 11, 5, 7, 8, 11]
label = ['J,V,C', 'J', 'J,C', 'J,V,C', 'J', 'J', 'V', 'V', 'V', 'C']
displ = [0.15, 0.1, 0.4, 0.25, 0.25, 0.23, 0.3, 0.2, 0.2, 0.25]

for i in range(len(label)):
    n1, n2 = labels.index('AS%.2d'%nas[i]), labels.index('VA%.2d'%nva[i])
    Gapjunction(ax, n1, n2, Rgj, pal[1], displ[i], label[i])

#DA-AS
nda = [2, 2, 3, 5, 7, 8, 9, 2]
nas = [2, 4, 4, 6, 10, 11, 10, 3]
label = ['J,V,C', 'J', 'J,V,C', 'J,V,C', 'J', 'J,C', 'J', 'V,C']
displ = [0.5, 0.3, 0.4, 0.5, 0.4, 0.5, 0.3, 0.5]

for i in range(len(label)):
    n1, n2 = labels.index('DA%.2d'%nda[i]), labels.index('AS%.2d'%nas[i])
    Gapjunction(ax, n1, n2, Rgj, pal[2], displ[i], label[i])

##################################
##### DRAW THE CIRCLES ###########
##################################
circles = ['c%i' for i in range(len(labels))]
for n in np.arange(len(labels)):
    circles[n] = plt.Circle(positions[n], radius = R, fc = colors[n], edgecolor = 'k', lw = 1)
    circles[n].set_alpha(alpha_foreground + 0.15)
    ax.text(positions[n][0], positions[n][1], labels[n][-2:], color = 'w', ha = 'center', va = 'center', fontsize = object_fontsize, weight='bold')
    ax.add_patch(circles[n])

'''######################################################################
## MALE
######################################################################'''
#VB-DB
nvb = [3, 4, 5, 5, 7, 7, 8, 8, 8, 9, 10, 10, 10, 11, 11]
ndb = [1, 2, 3, 5, 3, 5, 3, 5, 7, 4, 4, 5, 6, 6, 7]
label = ['J', 'J', 'J,C', 'J,C', 'J,C', 'J,C', 'J,C', 'J,C', 'J,C', 'J,C', 'J,C', 'J,C', 'J,C', 'J,C', 'J,C']
displ = [0.15, 0.2, 0.2, 0.2, 0.2, 0.2, 0.15, 0.75, 0.12, 0.1, 0.2, 0.2, 0.2, 0.15, 0.2]

for i in range(len(label)):
    n1, n2 = labels.index('VB%.2d'%nvb[i]), labels.index('DB%.2d'%ndb[i])
    Gapjunction(axmale, n1, n2, Rgj, pal[0], displ[i], label[i])

#AS-VA
nas = [1, 1, 2, 2, 8, 10]
nva = [3, 6, 3, 4, 12, 12]
label = ['J', 'J', 'J', 'J', 'J,C', 'J,C']
displ = [0.2, 0.12, 0.32, 0.2, 0.25, 0.2]

for i in range(len(label)):
    n1, n2 = labels.index('AS%.2d'%nas[i]), labels.index('VA%.2d'%nva[i])
    Gapjunction(axmale, n1, n2, Rgj, pal[1], displ[i], label[i])

#DA-AS
nda = [2, 2, 3, 9]
nas = [2, 4, 4, 11]
label = ['J', 'J', 'J', 'J,C']
displ = [0.5, 0.3, 0.4, 0.5]

for i in range(len(label)):
    n1, n2 = labels.index('DA%.2d'%nda[i]), labels.index('AS%.2d'%nas[i])
    Gapjunction(axmale, n1, n2, Rgj, pal[2], displ[i], label[i])

###################################
###### DRAW THE CIRCLES ###########
##################################
circles = ['c%i' for i in range(len(labels))]
for n in np.arange(len(labels)):
    circles[n] = plt.Circle(positions[n], radius = R, fc = colors[n], edgecolor = 'k', lw = 1)
    circles[n].set_alpha(alpha_foreground + 0.15)
    axmale.text(positions[n][0], positions[n][1], labels[n][-2:], color = 'w', ha = 'center', va = 'center', fontsize = object_fontsize, weight='bold')
    axmale.add_patch(circles[n])
#####################
#####################
#### Left Labels  and axes formats ##
#####################
lab_y = [18, 16, 13, 10, 7, 4, 1]
lab_y = [18, 15, 12, 9, 6, 3, 0]
lab_x, lab_xf = -4.5, 62.5
lab =['AS', 'DA', 'DB', 'DD', 'VD', 'VB', 'VA']
for i in range(7):
    ax.text(lab_x, lab_y[i], lab[i], fontsize = 60, va = 'center')
    ax.text(lab_xf, lab_y[i], lab[i], fontsize = 60, va = 'center')
    axmale.text(lab_x, lab_y[i], lab[i], fontsize = 60, va = 'center')
    axmale.text(lab_xf, lab_y[i], lab[i], fontsize = 60, va = 'center')

for axis in [ax, axmale]:
    axis.set_xlim(-6, 64)
    axis.set_ylim(-2 , 20)
    axis.axis('off')
    axis.set_aspect('equal')

#'''####################
##### save figure   ####
#####################'''
plt.figtext(0.01, 0.97, 'A', fontsize = 110)
plt.figtext(0.01, 0.47, 'B', fontsize = 110)
fig.subplots_adjust(left=-0.03, bottom=0.02, right=1.03, top=0.98)
plt.savefig('connections.png')
plt.savefig('connections.pdf')
plt.clf()
