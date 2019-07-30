import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
######################################################################
##### Simulation
 ###########
######################################################################
#oscillations = np.zeros((15, 3))
#selected = np.loadtxt('list_selected')
#os.chdir('script')
#for k, sel in enumerate(selected):
#    o = np.zeros((9*9*9, 3))
#    l=0
#    os.system('cp ../selected/%d/best.gen.dat ./'%(sel))
#    for ext1 in np.linspace(-5, 5, 9):
#        for ext2 in np.linspace(-5, 5, 9):
#            for ext3 in np.linspace(-5, 5, 9):
#                os.system('./main %.5f %.5f %.5f'%(ext1, ext2, ext3))
#                o[l] = np.loadtxt('osc.dat')
#                l=l+1
#    oscillations[k] = np.max(o, axis = 0)
#np.savetxt('../Oscillation.dat', oscillations)
#os.chdir('../')

oscillations = np.loadtxt('Oscillation.dat')

######################################################################
##### Plot ###########
######################################################################
plt.close('all')
fig = plt.figure(figsize = [32, 10])
gm = gridspec.GridSpec(45, 217)
axnet = plt.subplot(gm[1:, :128])
ax = plt.subplot(gm[1:-6, 140:-2])
######################################################################
neurons = ['AS', 'DA', 'DB', 'DD', 'VD', 'VB', 'VA']
[AS, DA, DB, DD, VD, VB, VA] = [0, 1, 2, 3, 4, 5, 6]
xloc = [5, 1, 9, 7.5, 2.5, 8.5, 1.5]
yloc = [18, 13, 13, 10, 7.5, 2, 1]
labels, positions = [], []
for seg in [0, 1, 2]:
    for i in range(7):
        labels.append('%s'%(neurons[i]))
        positions.append([13 * seg +2+ xloc[i], yloc[i]])
######################################################################
######################################################################
def get_alpha(From, To):
    if To[0] == From[0]: alpha = np.pi/4 - np.pi/4 * ( 1- np.sign(To[1] - From[1]))
    if (To[0] >= From[0]) & (To[1] >= From[1]) : alpha = np.arctan((To[1] - From[1])/(To[0] - From[0]))
    if (To[0] <= From[0]) & (To[1] <= From[1]) : alpha = np.pi + np.arctan((To[1] - From[1])/(To[0] - From[0]))
    if (To[0] < From[0]) & (To[1] > From[1]) : alpha = np.pi + np.arctan((To[1] - From[1])/(To[0] - From[0]))
    if (To[0] > From[0]) & (To[1] < From[1]) : alpha = np.arctan((To[1] - From[1])/(To[0] - From[0]))
    return (alpha)

def connect(axe, n1, n2, shape, radio, color):
    From, To = positions[n1], positions[n2]
    alpha = get_alpha(From, To)
    To_edge     = (To[0] - radio * np.cos(alpha), To[1] - radio * np.sin(alpha))
    From_edge   = (From[0] + radio * np.cos(alpha), From[1] + radio * np.sin(alpha))
    arrow_shape = (1.0, 2.5, 2.5)
    arrowprops=dict(arrowstyle='simple, tail_width=%.2f,head_width=%.2f,head_length=%.2f'%(arrow_shape), facecolor = color, edgecolor = 'k', linewidth = 0.6, connectionstyle = shape)
    axe.annotate("",  xy= To_edge, xycoords='data', xytext = From_edge, textcoords = 'data', arrowprops=arrowprops, zorder = 2)

def Gapjunction(axe, n1, n2, shape, radio, color):
    From, To = positions[n1], positions[n2]
    alpha = get_alpha(From, To)
    To_edge     = (To[0] - radio * np.cos(alpha), To[1] - radio * np.sin(alpha))
    From_edge   = (From[0] + radio * np.cos(alpha), From[1] + radio * np.sin(alpha))
    arrowprops=dict(arrowstyle='|-|, widthA=1.08,widthB=1.08', linewidth = 9, edgecolor = 'k', connectionstyle =  shape)
    axe.annotate("",  xy= To_edge, xycoords='data', xytext = From_edge, textcoords = 'data', arrowprops=arrowprops, zorder=2)
    arrowprops=dict(arrowstyle='|-|, widthA=1.04,widthB=1.04', linewidth = 7, edgecolor = color, connectionstyle =  shape)
    axe.annotate("",  xy= To_edge, xycoords='data', xytext = From_edge, textcoords = 'data', arrowprops=arrowprops, zorder=2)
##################################
#####   Plot Network   ###########
##################################
R = 1.2
Rs = 1.5
colors = ["r", "c", "y"]
##### DRAW THE CIRCLES ###########
circles = ['c%i' for i in range(len(labels))]
fcolor = ['0.95'] * 21
for n in [0, 1, 2]: fcolor[n] = colors[0]
for n in [10, 11, 13]: fcolor[n] = colors[1]
for n in [17, 18, 19]: fcolor[n] = colors[2]
textcolor = ['0.01'] * 21
for n in [0, 1, 2, 10, 11, 13, 17, 18, 19]: textcolor[n] = '0.90'

for n in np.arange(len(labels)):
    circles[n] = plt.Circle(positions[n], radius = R, fc = fcolor[n], edgecolor = 'k', lw = 1, zorder = 20)
    circles[n].set_alpha(1.00 )
    if n in [0, 1, 2, 10, 11, 13, 17, 18, 19]:
        axnet.text(positions[n][0], positions[n][1], labels[n][-4:], color = textcolor[n], ha = 'center', va = 'center', fontsize = 36, weight = 'bold', zorder = 22)
    else:
        axnet.text(positions[n][0], positions[n][1], labels[n][-4:], color = textcolor[n], ha = 'center', va = 'center', fontsize = 36, zorder = 22)
    
    axnet.add_patch(circles[n])
##### DRAW CONNECTIONS ###########
ablated_fitness = ['0.99'] * 10
for n in [0, 2, 3]: ablated_fitness[n] = colors[0]
seg = 0
connect(axnet, AS+seg, DA+seg, "arc3,rad= 0.00", Rs, ablated_fitness[0])
connect(axnet, AS+seg, VD+seg, "arc3,rad= 0.00", Rs, ablated_fitness[1])
connect(axnet, DA+seg, DB+seg, "arc3,rad= 0.00", Rs, ablated_fitness[2])
connect(axnet, DB+seg, AS+seg, "arc3,rad= 0.00", Rs, ablated_fitness[3])
connect(axnet, VD+seg, VA+seg, "arc3,rad= 0.00", Rs, ablated_fitness[4])
connect(axnet, VD+seg, VB+seg, "arc3,rad= 0.00", Rs, ablated_fitness[5])
connect(axnet, DA+seg, DD+seg, "arc3,rad= 0.00", Rs, ablated_fitness[6])
connect(axnet, VB+seg, DD+seg, "arc3,rad= 0.00", Rs, ablated_fitness[7])
connect(axnet, VA+seg, DD+seg, "arc3,rad= 0.00", Rs, ablated_fitness[8])
Gapjunction(axnet, VD+seg, DD+seg, "arc3,rad= 0.0", Rs, ablated_fitness[9])

ablated_fitness = ['0.99'] * 10
for n in [4, 8, 9]: ablated_fitness[n] = colors[1]
seg = 7
connect(axnet, AS+seg, DA+seg, "arc3,rad= 0.00", Rs, ablated_fitness[0])
connect(axnet, AS+seg, VD+seg, "arc3,rad= 0.00", Rs, ablated_fitness[1])
connect(axnet, DA+seg, DB+seg, "arc3,rad= 0.00", Rs, ablated_fitness[2])
connect(axnet, DB+seg, AS+seg, "arc3,rad= 0.00", Rs, ablated_fitness[3])
connect(axnet, VD+seg, VA+seg, "arc3,rad= 0.00", Rs, ablated_fitness[4])
connect(axnet, VD+seg, VB+seg, "arc3,rad= 0.00", Rs, ablated_fitness[5])
connect(axnet, DA+seg, DD+seg, "arc3,rad= 0.00", Rs, ablated_fitness[6])
connect(axnet, VB+seg, DD+seg, "arc3,rad= 0.00", Rs, ablated_fitness[7])
connect(axnet, VA+seg, DD+seg, "arc3,rad= 0.00", Rs, ablated_fitness[8])
Gapjunction(axnet, VD+seg, DD+seg, "arc3,rad= 0.0", Rs, ablated_fitness[9])

ablated_fitness = ['0.99'] * 10
for n in [5, 7, 9]: ablated_fitness[n] = colors[2]
seg = 14
connect(axnet, AS+seg, DA+seg, "arc3,rad= 0.00", Rs, ablated_fitness[0])
connect(axnet, AS+seg, VD+seg, "arc3,rad= 0.00", Rs, ablated_fitness[1])
connect(axnet, DA+seg, DB+seg, "arc3,rad= 0.00", Rs, ablated_fitness[2])
connect(axnet, DB+seg, AS+seg, "arc3,rad= 0.00", Rs, ablated_fitness[3])
connect(axnet, VD+seg, VA+seg, "arc3,rad= 0.00", Rs, ablated_fitness[4])
connect(axnet, DA+seg, DD+seg, "arc3,rad= 0.00", Rs, ablated_fitness[6])
connect(axnet, VA+seg, DD+seg, "arc3,rad= 0.00", Rs, ablated_fitness[8])
connect(axnet, VD+seg, VB+seg, "arc3,rad= 0.00", Rs, ablated_fitness[5])
connect(axnet, VB+seg, DD+seg, "arc3,rad= 0.00", Rs, ablated_fitness[7])
Gapjunction(axnet, VD+seg, DD+seg, "arc3,rad= 0.0", Rs, ablated_fitness[9])

axnet.set_xlim(0.5, 40.5)
axnet.set_ylim(-1.0 , 19.5)
axnet.axis('off')
axnet.set_aspect('equal')
####################################################################
import pandas as pd
import seaborn as sns

df = pd.DataFrame(columns=['osc', 'circ'])
for i in range(15):
    df.loc[i] = [oscillations[i, 0], 'AS-DA-DB']
    df.loc[i+15] = [oscillations[i, 1], 'VD-VA-DD']
    df.loc[i+30] = [oscillations[i, 2], 'VD-VB-DD']

import seaborn as sns
sns.set(style="whitegrid", palette="muted")
#sns.set(font_scale=7)#sns.set_context("poster")
b = sns.swarmplot(x="circ", y="osc", palette=["r", "c", "y"], data=df, ax = ax, s= 22)
ax.set_ylabel('Oscillatory activity',size=35)
ax.set_xlabel('Recurrent subcircuit',size=35)
b.tick_params(labelsize=25)

#ax.axis('equal')
####################################################################
plt.figtext(0.015, 0.94, 'A', ha = 'center', va = 'center', fontsize = 88)
plt.figtext(0.59, 0.94, 'B', ha = 'center', va = 'center', fontsize = 88)

fig.subplots_adjust(left=0.00, bottom=0.0, right=1.0, top=1.00)
plt.savefig('origin_oscillation.png', dpi = 100)
plt.savefig('origin_oscillation.pdf')
#plt.savefig('../../Figures/origin_oscillation.pdf')
#plt.savefig('../../Figures/origin_oscillation.png', dpi = 300)
