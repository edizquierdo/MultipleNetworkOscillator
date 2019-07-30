import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os

######################################################################
#####                        Simulation                    ###########
######################################################################
selected = np.loadtxt('list_selected')
######################################################################
#oscillations = np.zeros((15, 6))
#os.chdir('script')
#for k, sel in enumerate(selected):
#    o = np.zeros((9*9*9, 6))
#    l=0
#    os.system('cp ../selected/%d/best.gen.dat ./'%(sel))
#    for ASext in np.linspace(-5, 5, 9):
#        for DAext in np.linspace(-5, 5, 9):
#            for DBext in np.linspace(-5, 5, 9):
#                os.system('./main %.5f %.5f %.5f'%(ASext, DAext, DBext))
#                o[l] = np.loadtxt('osc.dat')
#                l=l+1
#    oscillations[k] = np.max(o, axis = 0)
#    print(k)
#np.savetxt('../Oscillation.dat', oscillations)
#os.chdir('../')

oscillations = np.loadtxt('Oscillation.dat')
######################################################################
#####                        FIGURE STUFF                  ###########
######################################################################
labels = ['AS'+ r'$\to$' +'VD\nablated', 'Full\nUnit', 'DA'+ r'$\to$' +'DD\nablated']

colors = ["r", "c", "y"]

plt.close('all')
fig = plt.figure(figsize = [30, 10])
gm = gridspec.GridSpec(120, 225)

ax1 = plt.subplot(gm[5:, 10:35])

axnet1 = plt.subplot(gm[5:, 42:82])
ax2 = plt.subplot(gm[5:, 91:128])

axnet2 = plt.subplot(gm[5:, 135:175])
ax3 = plt.subplot(gm[5:, 184:221])

######################################################################
######################################################################
osc_solutions = [1,2,4,5,8,9,10,11,12,13,14]
non_osc_sol = [0, 3, 6, 7]
print('Non oscillating solutions:')
for non in non_osc_sol: print(selected[non])

B_osc = oscillations[:, [1,3,5]]
fz = 35

for i in osc_solutions:
    ax2.plot(B_osc[i,[0,1]], '--ok', markerfacecolor = 'w', markeredgewidth = 5, ms = 22)
    ax3.plot(B_osc[i,[0,2]], '--ok', markerfacecolor = 'w', markeredgewidth = 5, ms = 22)
    ax2.plot([1], B_osc[i,[1]], 'o', markerfacecolor = colors[1], markeredgewidth = 0, ms = 26)
    ax3.plot([1], B_osc[i,[2]], 'o', markerfacecolor = colors[2], markeredgewidth = 0, ms = 26)

import pandas as pd
import seaborn as sns

df = pd.DataFrame(columns=['osc', 'circuit'])
for i in osc_solutions:
    df.loc[i] = [B_osc[i, 0], labels[1]]

sns.set(style="whitegrid", palette="muted")
b = sns.swarmplot(x="circuit", y="osc", palette=["w", "c", "y"], data=df, edgecolor = 'k', linewidth = 5, ax = ax1, s= 22)
#b.set_xticklabels([])
#b.set_xlabel(labels[1],fontsize=fz)

df = pd.DataFrame(columns=['osc', 'circuit'])
for i in non_osc_sol:
    df.loc[i] = [B_osc[i, 0], labels[1]]

sns.set(style="whitegrid", palette="muted")
b = sns.swarmplot(x="circuit", y="osc", palette=["k", "c", "y"], data=df, edgecolor = 'k', linewidth = 5, ax = ax1, s= 22)
b.set_xticklabels([])
b.set_xlabel(labels[1],fontsize=fz)


for ax in [ax1, ax2, ax3]:
    ax.set_ylabel('VB oscillation', fontsize = fz)
    ax.set_ylim(-0.02, 1.05)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(('%.1f'%a for a in [0, 0.2, 0.4, 0.6, 0.8, 1.0]), fontsize = 25)

ax2.text(1.0, 0.67, '*', fontsize = 38, ha = 'center', va = 'center')

for k, ax in enumerate([ax2, ax3]):
    ax.set_xticks([0, 1])
    ax.set_xlim(-0.2, 1.2)
ax2.set_xticklabels([labels[i] for i in [1, 0]], fontsize = fz)
ax3.set_xticklabels([labels[i] for i in [1, 2]], fontsize = fz)

######################################################################
##### Networks ###########
######################################################################
neurons = ['AS', 'DA', 'DB', 'DD', 'VD', 'VB', 'VA']
[AS, DA, DB, DD, VD, VB, VA] = [0, 1, 2, 3, 4, 5, 6]
xloc = [5, 1, 9, 7.5, 2.5, 8.5, 1.5]
yloc = [18, 13, 13, 10, 7.5, 2, 1]
labels, positions = [], []
for i in range(7):
    labels.append('%s'%(neurons[i]))
    positions.append([xloc[i], yloc[i]])
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
    arrow_shape = (0.8, 2.1, 2.1)
    arrowprops=dict(arrowstyle='simple, tail_width=%.2f,head_width=%.2f,head_length=%.2f'%(arrow_shape), facecolor = color, edgecolor = '0.25', linewidth = 0.6, connectionstyle = shape)
    line = axe.annotate("",  xy= To_edge, xycoords='data', xytext = From_edge, textcoords = 'data', arrowprops=arrowprops, zorder = 2)
    return line

def Gapjunction(axe, n1, n2, shape, radio, color):
    From, To = positions[n1], positions[n2]
    alpha = get_alpha(From, To)
    To_edge     = (To[0] - radio * np.cos(alpha), To[1] - radio * np.sin(alpha))
    From_edge   = (From[0] + radio * np.cos(alpha), From[1] + radio * np.sin(alpha))
    arrowprops=dict(arrowstyle='|-|, widthA=1.08,widthB=1.08', linewidth = 8, edgecolor = '0.25', connectionstyle =  shape)
    axe.annotate("",  xy= To_edge, xycoords='data', xytext = From_edge, textcoords = 'data', arrowprops=arrowprops, zorder=2)
    arrowprops=dict(arrowstyle='|-|, widthA=1.04,widthB=1.04', linewidth = 6, edgecolor = color, connectionstyle =  shape)
    axe.annotate("",  xy= To_edge, xycoords='data', xytext = From_edge, textcoords = 'data', arrowprops=arrowprops, zorder=2)
##################################
#####   Plot Network   ###########
##################################
R = 1.2
Rs = 1.5
##### DRAW THE CIRCLES ###########
circles1 = ['c%i' for i in range(len(labels))]
circles2 = ['c%i' for i in range(len(labels))]
fcolor = ['0.95'] * 14
textcolor = ['0.05'] * 14
#for n in [0, 1, 2, 7, 8, 9, 10, 11, 12, 13]: textcolor[n] = '0.95'

for n in np.arange(len(labels)):
    circles1[n] = plt.Circle(positions[n], radius = R, fc = fcolor[n], edgecolor = 'k', lw = 1, zorder = 20)
    circles1[n].set_alpha(1.00 )
    circles2[n] = plt.Circle(positions[n], radius = R, fc = fcolor[n], edgecolor = 'k', lw = 1, zorder = 20)
    circles2[n].set_alpha(1.00 )
    for ax in [axnet1, axnet2]:
        ax.text(positions[n][0], positions[n][1], labels[n][-4:], color = textcolor[n], ha = 'center', va = 'center', fontsize = 32, zorder = 22)
    axnet1.add_patch(circles1[n])
    axnet2.add_patch(circles2[n])
#    axnet2.add_patch(circles[n])
##### DRAW CONNECTIONS ###########
ablated_fitness = [colors[1]] * 10
#ablated_fitness[1] = '1.0'
ablated_fitness[1] = 'red'
#ablated_fitness[6] = '#cc7a00'
l = connect(axnet1, AS, VD, "arc3,rad= 0.00", Rs, ablated_fitness[1])#l.set_linestyle('--')
connect(axnet1, AS, DA, "arc3,rad= 0.00", Rs, ablated_fitness[0])
connect(axnet1, DA, DB, "arc3,rad= 0.00", Rs, ablated_fitness[2])
connect(axnet1, DB, AS, "arc3,rad= 0.00", Rs, ablated_fitness[3])
connect(axnet1, VD, VA, "arc3,rad= 0.00", Rs, ablated_fitness[4])
connect(axnet1, VD, VB, "arc3,rad= 0.00", Rs, ablated_fitness[5])
connect(axnet1, DA, DD, "arc3,rad= 0.00", Rs, ablated_fitness[6])
connect(axnet1, VB, DD, "arc3,rad= 0.00", Rs, ablated_fitness[7])
connect(axnet1, VA, DD, "arc3,rad= 0.00", Rs, ablated_fitness[8])
Gapjunction(axnet1, VD, DD, "arc3,rad= 0.0", Rs, ablated_fitness[9])


ablated_fitness = [colors[2]] * 10
ablated_fitness[6] = 'red'
#ablated_fitness[6] = '1.0'
#ablated_fitness[1] = '#cc7a00'
connect(axnet2, DA, DD, "arc3,rad= 0.00", Rs, ablated_fitness[6])
connect(axnet2, AS, DA, "arc3,rad= 0.00", Rs, ablated_fitness[0])
connect(axnet2, AS, VD, "arc3,rad= 0.00", Rs, ablated_fitness[1])
connect(axnet2, DA, DB, "arc3,rad= 0.00", Rs, ablated_fitness[2])
connect(axnet2, DB, AS, "arc3,rad= 0.00", Rs, ablated_fitness[3])
connect(axnet2, VD, VA, "arc3,rad= 0.00", Rs, ablated_fitness[4])
connect(axnet2, VD, VB, "arc3,rad= 0.00", Rs, ablated_fitness[5])
connect(axnet2, VB, DD, "arc3,rad= 0.00", Rs, ablated_fitness[7])
connect(axnet2, VA, DD, "arc3,rad= 0.00", Rs, ablated_fitness[8])
Gapjunction(axnet2, VD, DD, "arc3,rad= 0.0", Rs, ablated_fitness[9])

for ax in [axnet1, axnet2]:
    ax.set_xlim(-0.5, 10.5)
    ax.set_ylim(-1.0 , 19.5)
    ax.axis('off')
    ax.set_aspect('equal')
####################################################################

####################################################################
fz = 88
plt.figtext(0.02, 0.94, 'A', ha = 'center', va = 'center', fontsize = fz)
plt.figtext(0.22, 0.94, 'B', ha = 'center', va = 'center', fontsize = fz)
plt.figtext(0.62, 0.94, 'C', ha = 'center', va = 'center', fontsize = fz)


fig.subplots_adjust(left=0.02, bottom=0.15, right=1.00, top=1.00)
plt.savefig('dv_propagation.png', dpi = 100)
plt.savefig('dv_propagation.pdf')
#######################################################################
#######################################################################
