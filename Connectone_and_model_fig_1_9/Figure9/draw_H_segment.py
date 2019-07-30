import numpy as np
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
#####################################
########### correlative   ###########
#####################################

object_fontsize = 32
alpha_foreground = 0.72
[AS, DA, DB, DD,  VD, VB, VA] = [0, 1, 2, 3, 4, 5, 6]
neurons = ['AS', 'DA', 'DB', 'DD', 'VD', 'VB', 'VA']
n_classes = [11, 9, 7, 6, 13, 11, 12]
col = ['#cc66cc','#d8b365', '#5ab4ac', 'blue', '#99ccff', '#01665e', '#8c510a']
#######################################################
################ Circles locations ####################
xo = [1.1, -1, 2.9, 3, 0.9, 2, 0]
yloc = [18, 15, 12, 9, 6, 3, 0]
xrange = [80.2, 78, 80, 76, 82, 80, 78]
xrange = [60.2, 58, 60, 56, 61, 60, 58]

lab = [['%s%.2d'%(neurons[i], numb) for numb in range(1, 1+n_classes[i])] for i in range(7)]
pos = [[[x, yloc[i]] for x in np.linspace(xo[i], xrange[i] -xo[i], n_classes[i])] for i in range(7)]
colo = [[col[i] for x in range(n_classes[i])] for i in range(7)]
labels, positions, colors = [], [], []

for i in range(7):
    for j in range(n_classes[i]):
        labels.append(lab[i][j])
        positions.append(pos[i][j])
        colors.append(colo[i][j])

######################################################################
######################################################################
def get_alpha(From, To):
    if To[0] == From[0]:
        alpha = np.pi/4 - np.pi/4 * ( 1- np.sign(To[1] - From[1]))
#        print (alpha)
        return (alpha)
    if (To[0] >= From[0]) & (To[1] >= From[1]) :
        alpha = np.arctan((To[1] - From[1])/(To[0] - From[0]))
    if (To[0] <= From[0]) & (To[1] <= From[1]) :
        alpha = np.pi + np.arctan((To[1] - From[1])/(To[0] - From[0]))
    if (To[0] < From[0]) & (To[1] > From[1]) :
        alpha = np.pi + np.arctan((To[1] - From[1])/(To[0] - From[0]))
    if (To[0] > From[0]) & (To[1] < From[1]) :
        alpha = np.arctan((To[1] - From[1])/(To[0] - From[0]))
    return (alpha)
######################################################################

#props = dict(boxstyle='round', facecolor='wheat', alpha=1.0)


def Gapjunction(axe, n1, n2, radio, color, displ, lab):
    props = dict(boxstyle='round', facecolor=color, alpha=0.0, linewidth = 0.0)
    From, To = positions[n1], positions[n2]
    alpha = get_alpha(From, To)
    xt, yt = From[0] + displ*(To[0]-From[0]) -0.52, From[1] + displ*(To[1]-From[1]) -0.0
    axe.text(xt, yt , lab, fontsize = 24, ha = 'center', va = 'center', bbox=props)
    To_edge     = (To[0] - radio * np.cos(alpha), To[1] - radio * np.sin(alpha))
    From_edge   = (From[0] + radio * np.cos(alpha), From[1] + radio * np.sin(alpha))
    arrowprops=dict(arrowstyle='|-|, widthA=1.4,widthB=1.4', linewidth = 12, edgecolor = color, connectionstyle = "arc3, rad = %.3f"%0.0)
    axe.annotate("",  xy= To_edge, xycoords='data', xytext = From_edge, textcoords = 'data', arrowprops=arrowprops, zorder=0)






