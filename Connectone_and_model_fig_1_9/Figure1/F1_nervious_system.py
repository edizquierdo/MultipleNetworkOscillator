import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as path_effects
#######################################################################
######################################################################
[AS, DA, DB, DD, VD, VB, VA] = [0, 1, 2, 3, 4, 5, 6]

object_fontsize = 45

neurons = ['AS', 'DA', 'DB', 'DD', 'VD', 'VB', 'VA']
n_classes = [7, 7, 7, 7, 7, 7, 7]

col = ['#cc66cc','#d8b365', '#5ab4ac', '#42a4d1', '#99ccff', '#01665e', '#8c510a']

xo = [5, 1, 9, 7.5, 2.5, 8.5, 1.5]

yloc = [18, 13, 13, 10, 7.5, 3, 2]
unitpos = [-5, 2,     17.0, 17.0, 17.0,         17.3, 14.5]

labels, positions, colors = [], [], []
lab = [['%s'%(neurons[i]) for numb in range(1, 1+n_classes[i])] for i in range(7)]
pos = [[[x, yloc[i]] for x in xo[i]+np.cumsum(unitpos)] for i in range(7)]
colo = [[col[i] for x in range(n_classes[i])] for i in range(7)]
for j in range(7):
    for i in range(7):
        labels.append(lab[i][j])
        positions.append(pos[i][j])
        colors.append(colo[i][j])

mus_loc_x = np.linspace(-26.5, 104, 24)
muscle_location = []
muscle_names = []
for i in range(12):
    muscle_location.append([mus_loc_x[i*2], 20.70, 'DM %.2d'%(i*2 + 1)])
    muscle_names.append('D %.2d'%(i*2 + 1))
    muscle_location.append([mus_loc_x[i*2], -0.70, 'VM %.2d'%(i*2 + 1)])
    muscle_names.append('V %.2d'%(i*2 + 1))
    muscle_location.append([mus_loc_x[i*2+1], 22.00, 'DM %.2d'%(i*2 + 2)])
    muscle_names.append('D %.2d'%(i*2 + 2))
    muscle_location.append([mus_loc_x[i*2+1], -2.00, 'VM %.2d'%(i*2 + 2)])
    muscle_names.append('V %.2d'%(i*2 + 2))

muscles_names = list(np.array(muscle_location)[:,2])
#######################################################################
######################################################################
def get_alpha(From, To):
    if To[0] == From[0]: alpha = np.pi/4 - np.pi/4 * ( 1- np.sign(To[1] - From[1]))
    if (To[0] >= From[0]) & (To[1] >= From[1]) : alpha = np.arctan((To[1] - From[1])/(To[0] - From[0]))
    if (To[0] <= From[0]) & (To[1] <= From[1]) : alpha = np.pi + np.arctan((To[1] - From[1])/(To[0] - From[0]))
    if (To[0] < From[0]) & (To[1] > From[1]) : alpha = np.pi + np.arctan((To[1] - From[1])/(To[0] - From[0]))
    if (To[0] > From[0]) & (To[1] < From[1]) : alpha = np.arctan((To[1] - From[1])/(To[0] - From[0]))
    return (alpha)

def conectS_bold(axe, n1, n2, radio, color, curv):
    From, To = positions[n1], positions[n2]
    alpha = get_alpha(From, To)
    To_edge     = (To[0] - radio * np.cos(alpha), To[1] - radio * np.sin(alpha))
    From_edge   = (From[0] + radio * np.cos(alpha), From[1] + radio * np.sin(alpha))
    arrow_shape = (1.7, 3.7, 3.7)
    arrowprops=dict(arrowstyle='simple, tail_width=%.2f,head_width=%.2f,head_length=%.2f'%(arrow_shape), facecolor = color, edgecolor = 'k', linewidth = 0.0, connectionstyle = "arc3, rad = %.3f"%curv)
    axe.annotate("",  xy= To_edge, xycoords='data', xytext = From_edge, textcoords = 'data', arrowprops=arrowprops, zorder = 2)

def Gapjunction(axe, n1, n2, radio, color, curv):
    From, To = positions[n1], positions[n2]
    alpha = get_alpha(From, To)
    To_edge     = (To[0] - radio * np.cos(alpha), To[1] - radio * np.sin(alpha))
    From_edge   = (From[0] + radio * np.cos(alpha), From[1] + radio * np.sin(alpha))
    arrowprops=dict(arrowstyle='|-|, widthA=1.58,widthB=1.58', linewidth = 16, edgecolor = color, connectionstyle = "arc3, rad = %.3f"%curv)
    axe.annotate("",  xy= To_edge, xycoords='data', xytext = From_edge, textcoords = 'data', arrowprops=arrowprops, zorder=0)
#    arrowprops=dict(arrowstyle='|-|, widthA=1.00,widthB=1.00', linewidth = 2, edgecolor = 'w', connectionstyle = "arc3, rad = %.3f"%curv)
#    axe.annotate("",  xy= To_edge, xycoords='data', xytext = From_edge, textcoords = 'data', arrowprops=arrowprops, zorder=0)
##################################

######################################################################
#### Figure and axes definition
######################################################################
plt.close('all')
fig = plt.figure(figsize = [41, 20])
R, Rr, Rgj = 1.2, 1.23, 1.4
gmain = gridspec.GridSpec(100, 100)
ax = plt.subplot(gmain[:,:])
######################################################################
##########################

#### intraunit synapses ####
curv = 0.00
intraunit_color = 'k'
for seg in [3]:#range(7):
    conectS_bold(ax, AS+seg*7, DA+seg*7, Rr, intraunit_color, curv)
    conectS_bold(ax, AS+seg*7, VD+seg*7, Rr, intraunit_color, -0.02)
    conectS_bold(ax, DA+seg*7, DB+seg*7, Rr, intraunit_color, curv)
    conectS_bold(ax, DB+seg*7, AS+seg*7, Rr, intraunit_color, curv)
    conectS_bold(ax, VD+seg*7, VA+seg*7, Rr, intraunit_color, curv)
    conectS_bold(ax, VD+seg*7, VB+seg*7, Rr, intraunit_color, curv)
    conectS_bold(ax, DA+seg*7, DD+seg*7, Rr, intraunit_color, curv)
    conectS_bold(ax, VB+seg*7, DD+seg*7, Rr, intraunit_color, curv)
    conectS_bold(ax, VA+seg*7, DD+seg*7, Rr, intraunit_color, curv)
    
    Gapjunction(ax, VD+seg*7, DD+seg*7, Rgj, intraunit_color, curv)

#### interunit synapses ####
interunit_color = 'r'
for seg in [2, 3]:#range(6):
    conectS_bold(ax, DB+seg*7, DD+(seg+1)*7, Rr, interunit_color, curv)
    conectS_bold(ax, VA+(seg+1)*7, DD+seg*7, Rr, interunit_color, curv)
    Gapjunction(ax, AS+seg*7, VA+(seg+1)*7, Rgj, interunit_color, -0.2)
    Gapjunction(ax, DA+seg*7, AS+(seg+1)*7, Rgj, interunit_color, curv)
    Gapjunction(ax, VB+seg*7, DB+(seg+1)*7, Rgj, interunit_color, -0.10)

##################################
##### DRAW THE CIRCLES ###########
##################################
fcolors = ['w','k', 'k', 'w', 'k', 'w', 'w']
circles = ['c%i' for i in range(len(labels))]
alphas = np.ones(49)*0.3
for i in range(21, 28): alphas[i] = 1.0
for n in np.arange(len(labels)):
    circles[n] = plt.Circle(positions[n], radius = R, fc = colors[n], edgecolor = 'k', lw = 1, alpha = alphas[n])
#    circles[n].set_alpha(alphas[i])
    t = ax.text(positions[n][0], positions[n][1], labels[n], color = fcolors[n%7], ha = 'center', va = 'center', fontsize = 65, alpha = alphas[n])
#    t.set_path_effects([path_effects.Stroke(linewidth=5, foreground='black'),
#                       path_effects.Normal()])
    ax.add_patch(circles[n])
    

####################
### Draw muscles  ##
####################
#for m in [12, 13, 32, 33]: muscle_location[m][2] = ''
blue = ['#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b']
mfc, mec = '0.20', '0.50'
mfc, mec = '#a6611a', '#018571'
mfc, mec = '#a63603', '#fd8d3c'
mfc, mec = blue[6], blue[8]

mwidth, mhight = 11.1, 1.25
for i in [4, 5, 6]: colors[i] = '0.15'
for k, muscle in enumerate(muscle_location):
    if ((k>13)*(k<32)):
        if int(k/2) in [4, 5, 6, 10, 11, 12, 16, 17, 18, 19]:
            muscle.append(mfc)
        else:
            muscle.append(mec)
        ax.add_patch(mpatch.Ellipse(xy=[muscle[0], muscle[1]], width= mwidth, height=mhight, facecolor = muscle[3], edgecolor = 'k', alpha = 1.0))
        ax.text(muscle[0], muscle[1], muscle[2], color = 'w', fontsize = 50, ha = 'center', va = 'center')


ax.set_xlim(7.5, 64.3)
ax.set_ylim(-3.5 , 23.5)
ax.set_aspect('equal')
ax.axis('off')

ds = np.diff(mus_loc_x)[0]/2. 
ds0 = -12.16
dsp = 0

segment_color = '0.45'

for i in range(50):
    ax.plot([ds0 + dsp + 0.2, ds0 + dsp + ds - 0.4], [-3.2, -3.2], lw = 19, color = segment_color)
    ax.plot([ds0 + dsp + 0.2, ds0 + dsp + ds - 0.4], [23.2, 23.2], lw = 19, color = segment_color)
    dsp += ds

print(dsp)

ax.text(9.5, 10, 'Head', fontsize = 56, ha = 'center', va = 'center')
ax.text(62.3, 10, 'Tail', fontsize = 56, ha = 'center', va = 'center')

'''####################
#### save figure   ####
#####################'''
fig.subplots_adjust(left=0, bottom=0.00, right=1.0, top=1.0)
plt.savefig('model.png')
plt.savefig('model.pdf')
plt.clf()
