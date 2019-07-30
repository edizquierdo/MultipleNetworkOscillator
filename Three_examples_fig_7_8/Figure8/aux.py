import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
##############################################################
######################################################################
neurons = ['AS', 'DA', 'DB', 'DD', 'VD', 'VB', 'VA']
[AS, DA, DB, DD, VD, VB, VA] = [0, 1, 2, 3, 4, 5, 6]

xloc = [5, 1, 9, 7.5, 2.5, 8.5, 1.5]
yloc = [18, 13, 13, 10, 7.5, 2, 1]
labels, positions = [], []
for seg in [0, 1]:
    for i in range(7):
        labels.append('%s'%(neurons[i]))
        positions.append([20 * seg + xloc[i], yloc[i]])
        
R = 1.2
Rs = 1.35
######################################################################
######################################################################
def get_alpha(From, To):
    if To[0] == From[0]: alpha = np.pi/4 - np.pi/4 * ( 1- np.sign(To[1] - From[1]))
    if (To[0] >= From[0]) & (To[1] >= From[1]) : alpha = np.arctan((To[1] - From[1])/(To[0] - From[0]))
    if (To[0] <= From[0]) & (To[1] <= From[1]) : alpha = np.pi + np.arctan((To[1] - From[1])/(To[0] - From[0]))
    if (To[0] < From[0]) & (To[1] > From[1]) : alpha = np.pi + np.arctan((To[1] - From[1])/(To[0] - From[0]))
    if (To[0] > From[0]) & (To[1] < From[1]) : alpha = np.arctan((To[1] - From[1])/(To[0] - From[0]))
    return (alpha)
######################################################################
def connect(axe, n1, n2, shape, radio, color, order, sign):
    radio2 = radio
    if sign == -1: radio2 = 1.65
    From, To = positions[n1], positions[n2]
    alpha = get_alpha(From, To)
    To_edge     = (To[0] - radio2 * np.cos(alpha), To[1] - radio2 * np.sin(alpha))
    From_edge   = (From[0] + radio * np.cos(alpha), From[1] + radio * np.sin(alpha))
    arrow_shape = (0.8, 2, 2)
    if sign == -1:
        l = axe.annotate("",  xy= To_edge, xycoords='data', xytext = From_edge, textcoords = 'data', arrowprops=dict(arrowstyle='simple, tail_width = %s, head_width = %s, head_length = %s'%(0.8, 0.1, 0.1), facecolor = color, edgecolor = 'k', connectionstyle = shape), alpha = 1.0, zorder = order)
        c = plt.Circle(To_edge, radius = 0.35, fc = color, edgecolor = 'k', alpha = 1.0, zorder = order+2)
        cir = axe.add_patch(c)
    else:
        arrowprops=dict(arrowstyle='simple, tail_width=%.2f,head_width=%.2f,head_length=%.2f'%(arrow_shape), facecolor = color, edgecolor = 'k', linewidth = 0.6, connectionstyle = shape)
        axe.annotate("",  xy= To_edge, xycoords='data', xytext = From_edge, textcoords = 'data', arrowprops=arrowprops, zorder = order)
######################################################################
def Gapjunction(axe, n1, n2, shape, radio, color, order):
    From, To = positions[n1], positions[n2]
    alpha = get_alpha(From, To)
    To_edge     = (To[0] - radio * np.cos(alpha), To[1] - radio * np.sin(alpha))
    From_edge   = (From[0] + radio * np.cos(alpha), From[1] + radio * np.sin(alpha))
    arrowprops=dict(arrowstyle='|-|, widthA=1.08,widthB=1.08', linewidth = 9, facecolor = color, edgecolor = 'k', connectionstyle =  shape)
    axe.annotate("",  xy= To_edge, xycoords='data', xytext = From_edge, textcoords = 'data', arrowprops=arrowprops, zorder=order)
    arrowprops=dict(arrowstyle='|-|, widthA=1.02,widthB=1.02', linewidth = 8, facecolor = color, edgecolor = color, connectionstyle =  shape)
    axe.annotate("",  xy= To_edge, xycoords='data', xytext = From_edge, textcoords = 'data', arrowprops=arrowprops, zorder=order)

######################################################################
######################################################################
def plot_segments(axnet, color, neurons, synapses, phen):
    circles = ['c%i' for i in range(len(labels))]
    fcolor = ['0.95'] * 15
    textcolor = ['0.65'] * 15
    for neu in neurons: 
        fcolor[neu] = color
        textcolor[neu] = '0.95'
    for n in np.arange(len(labels)):
        circles[n] = plt.Circle(positions[n], radius = R, fc = fcolor[n], edgecolor = 'k', lw = 1, zorder = 20)
        circles[n].set_alpha(1.00 )
        if n in neurons:
            axnet.text(positions[n][0], positions[n][1], labels[n][-4:], color = textcolor[n], ha = 'center', va = 'center', fontsize = 32, weight='bold', zorder = 22)
        else:
            axnet.text(positions[n][0], positions[n][1], labels[n][-4:], color = textcolor[n], ha = 'center', va = 'center', fontsize = 32, zorder = 22)
        axnet.add_patch(circles[n])
    ##### DRAW CONNECTIONS ###########
    ablated_fitness = [color] * 15
    zorder = [10] * 15
    for i in synapses: ablated_fitness[i] = color 
    for i in synapses: zorder[i] = 12 
    for syn in synapses:
        seg = 0
        if (syn == 0):
            connect(axnet, AS+seg, DA+seg, "arc3,rad= 0.00", Rs, ablated_fitness[0], zorder[0], np.sign(phen[0]))
        if (syn == 1):
            connect(axnet, AS+seg, VD+seg, "arc3,rad= 0.00", Rs, ablated_fitness[1], zorder[1], np.sign(phen[1]))
        if (syn == 2):
            connect(axnet, DA+seg, DB+seg, "arc3,rad= 0.00", Rs, ablated_fitness[2], zorder[2], np.sign(phen[2]))
        if (syn == 3):
            connect(axnet, DB+seg, AS+seg, "arc3,rad= 0.00", Rs, ablated_fitness[3], zorder[3], np.sign(phen[3]))
        if (syn == 4):
            connect(axnet, VD+seg, VA+seg, "arc3,rad= 0.00", Rs, ablated_fitness[4], zorder[4], np.sign(phen[4]))
        if (syn == 5):
            connect(axnet, VD+seg, VB+seg, "arc3,rad= 0.00", Rs, ablated_fitness[5], zorder[5], np.sign(phen[5]))
        seg = 7
        if (syn == 6):
            connect(axnet, AS+seg, DA+seg, "arc3,rad= 0.00", Rs, ablated_fitness[0], zorder[0], np.sign(phen[0]))
        if (syn == 7):
            connect(axnet, AS+seg, VD+seg, "arc3,rad= 0.00", Rs, ablated_fitness[1], zorder[1], np.sign(phen[1]))
        if (syn == 8):
            connect(axnet, DA+seg, DB+seg, "arc3,rad= 0.00", Rs, ablated_fitness[2], zorder[2], np.sign(phen[2]))
        if (syn == 9):
            connect(axnet, DB+seg, AS+seg, "arc3,rad= 0.00", Rs, ablated_fitness[3], zorder[3], np.sign(phen[3]))
        if (syn == 10):
            connect(axnet, VD+seg, VA+seg, "arc3,rad= 0.00", Rs, ablated_fitness[4], zorder[4], np.sign(phen[4]))
        if (syn == 11):
            connect(axnet, VD+seg, VB+seg, "arc3,rad= 0.00", Rs, ablated_fitness[5], zorder[5], np.sign(phen[5]))

        if (syn == 12):
            Gapjunction(axnet, AS, VA+seg, "arc3,rad=-0.15", Rs, ablated_fitness[12], zorder[12])
        if (syn == 13):
            Gapjunction(axnet, DA, AS+seg, "arc3,rad= -0.03", Rs, ablated_fitness[13], zorder[13])
        if (syn == 14):
            Gapjunction(axnet, VB, DB+seg, "arc3,rad=-0.08", Rs, ablated_fitness[14], zorder[14])

    axnet.set_xlim(-0.5, 30.5)
    axnet.set_ylim(-0.5 , 19.5)
    axnet.axis('off')
    axnet.set_aspect('equal')


