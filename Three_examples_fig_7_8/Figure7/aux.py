import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
nrods = 51

##############################################################
##### Curvature loocking for the radio of cincurscribed circle
def get_curvature_angle(p1, p2, p3):
    p1, p2, p3 = np.array(p1), np.array(p2), np.array(p3)
    a = np.linalg.norm(p1 - p2)
    b = np.linalg.norm(p2 - p3)
    c = np.linalg.norm(p3 - p1)
    avec, bvec = p1 - p2, p3 - p2
    cosangle = np.dot(avec, bvec)/(a * b)
    sign = np.sign(avec[0]*bvec[1] - avec[1]*bvec[0])
    if (-1 < cosangle < 1):
        return sign *((a+b+c) * (-a+b+c) * (a-b+c) * (a+b-c))**(0.5) / (a*b*c)
    else:
        return 0.0
##############################################################
def curvature_body(body):
    nt = len(body)    
    curv_rad = np.zeros((nt, 49))
    for t in range(nt):
        x = body[t, range(1, 154, 3)]
        y = body[t, range(2, 154, 3)]
        for s in range(49):
            p1, p2, p3 = (x[s], y[s]), (x[s+1], y[s+1]),(x[s+2], y[s+2])
            curv_rad[t, s] = get_curvature_angle(p1, p2, p3)
    return curv_rad
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
    
#    
#    arrowprops=dict(arrowstyle='simple, tail_width=%.2f,head_width=%.2f,head_length=%.2f'%(arrow_shape), facecolor = color, edgecolor = 'k', linewidth = 0.6, connectionstyle = shape)
#    axe.annotate("",  xy= To_edge, xycoords='data', xytext = From_edge, textcoords = 'data', arrowprops=arrowprops, zorder = order)

def Gapjunction(axe, n1, n2, shape, radio, color, order):
    From, To = positions[n1], positions[n2]
    alpha = get_alpha(From, To)
    To_edge     = (To[0] - radio * np.cos(alpha), To[1] - radio * np.sin(alpha))
    From_edge   = (From[0] + radio * np.cos(alpha), From[1] + radio * np.sin(alpha))
    arrowprops=dict(arrowstyle='|-|, widthA=1.08,widthB=1.08', linewidth = 9, facecolor = color, edgecolor = 'k', connectionstyle =  shape)
    axe.annotate("",  xy= To_edge, xycoords='data', xytext = From_edge, textcoords = 'data', arrowprops=arrowprops, zorder=order)
    arrowprops=dict(arrowstyle='|-|, widthA=1.02,widthB=1.02', linewidth = 8, facecolor = color, edgecolor = color, connectionstyle =  shape)
    axe.annotate("",  xy= To_edge, xycoords='data', xytext = From_edge, textcoords = 'data', arrowprops=arrowprops, zorder=order)

def plot_segments(axnet, color, neurons, synapses, phen):
    circles = ['c%i' for i in range(len(labels))]
    fcolor = ['0.95'] * 15
    textcolor = ['0.05'] * 15
    for neu in neurons: 
        fcolor[neu] = color
        textcolor[neu] = '0.95'
    for n in np.arange(len(labels)):
        circles[n] = plt.Circle(positions[n], radius = R, fc = fcolor[n], edgecolor = 'k', lw = 1, zorder = 20)
        circles[n].set_alpha(1.00 )
        if n in neurons:
            axnet.text(positions[n][0], positions[n][1], labels[n][-4:], color = textcolor[n], ha = 'center', va = 'center', fontsize = 32, weight = 'bold', zorder = 22)
        else:
            axnet.text(positions[n][0], positions[n][1], labels[n][-4:], color = textcolor[n], ha = 'center', va = 'center', fontsize = 32, zorder = 22)
        axnet.add_patch(circles[n])
    ##### DRAW CONNECTIONS ###########
    ablated_fitness = ['0.99'] * 15
    zorder = [10] * 15
    for i in synapses: ablated_fitness[i] = color 
    for i in synapses: zorder[i] = 12 

    for seg in [0, 7]:
        connect(axnet, AS+seg, DA+seg, "arc3,rad= 0.00", Rs, ablated_fitness[0], zorder[0], np.sign(phen[0]))
        connect(axnet, AS+seg, VD+seg, "arc3,rad= 0.00", Rs, ablated_fitness[1], zorder[1], np.sign(phen[1]))
        connect(axnet, DA+seg, DB+seg, "arc3,rad= 0.00", Rs, ablated_fitness[2], zorder[2], np.sign(phen[2]))
        connect(axnet, DB+seg, AS+seg, "arc3,rad= 0.00", Rs, ablated_fitness[3], zorder[3], np.sign(phen[3]))
        connect(axnet, VD+seg, VA+seg, "arc3,rad= 0.00", Rs, ablated_fitness[4], zorder[4], np.sign(phen[4]))
        connect(axnet, VD+seg, VB+seg, "arc3,rad= 0.00", Rs, ablated_fitness[5], zorder[5], np.sign(phen[5]))
        connect(axnet, DA+seg, DD+seg, "arc3,rad= 0.00", Rs, ablated_fitness[6], zorder[6], np.sign(phen[6]))
        connect(axnet, VB+seg, DD+seg, "arc3,rad= 0.00", Rs, ablated_fitness[7], zorder[7], np.sign(phen[7]))
        connect(axnet, VA+seg, DD+seg, "arc3,rad= 0.00", Rs, ablated_fitness[8], zorder[8], np.sign(phen[8]))
        Gapjunction(axnet, VD+seg, DD+seg, "arc3,rad= 0.0", Rs, ablated_fitness[9], zorder[9])

    connect(axnet, DB, DD+seg, "arc3,rad= 0.05", Rs, ablated_fitness[10], zorder[10], np.sign(phen[10]))
    connect(axnet, VA+seg, DD, "arc3,rad= 0.00", Rs, ablated_fitness[11], zorder[11], np.sign(phen[11]))
    Gapjunction(axnet, AS, VA+seg, "arc3,rad=-0.15", Rs, ablated_fitness[12], zorder[12])
    Gapjunction(axnet, DA, AS+seg, "arc3,rad= -0.03", Rs, ablated_fitness[13], zorder[13])
    Gapjunction(axnet, VB, DB+seg, "arc3,rad=-0.08", Rs, ablated_fitness[14], zorder[14])

    axnet.set_xlim(-0.5, 30.5)
    axnet.set_ylim(-0.5 , 19.5)
    axnet.axis('off')
    axnet.set_aspect('equal')


r = 2* 40 *10**-6 * np.abs(np.sin(np.arccos((np.arange(51) - 25)/(25+0.2))))

def plot_worm(ax0,ax1, body):    
    fzl = 28
    x = body[:, range(1, 154, 3)]
    xc = np.sum(x, axis = 1)/nrods
    y = body[:, range(2, 154, 3)]
    yc = np.sum(y, axis = 1)/nrods
    a = body[:, range(3, 154, 3)]
    v = np.sqrt((np.diff(xc))**2 + (np.diff(yc))**2)/0.05
    curv = curvature_body(body)
    ############ Curvature  #######
    ###############################
    ax1.imshow(curv.T/1000., cmap=plt.get_cmap('seismic'), aspect='auto', vmin = -10, vmax = 10, origin = 'lower')
    ax1.set_xticks([0, 60, 120, 180, 240, 300])
    ax1.set_xticklabels(['0', '3', '6', '9', '12', '15'], fontsize = 24)
    ax1.set_yticks([-0.5, 47.5])
    ax1.set_yticklabels(['H', 'T'], fontsize = 26)
    ax1.set_xlabel('Time (s)', fontsize = 28)


