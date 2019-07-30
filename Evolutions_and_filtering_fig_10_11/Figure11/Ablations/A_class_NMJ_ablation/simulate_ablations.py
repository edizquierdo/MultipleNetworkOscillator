import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
nrods = 51

for ind in range(104):
    ######################################################################
    os.system('cp code/main selected/%d'%ind)
    os.chdir('selected/%d'%ind)
    os.system('./main')
    os.chdir('../../')
    #################     PLOT   ################
    ############## neural activity  ################
    [AS, DA, DB, DD,  VD, VB, VA] = [0, 1, 2, 3, 4, 5, 6]
    colors = ['b', 'g', 'k', 'r', 'y', 'orange', 'brown']
    activity = np.loadtxt('./selected/%d/act.dat'%ind)[:201, 1:50]

    plt.close('all')
    fig = plt.figure(figsize = [30,14])
    #axes
    gm = gridspec.GridSpec(100, 190)
    ax0 = plt.subplot(gm[:, 45:115])
    ax1 = plt.subplot(gm[10:52, :40])
    ax2 = plt.subplot(gm[58:, :40])

    gt = gridspec.GridSpecFromSubplotSpec(70, 90, subplot_spec = gm[5:-5, 130:])
    axesA = [plt.subplot(gt[i*10: i*10 + 8, 0:25]) for i in range(7)]
    axesB = [plt.subplot(gt[i*10: i*10 + 8, 30:55]) for i in range(7)]
    axesC = [plt.subplot(gt[i*10: i*10 + 8, 60:85]) for i in range(7)]

    for i in range(7):
        axesA[i].plot(activity[:,DB + i*7], color = colors[DB]) #DB head
        axesB[i].plot(activity[:,DA + i*7], color = colors[DA]) #DA head
        axesC[i].plot(activity[:,AS + i*7], color = colors[AS]) #AS head
        axesC[i].plot(activity[:,DD + i*7], color = colors[DD]) #DD head

        axesA[i].plot(activity[:,VB + i*7], color = colors[VB]) #VB head
        axesB[i].plot(activity[:,VA + i*7], color = colors[VA]) #VA head
        axesC[i].plot(activity[:,VD + i*7], color = colors[VD]) #VD head

    for ax in axesA:
        ax.axis('off')
        ax.set_ylim(-0.01, 1.1)
    for ax in axesB:
        ax.axis('off')
        ax.set_ylim(-0.01, 1.1)
    for ax in axesC:
        ax.axis('off')
        ax.set_ylim(-0.01, 1.1)

    fz = 21
    for k in [0.08, 0.93]:
        plt.figtext(0.71, k, 'DB', fontsize = fz, color = colors[DB])
        plt.figtext(0.74, k, 'VB', fontsize = fz, color = colors[VB])
        
        plt.figtext(0.82, k, 'DA', fontsize = fz, color = colors[DA])
        plt.figtext(0.85, k, 'VA', fontsize = fz, color = colors[VA])
        
        plt.figtext(0.90, k, 'AS', fontsize = fz, color = colors[AS])
        plt.figtext(0.93, k, 'DD', fontsize = fz, color = colors[DD])
        plt.figtext(0.96, k, 'VD', fontsize = fz, color = colors[VD])

    segment = ['Head', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6']
    for i in range(7):
        plt.figtext(0.67, 0.85 - .115 * i,  segment[i], fontsize = fz, color = 'k', ha = 'center')

    ##################     CALCULATE   ################
    ###################### CURVATURE  ################
    body = np.loadtxt('selected/%s/body.dat'%ind)[:201, :] ## first 50 seconds of simulation
    x = body[:, range(1, 154, 3)]
    xc = np.sum(x, axis = 1)/nrods
    y = body[:, range(2, 154, 3)]
    yc = np.sum(y, axis = 1)/nrods
    v = np.sqrt((np.diff(xc))**2 + (np.diff(yc))**2)/0.05

    curv = np.transpose(np.loadtxt('selected/%s/curv.dat'%ind))[1:,800:1001]
    print (np.max(curv), np.min(curv))
    ############ Trayectory #######
    ###############################
    k = 0
    for i in range(len(x)):
        k += 1
        if k == 5:
            ax0.plot(x[i], y[i], 'k', lw = 0.25, alpha = 0.5)
            ax0.plot(x[i][0], y[i][0], 'or', markeredgecolor = 'r', ms = 1)
            k = 0

    k=10
    for i in range(len(x)):
        k+=1
        if k == 60:
            ax0.plot(x[i], y[i], 'g', lw = 3.5)
            ax0.plot(x[i][0], y[i][0], 'or', ms = 6)
            ax0.plot(x[i][-1], y[i][-1], 'og', ms = 4)
            k=0

    ax0.set_xlim(np.min(x), np.max(x))
    ax0.set_ylim(np.min(y), np.max(y))
    ax0.set_xticks([])
    ax0.set_yticks([])
    ax0.set_aspect('equal')
    ax0.set_ylabel('Y')
    ax0.set_xlabel('X')

    ############ Curvature  #######
    ###############################
    ax1.imshow(curv, cmap=plt.get_cmap('seismic'), aspect='auto', vmin = -10, vmax = 10, origin = 'lower')
    ax1.set_xticks([])
    ax1.set_yticks([-0.5, 22.5])
    ax1.set_yticklabels(['H', 'T'], fontsize = 22)

    ############ Velocity #######
    ###############################
    ax2.plot(np.linspace(0, 10, len(v)), 1000*v, 'k')
    ax2.axhline(y = 0.22, linestyle =  '--', color = 'r')
    ax2.set_ylim(0, 0.5)
    ax2.set_xlim(0, 10)
    ax2.set_ylabel('Velocity mm/s')
    ax2.set_xlabel('Time')

    fig.subplots_adjust(left=0.05, bottom=0.04, right=0.98, top=0.96)
    plt.savefig('kinogram_figures/plot_sel_%d.png'%ind)
