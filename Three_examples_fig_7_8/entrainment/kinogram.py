import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
nrods = 51
import time
t = time.time()
ltrace = 400

################################################
def get_curvature(p1, p2, p3):
    p1, p2, p3 = np.array(p1), np.array(p2), np.array(p3)
    a = np.linalg.norm(p1 - p2)
    b = np.linalg.norm(p2 - p3)
    c = np.linalg.norm(p3 - p1)
    avec, bvec = p1 - p2, p3 - p2
    cosangle = np.dot(avec, bvec)/(a * b)
    sign = np.sign(avec[0]*bvec[1] - avec[1]*bvec[0])
    if (-1 < cosangle < 1):
        return sign *((a+b+c) * (-a+b+c) * (a-b+c) * (a+b-c))**(0.5) / (1000*a*b*c)
    else:
        return 0
def get_kinogram(body):
    curv = np.zeros((ltrace, 49))
    x = body[:, range(1, 154, 3)]
    y = body[:, range(2, 154, 3)]
    ##############################
    for t in range(ltrace):
        x = body[t, range(1, 154, 3)]
        y = body[t, range(2, 154, 3)]
        for s in range(49):
            p1, p2, p3 = (x[s], y[s]), (x[s+1], y[s+1]),(x[s+2], y[s+2])
            curv[t, s] = get_curvature(p1, p2, p3)
    return curv.T

labels = [r'$\Delta\phi = $' + angle for angle in ['0', '60', '120', '180', '240', '300', '360']]
bodylab = ['Head', '',  '',  'Midbody',  '',  'Tail',  'Tail']
delta_steps = [0, 8, 15, 23, 30, 38, 46]
def kinogram(ax, c1, c2, j):
    ############ Curvature  #######
    dstep = delta_steps[j]
    entrain = c1 - c2
    ax[0].imshow(c1[:, 20:180], cmap=plt.get_cmap('seismic'), aspect='auto', vmin = -10, vmax = 10, origin = 'lower')
    ax[1].imshow(entrain[:, 20:180], cmap=plt.get_cmap('BrBG'), aspect='auto', vmin = -10, vmax = 10, origin = 'lower')
    
    ax[0].set_yticks([4, 44])
    ax[0].set_yticklabels(['Head', 'Tail'], fontsize = 16)
    ax[0].axvline(x = 80, lw = 2, ls= '--', color = 'k')
    ax[0].set_xticks([0, 40, 80, 120])
    ax[0].set_xticklabels(['0', '2', '4', '6'], fontsize = 16)
    ax[0].set_ylabel(labels[j], fontsize = 23)

    ax[1].axis('off')
    if j == 0:
        ax[1].text(80, 25, 'Difference to\nfully entrained', ha='center', va='center', fontsize=21)
    else:
        ax[1].text(80, 53, 'Mean difference = %.1f'%np.mean(np.abs(entrain)), ha='center', va='center', fontsize=21)
    return (np.mean(np.abs(entrain)))
################################################
############# TAIL ANALYSIS AND FIGURES  #######
for ind in [30, 101, 103]:
    seg = 0 # 5= entrain tail; 0 = entrail
    evaluate = 'entrain_head'
    plt.close('all')
    fig = plt.figure(figsize = [12,15])
    gm = gridspec.GridSpec(700, 200)
    axes = [[plt.subplot(gm[10+i%7*100:90+i%7*100, 8+j*100:92+j*100]) for j in range(2)]  for i in range(7)]
    curvatures = [get_kinogram(np.loadtxt('%s/play/body_%d_%d_%d.dat'%(evaluate, ind, seg, rep))) for rep in range(6)]
    print(curvatures[0].shape)
    data_entrainment = np.zeros(7)
    
    for i in range(7):
        c = curvatures[i%6]
        data_entrainment[i] = kinogram(axes[i], curvatures[i%6], curvatures[0], i)
    ##############################
    np.savetxt('figures/data_%s_%d'%(evaluate, ind), data_entrainment)
    ##############################
    fig.subplots_adjust(left=0.06, bottom=0.06, right=0.98, top=0.88)
    plt.figtext(0.5,0.94,'Individual %d\nReplacing neural activity in the %s'%(ind, bodylab[seg]), ha='center', va='center', fontsize=23)
    plt.savefig('./figures/kinogram_%d_%s_seg_%d.png'%(ind, evaluate, seg+1))


#################################################
############## HEAD ANALYSIS AND FIGURES  #######
for ind in [30, 101, 103]:
    seg = 5 # 5= entrain tail; 0 = entrail
    evaluate = 'entrain_tail'
    plt.close('all')
    fig = plt.figure(figsize = [12,15])
    gm = gridspec.GridSpec(700, 200)
    axes = [[plt.subplot(gm[10+i%7*100:90+i%7*100, 8+j*100:92+j*100]) for j in range(2)]  for i in range(7)]
    curvatures = [get_kinogram(np.loadtxt('%s/play/body_%d_%d_%d.dat'%(evaluate, ind, seg, rep))) for rep in range(6)]
    print(curvatures[0].shape)
    data_entrainment = np.zeros(7)
    
    for i in range(7):
        data_entrainment[i] = kinogram(axes[i], curvatures[i%6], curvatures[0], i)
    ##############################
    np.savetxt('figures/data_%s_%d'%(evaluate, ind), data_entrainment)
    ##############################
    fig.subplots_adjust(left=0.06, bottom=0.06, right=0.98, top=0.88)
    plt.figtext(0.5,0.94,'Individual %d\nReplacing neural activity in the %s'%(ind, bodylab[seg]), ha='center', va='center', fontsize=23)
    plt.savefig('./figures/kinogram_%d_%s_seg_%d.png'%(ind, evaluate, seg+1))

