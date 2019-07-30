import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
nrods = 51



##### Curvature loocking for the radio of cincurscribed circle ############
def get_curvature_angle(p1, p2, p3):
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
########################################################################


selected = np.loadtxt('list_selected')
os.chdir('./script')
os.system('make clean')
os.system('make')


for ind in selected:#[[4]:#4]:#
    os.system('cp ../../../filtered_agents/%d/best.gen.dat ./'%(ind))
    os.system('./main')
    
    curv = np.zeros((602, 49))
    body = np.loadtxt('./body.dat') ## first 100 seconds of simulation
    for t in range(602):
        x = body[t, range(1, 154, 3)]
        y = body[t, range(2, 154, 3)]
        for s in range(49):
            p1, p2, p3 = (x[s], y[s]), (x[s+1], y[s+1]),(x[s+2], y[s+2])
            curv[t, s] = get_curvature_angle(p1, p2, p3)
#        print(t)
    np.savetxt('curv_%d'%ind, curv)
    plt.close('all')
    fig = plt.figure(figsize = [15,8])
    #axes
    gm = gridspec.GridSpec(100, 100)
    ax1 = plt.subplot(gm[10:-10, 10:-10])
    
    print (ind)
    ############ Curvature  #######
    ax1.imshow(curv.T, cmap=plt.get_cmap('seismic'), aspect='auto', vmin = -10, vmax = 10, origin = 'lower', interpolation = 'none')
    ax1.set_xticks([])
    ax1.set_yticks([6, 16, 41])
    ax1.set_yticklabels(['Head', 'Neck', 'Tail'], fontsize = 22)
    ax1.set_xticks([0, 200, 400, 600])
    ax1.set_xticklabels(['0', '10', '20', '30'], fontsize = 22)
    ax1.set_xlim(0, 600)
    ax1.set_ylim(-0.5, 48.5)
    ax1.set_xlabel('Time (s)', fontsize = 22)
#    ax1.set_ylabel('Body axis', fontsize = 22)
    plt.figtext(0.5, 0.955, 'Inactivated neck BWM2, Ind %d'%ind, fontsize = 24, va = 'center', ha = 'center')
    
    #calculate amplitud as in Fouad 2017
    amplitude_fouad = np.zeros(3)
    head, neck, tail = 7, 17, 42
    L = 100.
    p = 0
    for position in [head, neck, tail]:
        dt = 1/20.
        dk = np.diff(curv[:,position])
        amplitude_fouad[p] = np.sqrt(L*np.mean(dk**2)/dt)
        p+=1
    np.savetxt('../../data/amplitude_fouad_bwm_%d'%ind, amplitude_fouad)
    
    
    plt.figtext(0.5, 0.90, 'amplitude Head = %.1f, Neck = %.1f, Tail = %.1f'%(amplitude_fouad[0], amplitude_fouad[1], amplitude_fouad[2]), fontsize = 24, va = 'center', ha = 'center')

    fig.subplots_adjust(left=0.05, bottom=0.04, right=0.95, top=0.94)
    
    plt.savefig('../kinogram_%d'%ind)
os.chdir('../')



