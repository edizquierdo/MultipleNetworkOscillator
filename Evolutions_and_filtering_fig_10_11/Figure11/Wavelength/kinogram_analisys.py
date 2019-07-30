import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
nrods = 51

def plot_lreg(ax, cl):
    x = cl[:,0]
    y = cl[:,1]
    A = np.vstack([x, np.ones(len(x))]).T
    m, c = np.linalg.lstsq(A, y)[0]
    x = np.linspace(0,1200, 5000)
    ax.plot(x, m*x + c, '--k', lw = 2, zorder = 4)
    return m, c

wavelength = np.zeros(104)

for ind in range(1): #[[100]:#
        plt.close('all')
        fig = plt.figure(figsize = [12,6])
        #axes
        gm = gridspec.GridSpec(100, 100)
        ax1 = plt.subplot(gm[10:-10, 10:-10])
        
        curv = np.transpose(np.loadtxt('../behavior/selected/%s/curv.dat'%ind))[1:,1:800]
        print (ind)
    ############ Curvature  #######
    ###############################
        ax1.imshow(curv, cmap=plt.get_cmap('seismic'), aspect='auto', vmin = -10, vmax = 10, origin = 'lower')
        ax1.set_xticks([])
        ax1.set_yticks([-0.5, 22.5])
        ax1.set_yticklabels(['H', 'T'], fontsize = 22)

        ax1.set_ylim(-0.5, 22.5)
        
        ##### Curvature regression
        cs = np.array(np.where(curv>0.0))[0] # curvature segment data
        ct = np.array(np.where(curv>0.0))[1] # curvature time data
        
        ### Initial time of a complete curve
        classes = [[] for i in range(42)]
        test = 45.5 # expected bins in a cycle at a frequency of 0.44 hz
        changes = []
        ##########  linear regresion data  ####
        lstch = 0
        for i in range(23): ## get high curvature points in each body segment
            to = ct[np.where(cs == i)]
            sg = cs[np.where(cs == i)]
            dto = np.diff(to)
            changes = to[np.where(dto > 3)[0]]
            k = 0
            for ch in changes:
                if ((i == 0) and (ch < lstch + test/5)): continue
                if ((i > 0)  and (ch < lstch - test/5)): continue
                if k == 0: lstch = ch
                classes[k].append([ch, i])
                ax1.plot([ch], [i], '.g', ms = 4)
                k +=1

        wl= 0

        ceval = 10
        mset, cset = np.zeros(ceval), np.zeros(ceval)

        for i in range(ceval):
            mset[i], cset[i] = plot_lreg(ax1, np.array(classes[i]))

        for i in range(1, ceval):
            xeval = -cset[i]/mset[i]
            yeval = mset[i-1] * xeval + cset[i-1]
            plt.plot([xeval, xeval], [0, yeval], '.k', ms = 13)
            plt.plot([xeval, xeval], [0, yeval], '.g', ms = 9)
            wl += (yeval)/22./(ceval -1)

        wavelength[ind] = wl

        ax1.set_xlim(0,350)

#####################################
        plt.figtext(0.053, 0.93, 'wavelength ~= %.2f worm length'%(wl), ha = 'left', va = 'center', fontsize = 20)

        fig.subplots_adjust(left=0.05, bottom=0.04, right=0.98, top=0.96)
        plt.savefig('wavelength_figures/plot_sel_%d.png'%ind)
        
### The script fails in this 7 solution, they were calculated manually.
#wavelength[48] = 0.7
#wavelength[52] = 0.75
#wavelength[60] = 0.9
#wavelength[67] = 0.8
#wavelength[72] = 0.8
#wavelength[97] = 0.85
#wavelength[100] = 0.8

#np.savetxt('kinogram_wavelength', wavelength)


#wavelength = np.loadtxt('kinogram_wavelength')



#plt.close('all')

#fig = plt.figure(figsize = [12,6])

#gm = gridspec.GridSpec(100, 190)
#ax2 = plt.subplot(gm[:, :])

#h = ax2.hist(wavelength, bins=20, orientation="vertical", color = 'b')
##ax2.set_xlim(0,1)
##ax2.yaxis.tick_right()
##ax2.yaxis.set_label_position("right")
#ax2.set_ylabel('Counts', fontsize = 18, rotation = 90, labelpad = 13)
##ax2.set_ylim(0, 110)
#ax2.set_xlabel('wavelength (body length)', fontsize = 18, labelpad = 2)
##ax0.set_xticks([])
##ax2.set_yticks([25, 50, 75, 100])
##ax2.get_children()[19].set_facecolor('r')

#plt.savefig('wavelength_distribution.png')
