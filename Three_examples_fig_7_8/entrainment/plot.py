import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np


################################################
################     DATA     ##################
anterio30 = np.loadtxt('figures/data_entrain_head_30')
anterio101 = np.loadtxt('figures/data_entrain_head_101')
anterio103 = np.loadtxt('figures/data_entrain_head_103')

posterio30 = np.loadtxt('figures/data_entrain_tail_30')
posterio101 = np.loadtxt('figures/data_entrain_tail_101')
posterio103 = np.loadtxt('figures/data_entrain_tail_103')


################################################
plt.close('all')
fig = plt.figure(figsize = [16,4])
#axes
gm = gridspec.GridSpec(70, 290)
ax1 = plt.subplot(gm[:, 0:90])
ax2 = plt.subplot(gm[:, 100:190])
ax3 = plt.subplot(gm[:, 200:290])


ax1.plot(anterio101, '--dr', label = 'Head to tail entrainment')
ax1.plot(posterio101, '--og', label = 'Tail to head entrainment')


ax2.plot(anterio103, '--dr', label = 'Head to tail entrainment')
ax2.plot(posterio103, '--og', label = 'Tail to head entrainment')

ax3.plot(anterio30, '--dr', label = 'Head to tail entrainment')
ax3.plot(posterio30, '--og', label = 'Tail to head entrainment')

ind = [101, 103, 30]

for k, ax in enumerate([ax1, ax2, ax3]):
    ax.set_ylim(0, 360)
    ax.set_ylabel('Entrainment evaluation')
    ax.set_xticks([0, 1, 2, 3, 4, 5])
    ax.set_xticklabels(['0', '60', '120', '180', '240', '300'])
    ax.set_xlabel('Degrees entrainment displacement')
    ax.set_title('Individual %s'%(str(ind[k])))


ax3.legend()

plt.savefig('./entrainment.png')
