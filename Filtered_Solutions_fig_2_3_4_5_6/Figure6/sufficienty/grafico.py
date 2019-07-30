import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

selected = [5, 9, 30, 42, 45, 53, 4, 21, 23, 37, 63, 66, 68, 101, 103]

plt.close('all')
fig = plt.figure(0, figsize = (14, 7))

gs = gridspec.GridSpec(100, 100)

ax11 = plt.subplot(gs[50:, :30])
ax12 = plt.subplot(gs[50:, 35:65])
ax13 = plt.subplot(gs[50:, 70:100])
axes = [ax11, ax12, ax13]

ax1 = plt.subplot(gs[:50, :30])
ax2 = plt.subplot(gs[:50, 35:65])
ax3 = plt.subplot(gs[:50, 70:100])
axes2 = [ax1, ax2, ax3]

full = np.loadtxt('data/data_full')
groups = ['minimal_9', 'minimal_45', 'minimal_23']
groups_label = ['Minimal i', 'Minimal ii', 'Minimal iii']

data = np.array([np.mean(np.loadtxt('data/data_%s'%s), axis = 1)/np.mean(full, axis = 1) for s in groups])
data[np.where(data<0)] = 0 
cap = dict(ecolor='black', lw=2, capsize=5, capthick=2)

for k, dat in enumerate(data): 
    axes[k].bar([0], np.mean(dat[:3]),  yerr=np.std(dat[:3])/np.sqrt(3),  color = 'purple', error_kw=cap)
    axes[k].bar([1], np.mean(dat[3:6]), yerr=np.std(dat[3:6])/np.sqrt(3), color = 'red', error_kw=cap)
    axes[k].bar([2], np.mean(dat[6:]),  yerr=np.std(dat[6:])/np.sqrt(9),  color = 'green', error_kw=cap)
    
    axes[k].set_xlim(-0.2, 3)
    axes[k].set_ylim(0, 0.88)
    axes[k].set_xticks([0.4, 1.4, 2.4])
    axes[k].set_xticklabels(['i', 'ii', 'iii'], fontsize = 16)
#    axes[k].set_xlabel(groups_label[k], fontsize = 16)
    axes[k].set_xlabel("Coordination groups", fontsize = 16)

axes[0].set_ylabel('Locomotion fitness', fontsize = 16)
#################################################################
ax1.imshow(plt.imread('./networks/Network_9.png'))
ax2.imshow(plt.imread('./networks/Network_45.png'))
ax3.imshow(plt.imread('./networks/Network_23.png'))
for ax in axes2: ax.axis('Off')
#################################################################
#################################################################
fz = 32
plt.figtext(0.015, 0.965, 'A', ha = 'center', va = 'center', fontsize = fz)
plt.figtext(0.015, 0.57, 'B', ha = 'center', va = 'center', fontsize = fz)
plt.figtext(0.35, 0.965, 'C', ha = 'center', va = 'center', fontsize = fz)
plt.figtext(0.35, 0.57, 'D', ha = 'center', va = 'center', fontsize = fz)
plt.figtext(0.68, 0.965, 'E', ha = 'center', va = 'center', fontsize = fz)
plt.figtext(0.68, 0.57, 'F', ha = 'center', va = 'center', fontsize = fz)


fz = 20
plt.figtext(0.19, 0.97, 'Anterior' + r'$\leftarrow$' + 'Posterior', ha = 'center', va = 'center', fontsize = fz)
plt.figtext(0.52, 0.97, 'Anterior' + r'$\leftrightarrow$' + 'Posterior', ha = 'center', va = 'center', fontsize = fz)
plt.figtext(0.84, 0.97, 'Anterior' + r'$\rightarrow$' + 'Posterior', ha = 'center', va = 'center', fontsize = fz)
#################################################################
fig.subplots_adjust(left=0.06, bottom=0.11, right=0.99, top=0.98)
plt.savefig('resumen.png', dpi = 300)

##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
#import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec

#selected = [5, 9, 30, 42, 45, 53, 4, 21, 23, 37, 63, 66, 68, 101, 103]

#plt.close('all')
#fig = plt.figure(0, figsize = (18, 7))

#gs = gridspec.GridSpec(100, 140)

#ax11 = plt.subplot(gs[50:, :30])
#ax12 = plt.subplot(gs[50:, 35:65])
#ax13 = plt.subplot(gs[50:, 70:100])
#ax14 = plt.subplot(gs[50:, 105:135])
#axes = [ax11, ax12, ax13, ax14]

#ax1 = plt.subplot(gs[:50, :30])
#ax2 = plt.subplot(gs[:50, 35:65])
#ax3 = plt.subplot(gs[:50, 70:100])
#ax4 = plt.subplot(gs[:50, 105:135])
#axes2 = [ax1, ax2, ax3, ax4]

#full = np.loadtxt('data/data_full')
#groups = ['minimal_9', 'minimal_23', 'minimal_45a', 'minimal_45b']

#data = np.array([np.mean(np.loadtxt('data/data_%s'%s), axis = 1)/np.mean(full, axis = 1) for s in groups])
#data[np.where(data<0)] = 0 
#cap = dict(ecolor='black', lw=2, capsize=5, capthick=2)

#for k, dat in enumerate(data): 
#    axes[k].bar([0], np.mean(dat[:3]),  yerr=np.std(dat[:3])/np.sqrt(3),  color = 'purple', error_kw=cap)
#    axes[k].bar([1], np.mean(dat[3:6]), yerr=np.std(dat[3:6])/np.sqrt(3), color = 'red', error_kw=cap)
#    axes[k].bar([2], np.mean(dat[6:]),  yerr=np.std(dat[6:])/np.sqrt(9),  color = 'green', error_kw=cap)
#    
#    axes[k].set_xlim(-0.2, 3)
#    axes[k].set_ylim(0, 0.82)
#    axes[k].set_xticks([0.4, 1.4, 2.4])
#    axes[k].set_xticklabels(['i', 'ii', 'iii'], fontsize = 16)
#    axes[k].set_xlabel(groups[k], fontsize = 16)
#axes[0].set_ylabel('Locomotion fitness', fontsize = 16)
##################################################################
#ax1.imshow(plt.imread('../networks/Network_9.png'))
#ax2.imshow(plt.imread('../networks/Network_23.png'))
#ax3.imshow(plt.imread('../networks/Network_45a.png'))
#ax4.imshow(plt.imread('../networks/Network_45b.png'))
#for ax in axes2: ax.axis('Off')
##################################################################
##################################################################
#fig.subplots_adjust(left=0.06, bottom=0.11, right=0.99, top=0.98)
#plt.savefig('resumen.png', dpi = 300)
