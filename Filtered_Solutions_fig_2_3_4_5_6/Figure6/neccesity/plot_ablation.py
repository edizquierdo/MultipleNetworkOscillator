import numpy as np
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import os

selected = [4, 5, 9, 21, 23, 30, 37, 42, 45, 53, 63, 66, 68, 101, 103]
##################################################################
#connections = ['AS_VAn', 'DA_ASn', 'VB_DBn', 'DB_DDn', 'VAn_DD']

#max_fit_ablated = np.zeros((15, 5))
#for k, ind in enumerate(selected):
#    for j, con in enumerate(connections):
#        f = np.loadtxt('./selected/%d/fitness_%s.dat'%(ind, con))
#        mean = np.mean(f, axis = 1)
#        max_fit_ablated[k, j] = np.max(mean)
#np.savetxt('data_ablation.dat', max_fit_ablated)
##################################################################
labels_syn = ['AS'+ r'$\leftrightarrow$' +'VA' + r'$_{next}$'\
, 'DA'+ r'$\leftrightarrow$' +'AS' + r'$_{next}$'\
, 'DA'+ r'$\leftrightarrow$' +'AS' + r'$_{next}$'\
, 'DB'+ r'$\to$' +'DD' + r'$_{next}$'\
, 'VA' + r'$_{next}$'+ r'$\to$' +'DD']

max_fit_ablated = np.loadtxt('data_ablation.dat')

## Each individual intraunit
plt.close('all')
fig = plt.figure(figsize = [15, 6])
gm = gridspec.GridSpec(150, 140)
axi = plt.subplot(gm[15:-15, 0:130])
axc = plt.subplot(gm[15:-15, 132:135])
#########################

axi.set_xlabel('Individuals', fontsize = 17, labelpad = 14)
axi.xaxis.set_label_position('top')
axi.set_yticks(np.arange(3))
colores = ['purple', 'red', 'green', 'k', 'k']
for y in np.arange(5):
    axi.text(-1.5, y, labels_syn[y] + '\n ', fontsize = 17, color = colores[y], va = 'center', ha = 'center')
axi.set_yticklabels([])
axi.set_xticks(np.arange(15))
axi.set_xticklabels([str(j) for j in selected], fontsize = 14)
axi.xaxis.tick_top()
axi.set_ylim(-0.51, 4.51)
axi.set_xlim(-0.51, 14.51)
for xo in np.arange(16)-0.5:
    axi.axvline(x = xo, linewidth = 0.4, linestyle = '-', color = 'k')
for yo in np.arange(4)-0.5:
    axi.axhline(y = yo, linewidth = 0.4, linestyle = '-', color = 'k')
axi.spines['top'].set_visible(False)
axi.spines['right'].set_visible(False)
axi.spines['bottom'].set_visible(False)
axi.spines['left'].set_visible(False)
im = axi.imshow(max_fit_ablated.T, interpolation = 'none', cmap=plt.get_cmap('Purples'), vmin = 0, vmax = 1, aspect = 'auto')
barim = fig.colorbar(im, ax = axi, cax = axc)

###### Letters
#plt.figtext(0.012, 0.95, 'A', fontsize = 40, ha = 'center', va = 'center')
#plt.figtext(0.345, 0.95, 'B', fontsize = 40, ha = 'center', va = 'center')

fig.subplots_adjust(left=0.1, bottom=0.00, right=0.96, top=0.92)
axc.set_ylabel('Locomotion fitness', labelpad = 22, rotation = -90, fontsize = 17)
#axc.set_ylabelpad(15)

plt.savefig('interunit_neccesity.png', dpi = 300)
#plt.savefig('interunit_neccesity.pdf')
#plt.savefig('../../Figures/interunit_neccesity.pdf')
plt.savefig('../../Figures/interunit_neccesity.png', dpi = 300)


