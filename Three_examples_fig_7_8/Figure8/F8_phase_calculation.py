import numpy as np
import os
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
from aux import *

[AS, DA, DB, DD,  VD, VB, VA] = [0, 1, 2, 3, 4, 5, 6]
######################################################################
pal=["r", "c", "y", 'b', 'm']
## colors according to figure necc and suff
new_order = [2, 1, 0, 3, 4]
pal = [pal[i] for i in new_order]
######################################################################
import time
tt = time.time()
######################################################################
dt = 0.001
T=14
t = np.arange(0, T+dt, dt)
######################################################################
wrap = lambda x: ((x + np.pi)%(2*np.pi ) - np.pi)
#angw = lambda x: wrap(x) if wrap(x) > 0 else wrap(x) + 2*np.pi
angw = lambda x: wrap(x) if wrap(x) > 0 else wrap(x) + 2*np.pi
######################################################################
######################################################################
def reconstruction(trace):
    dtra = np.diff(trace)
    u, d, aux = np.ones(4), np.ones(4), 0
    for i in range(4):
        u[i] = int(aux + np.argmax(dtra[aux:aux+3000]))
        aux = int(u[i] + np.argmax(-dtra[u[i]:u[i]+3000]))
        d[i] = aux
    period = np.mean(np.diff([u,d])) #ms
    phase = angw((np.mean((u+d)/2%period)/period)*2*np.pi)
    return (phase)
######################################################################
##############  Simulations and hase calculations ####################
#######################################################################
#ind = 101
#os.system('cp script/%s ind/%d'%('VB_DB', ind))
#os.chdir('ind/%d'%ind)
#os.system('./VB_DB')
#act = np.loadtxt('act.dat').T
#os.chdir('../../')
#seg3 = act[22:29]
#seg4 = act[29:36]
#np.savetxt('./data/traces_seg_1_VB_DB', seg3)
#np.savetxt('./data/traces_seg_2_VB_DB', seg4)

#a1 = reconstruction(seg3[DB])
#a2 = reconstruction(seg3[AS])
#a3 = reconstruction(seg3[VD])
#a4 = reconstruction(seg3[VB])
#a5 = reconstruction(seg4[DB])
#a_VB_DB = [a1, a2, a3, a4, a5]
#np.savetxt('./data/a_VB_DB', a_VB_DB)
#######################################################################
#ind = 103
#os.system('cp script/%s ind/%d'%('DA_AS', ind))
#os.chdir('ind/%d'%ind)
#os.system('./DA_AS')
#act = np.loadtxt('act.dat').T
#os.chdir('../../')
#seg3 = act[22:29]
#seg4 = act[29:36]
#np.savetxt('./data/traces_seg_1_DA_AS', seg3)
#np.savetxt('./data/traces_seg_2_DA_AS', seg4)

#a1 = reconstruction(seg3[DB])
#a2 = reconstruction(seg3[AS])
#a3 = reconstruction(seg3[VD])
#a4 = reconstruction(seg3[VB])
#a5 = reconstruction(seg3[DA])
#a6 = reconstruction(seg4[AS])
#a7 = reconstruction(seg4[DA])
#a8 = reconstruction(seg4[DB])
#a_DA_AS = [a1, a2, a3, a4, a5, a6, a7, a8]
#np.savetxt('./data/a_DA_AS', a_DA_AS)
#######################################################################
#ind = 30
#os.system('cp script/%s ind/%d'%('AS_VA', ind))
#os.chdir('ind/%d'%ind)
#os.system('./AS_VA')
#act = np.loadtxt('act.dat').T
#os.chdir('../../')
#seg3 = act[22:29]
#seg4 = act[29:36]
#np.savetxt('./data/traces_seg_1_AS_VA', seg3)
#np.savetxt('./data/traces_seg_2_AS_VA', seg4)

#a1 = reconstruction(seg4[DB])
#a2 = reconstruction(seg4[AS])
#a3 = reconstruction(seg4[VD])
#a4 = reconstruction(seg4[VA])
#a5 = reconstruction(seg4[VB])
#a6 = reconstruction(seg3[AS])
#a7 = reconstruction(seg3[DA])
#a8 = reconstruction(seg3[DB])
#a_AS_VA = [a1, a2, a3, a4, a5, a6, a7, a8]
#np.savetxt('./data/a_AS_VA', a_AS_VA)
######################################################################
#####
######################################################################
##############################################################
plt.close('all')
fig = plt.figure(figsize = [40,9])
#axes
gm = gridspec.GridSpec(23, 290)
axnet1 = plt.subplot(gm[2:23, :90])
axnet2 = plt.subplot(gm[2:23, 100:190])
axnet3 = plt.subplot(gm[2:23, 200:290])
################################################################## 
############### Networks ###########################################
#########################   VB_DB  #############################################
ind = 101
phen = np.loadtxt('../network_kinograms/ind/%d/phenotype.dat'%ind)[[21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 39, 40, 41, 42, 43]]
neurons = [0, 2, 4, 5, 9]
synapses = [1, 3, 5, 14]
plot_segments(axnet1, pal[0], neurons, synapses, phen)
#########################   DA_AS  #############################################
ind = 103
phen = np.loadtxt('../network_kinograms/ind/%d/phenotype.dat'%ind)[[21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 39, 40, 41, 42, 43]]
neurons = [2, 0, 1, 7, 8, 9]
synapses = [3, 0, 13, 6, 8]
plot_segments(axnet2, pal[1], neurons, synapses, phen)
#########################   AS_VA  #############################################
ind = 30
phen = np.loadtxt('../network_kinograms/ind/%d/phenotype.dat'%ind)[[21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 39, 40, 41, 42, 43]]
neurons = [9, 7, 11, 13, 0, 1, 2]
synapses = [9, 7, 10, 12, 0, 2]
plot_segments(axnet3, pal[2], neurons, synapses, phen)
################################################################## 
################################################################## 
######################################################################
fz = 32
lstyle = "|-|, widthA=0.5, widthB=0.5"
lstyle_inter = 'simple, tail_width = %s, head_width = %s, head_length = %s'%(1.1, 2.1, 2.1)
######################################################################
a_VB_DB = np.loadtxt('./data/a_VB_DB')
#n_VB_DB = ['DB', 'AS', 'VD', 'VB', 'DBn']
[DB, AS, VD, VB, DBn] = [0, 1, 2, 3, 4]
axnet1.text(7, 16, '%d'%(angw(a_VB_DB[AS]-a_VB_DB[DB])*180/np.pi) + r'$^{o}$', fontsize = fz)
axnet1.text(4.2, 13, '%d'%(angw(a_VB_DB[VD]-a_VB_DB[AS])*180/np.pi) + r'$^{o}$', fontsize = fz)
axnet1.text(6, 5, '%d'%(angw(a_VB_DB[VB]-a_VB_DB[VD])*180/np.pi) + r'$^{o}$', fontsize = fz)
axnet1.text(16, 9, '%d'%(angw(a_VB_DB[DBn]-a_VB_DB[VB])*180/np.pi) + r'$^{o}$', fontsize = fz)

yloc = 1.01
axnet1.text(0.5, 1.05, r'$\Delta\phi$'+' = %d'%np.rint(angw(a_VB_DB[DBn]-a_VB_DB[DB])*180/np.pi) + r'$^{o}$', transform=axnet1.transAxes, fontsize = fz)
axnet1.annotate('', xy=(0.95, yloc), xycoords='axes fraction', xytext=(0.3, yloc),arrowprops=dict(arrowstyle=lstyle_inter, color='k', linewidth = 1))
######################################################################
a_DA_AS = np.loadtxt('./data/a_DA_AS')
#n_DA_AS = ['DB', 'AS', 'VD', 'VB', 'DA', 'ASn', 'DAn', 'DBn']
[DB, AS, DA, ASn, DAn, DBn] = [0, 1, 4, 5, 6, 7]
axnet2.text(7, 16, '%d'%(angw(a_DA_AS[AS]-a_DA_AS[DB])*180/np.pi) + r'$^{o}$', fontsize = fz)
axnet2.text(1, 16, '%d'%(angw(a_DA_AS[DA]-a_DA_AS[AS])*180/np.pi) + r'$^{o}$', fontsize = fz)
axnet2.text(15, 17, '%d'%(angw(a_DA_AS[ASn]-a_DA_AS[DA])*180/np.pi) + r'$^{o}$', fontsize = fz)
axnet2.text(21, 16, '%d'%(angw(a_DA_AS[DAn]-a_DA_AS[ASn])*180/np.pi) + r'$^{o}$', fontsize = fz)
axnet2.text(24, 14, '%d'%(angw(a_DA_AS[DBn]-a_DA_AS[DAn])*180/np.pi) + r'$^{o}$', fontsize = fz)

yloc = 1.01
axnet2.text(0.5, 1.05, r'$\Delta\phi$'+' = %d'%np.rint(angw(a_DA_AS[DBn]-a_DA_AS[DB])*180/np.pi) + r'$^{o}$', transform=axnet2.transAxes, fontsize = fz)
axnet2.annotate('', xy=(0.95, yloc), xycoords='axes fraction', xytext=(0.3, yloc),arrowprops=dict(arrowstyle=lstyle_inter, color='k', linewidth = 1))
######################################################################
a_AS_VA = np.loadtxt('./data/a_AS_VA')
#n_AS_VA = ['DBn', 'ASn', 'VDn', 'VAn', 'VBn', 'AS', 'DA', 'DB']
[DBn, ASn, VDn, VAn, AS, DA, DB] = [0, 1, 2, 3, 5, 6, 7]
axnet3.text(27, 16, '%d'%(angw(a_AS_VA[ASn]-a_AS_VA[DBn])*180/np.pi) + r'$^{o}$', fontsize = fz)
axnet3.text(24, 12, '%d'%(angw(a_AS_VA[VDn]-a_AS_VA[ASn])*180/np.pi) + r'$^{o}$', fontsize = fz)
axnet3.text(22.5, 4, '%d'%(angw(a_AS_VA[VAn]-a_AS_VA[VDn])*180/np.pi) + r'$^{o}$', fontsize = fz)


axnet3.text(15.5, 10, '%d'%(angw(a_AS_VA[AS]-a_AS_VA[VAn])*180/np.pi) + r'$^{o}$', fontsize = fz)

axnet3.text(0, 16, '%d'%(angw(a_AS_VA[DA]-a_AS_VA[AS])*180/np.pi) + r'$^{o}$', fontsize = fz)
axnet3.text(4, 14, '%d'%np.rint(angw(a_AS_VA[DB]-a_AS_VA[DA])*180/np.pi) + r'$^{o}$', fontsize = fz)

yloc = 1.01
axnet3.text(0.5, 1.05, r'$\Delta\phi$'+' = %d'%np.rint(angw(a_AS_VA[DBn]-a_AS_VA[DB])*180/np.pi) + r'$^{o}$', transform=axnet3.transAxes, fontsize = fz)
axnet3.annotate('', xy=(0.3, yloc), xycoords='axes fraction', xytext=(0.95, yloc),arrowprops=dict(arrowstyle=lstyle_inter, color='k', linewidth = 1))
################################################################## 
######################################################################
xpos = 0.015, 0.340, 0.675
ypos = [0.94, 0.43]
fzl = 70

letters = ['A', 'B', 'C']
[plt.figtext(xpos[i], ypos[0], letters[i], fontsize = fzl, ha = 'center', va = 'center') for i in range(3)]
#[plt.figtext(xpos[i], ypos[1], 'B%d'%(i+1), fontsize = fzl, ha = 'center', va = 'center') for i in range(3)]
######################################################################

fig.subplots_adjust(left=0.02, bottom=0.00, right=0.98, top=1.02)
plt.savefig('phase_delay.png'%(), dpi = 100)
plt.savefig('phase_delay.pdf')
#plt.savefig('phase_delay.svg')

