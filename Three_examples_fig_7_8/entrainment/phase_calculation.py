import numpy as np
import os
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec

[AS, DA, DB, DD,  VD, VB, VA] = [0, 1, 2, 3, 4, 5, 6]
#######################################################################
plt.close('all')
fig = plt.figure(figsize = [12,6])
gs = gridspec.GridSpec(140, 210)
ax = plt.subplot(gs[:,:])

######################################################################
import time
tt = time.time()
######################################################################
dt = 0.001
T=14
t = np.arange(0, T+dt, dt)
######################################################################
wrap = lambda x: ((x + np.pi)%(2*np.pi ) - np.pi)
angw = lambda x: wrap(x) if wrap(x) > 0 else wrap(x) + 2*np.pi
######################################################################
######################################################################
def reconstruction(trace):
    dtra = np.diff(trace)
    u, d, aux = np.ones(7), np.ones(7), 0
    for i in range(7):
        u[i] = int(aux + np.argmax(dtra[aux:aux+300]))
        aux = int(u[i] + np.argmax(-dtra[u[i]:u[i]+300]))
        d[i] = aux
    period = np.mean(np.diff([u,d])) #ms
    phase = angw((np.mean((u+d)/2%period)/period)*2*np.pi)
    return (phase)

for ind in [30, 101, 103]:
    for entrain in ['tail', 'head']:
        print ('ind = %d, entrain = %s'%(ind, entrain))
        act_reference = np.loadtxt('entrain_%s/play/act_played%d.dat'%(entrain, ind))
        phase_shift = np.zeros(6)
        for rep in range(6):
            act_evaluated = np.loadtxt('entrain_%s/play/act_%d_%d.dat'%(entrain, ind, rep))
            r = reconstruction(act_reference[:,DB])
            e = reconstruction(act_evaluated[:,DB])
            phase_shift[rep] = angw(wrap(r - e))*180/np.pi
            print(angw(wrap(r - e))*180/np.pi)
            np.savetxt('figures/data_entrain_%s_%d'%(entrain, ind), phase_shift)
######################################################################

#anterio30 = np.loadtxt('figures/data_entrain_head_30')
#anterio101 = np.loadtxt('figures/data_entrain_head_101')
#anterio103 = np.loadtxt('figures/data_entrain_head_103')

#posterio30 = np.loadtxt('figures/data_entrain_tail_30')
#posterio101 = np.loadtxt('figures/data_entrain_tail_101')
#posterio103 = np.loadtxt('figures/data_entrain_tail_103')


#######################################################################
#########  To Check phase calculation function
#######################################################################
#######################################################################
#plt.close('all')
#fig = plt.figure(figsize = [12,6])
#gs = gridspec.GridSpec(140, 210)
#ax = plt.subplot(gs[:,:])
#######################################################################
#wrap = lambda x: ((x + np.pi)%(2*np.pi ) - np.pi)
#angw = lambda x: wrap(x) if wrap(x) > 0 else wrap(x) + 2*np.pi
#######################################################################
#######################################################################
#def reconstruction(trace):
#    dtra = np.diff(trace)
#    ax.plot(dtra)
#    u, d, aux = np.ones(7), np.ones(7), 0
#    for i in range(7):
#        u[i] = int(aux + np.argmax(dtra[aux:aux+300]))
#        aux = int(u[i] + np.argmax(-dtra[u[i]:u[i]+300]))
#        d[i] = aux
#    plt.plot(u, np.zeros(7)+0.05, '.r')
#    plt.plot(d, np.zeros(7), '.k')
#    period = np.mean(np.diff([u,d])) #ms
#    phase = angw((np.mean((u+d)/2%period)/period)*2*np.pi)
#    return (phase)

#act = np.loadtxt('entrain_head/play/act_played30.dat')
#reconstruction(act[:,3])
