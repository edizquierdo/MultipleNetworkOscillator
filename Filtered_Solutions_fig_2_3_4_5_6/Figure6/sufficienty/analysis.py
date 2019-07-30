import numpy as np
import os


selected = [4, 5, 9, 21, 23, 30, 37, 42, 45, 53, 63, 66, 68, 101, 103]
######################################################################
#######################################################################

#minimals = ['full','minimal_9', 'minimal_23', 'minimal_45', 'chemical1', 'chemical2']
#for evaluated in minimals:
#    data = np.zeros((15, 10))
#    for k, ind in enumerate(selected):
#        os.system('cp scripts/%s ind/%d'%(evaluated, ind))
#        os.chdir('ind/%d'%ind)
#        os.system('./%s'%evaluated)
#        data[k] = np.loadtxt('fitness.dat')
#        os.chdir('../../')
#        print(ind)
#    np.savetxt('data/data_%s'%evaluated, data)
#######################################################################

full = np.loadtxt('data/data_full')
groups = ['minimal_9', 'minimal_45', 'minimal_23', 'chemical1', 'chemical2']
data = np.array([np.mean(np.loadtxt('data/data_%s'%s), axis = 1)/np.mean(full, axis = 1) for s in groups])
data[np.where(data<0)] = 0 
np.savetxt('data_suff.dat', data.T)

#connections = ['AS_VAn', 'DA_ASn', 'VB_DBn', 'DB_DDn', 'VAn_DD']

#max_fit_ablated = np.zeros((15, 5))
#for k, ind in enumerate(selected):
#    for j, con in enumerate(connections):
#        f = np.loadtxt('./selected/%d/fitness_%s.dat'%(ind, con))
#        mean = np.mean(f, axis = 1)
#        max_fit_ablated[k, j] = np.max(mean)
#np.savetxt('data_ablation.dat', max_fit_ablated)
