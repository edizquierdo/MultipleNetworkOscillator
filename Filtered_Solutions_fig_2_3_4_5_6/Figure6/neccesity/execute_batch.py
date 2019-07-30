import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
nrods = 51

from mpi4py import MPI
comm=MPI.COMM_WORLD
numproc=comm.size
rank = comm.Get_rank()

#[101, 21, 4, 9, 5, 66, 53, 45, 103, 30, 68, 37, 63, 23, 42]
connections = ['AS_VAn', 'DA_ASn', 'VB_DBn', 'DB_DDn', 'VAn_DD']
selected = [4, 5, 9, 21, 23, 30, 37, 42, 45, 53, 63, 66, 68, 101, 103]

for k, ind in enumerate(selected):
    if k%numproc == comm.rank:
        os.chdir('./selected/%d'%ind)
        for j, con in enumerate(connections):
            fitness = np.zeros((441, 8))
            os.system('cp ../../code_folders/%s/main ./'%(con))
            for rep in np.arange(441):
                try:
                    os.system('./main %d'%rep)
                    f = np.loadtxt('./fitness.dat')
                    fitness[rep] = f
                except:
                    print('something wrong in ind %d, con %s, rep %d'%(ind, con, rep))
            np.savetxt('fitness_%s.dat'%con, fitness)

