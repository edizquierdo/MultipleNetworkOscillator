import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
nrods = 51

sel = [4, 5, 9, 21, 23, 30, 37, 42, 45, 53, 63, 66, 68, 101, 103]

from mpi4py import MPI
comm=MPI.COMM_WORLD
numproc=comm.size
rank = comm.Get_rank()


for k, ind in enumerate(sel):#[23]:#
    if k%numproc == comm.rank:
        try:
            os.system('cp one/main filtered_agents/%d'%ind)
            os.chdir('filtered_agents/%d'%ind)
            os.system('./main')
            print('simulating individual %d'%ind)
            os.chdir('../../')
        except:
            print('\nsimulation %d fail\n'%ind)

