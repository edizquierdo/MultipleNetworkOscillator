import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os

connections = ['AS_VAn', 'DA_ASn', 'VB_DBn', 'DB_DDn', 'VAn_DD']

for con in connections:
    os.chdir('./%s'%con)
    os.system('make clean')
    os.system('make')
    os.chdir('../')
