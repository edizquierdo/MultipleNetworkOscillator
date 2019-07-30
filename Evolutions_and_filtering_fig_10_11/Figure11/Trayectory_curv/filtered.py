import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os

sel = np.loadtxt('list_selected')

tr = np.loadtxt('trayectory_curvature')

for s in sel:
    os.system('cp trayectory_figures/plot_sel_%d.png ./filtered/'%s)
    
a = np.array([1./tr[s] for s in sel])

plt.plot(a, '.')
plt.savefig('filtered.png')

