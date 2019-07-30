import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os

sel = np.loadtxt('list_selected')

tr = np.loadtxt('kinogram_wavelength')

for s in sel:
    os.system('cp wavelength_figures/plot_sel_%d.png ./filtered_figures/'%s)
    
a = np.array([tr[s] for s in sel])

plt.plot(a, '.')

plt.ylabel('Wavelength')
plt.ylim(0.38, 0.92)
plt.axhline(y = 0.4, color = 'r', linestyle = '--')
plt.axhline(y = 0.9, color = 'r', linestyle = '--')

plt.xlabel('Individuals')


plt.savefig('wavelength_filtered.png')


