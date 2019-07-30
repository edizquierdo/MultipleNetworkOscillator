import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
nrods = 51

import time
t = time.time()

os.chdir('./save')
for ind in [30, 101, 103]:
        os.system('cp ../../ind/%d/best.gen.dat ./'%(ind))
        os.system('./main %d'%ind)
os.chdir('../')



