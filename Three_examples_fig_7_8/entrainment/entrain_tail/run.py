import os

import time
t = time.time()

os.chdir('save')
os.system('make clean')
os.system('make')
os.chdir('../play')
os.system('make clean')
os.system('make')
os.chdir('../')
os.system('python save_traces.py')
os.system('python play_traces.py')

print (time.time() - t)
