import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os

best = np.zeros((160, 2002))
pop = np.zeros((160, 2002))


#Read data
for i in range(160):
    f = np.loadtxt('../evolution_results/%d/fitness.dat'%(i+1))
    ngen = len(f)
    if (ngen == 2002):
        best[i] = f[:,1]
        pop[i] = f[:,2]
    else:
        best[i][:ngen] = f[:,1]
        pop[i][:ngen] = f[:,2]
        for ne in range(ngen, 2002):
            best[i][ne] = f[-1:,1]
            pop[i][ne] = f[-1,2]
        print(i)

#best = np.loadtxt('best.txt')
#pop = np.loadtxt('pop.txt')

bstg1 = np.zeros((160, 300))
pstg1 = np.zeros((160, 300))


bstg2 = np.zeros((160, 2002))
pstg2 = np.zeros((160, 2002))

fstg = np.where(np.diff(best, axis = 1) < 0)[1]

for i in range(160):
    ch = np.min([299, fstg[i] ])
    bstg1[i][:ch] = best[i][:ch]
    bstg1[i][ch:] = best[i][ch]
    pstg1[i][:ch] = pop[i][:ch]
    pstg1[i][ch:] = pop[i][ch]
    
    last = 2002 - (ch)
    bstg2[i][:last] = best[i][ch:]
    bstg2[i][last:] = best[i][2001]
    plt.plot(bstg2[i])
    pstg2[i][:ch] = pop[i][:ch]
    pstg2[i][ch:] = pop[i][ch]
    
ms1, stds1 = np.mean(bstg1, axis = 0), np.std(bstg1, axis = 0)
pms1, pstds1 = np.mean(pstg1, axis = 0), np.std(pstg1, axis = 0)

ms2, stds2 = np.mean(bstg2, axis = 0), np.std(bstg2, axis = 0)
pms2, pstds2 = np.mean(pstg2, axis = 0), np.std(pstg2, axis = 0)

last_gen = bstg2[:,-1]
print('a')
plt.close('all')
print('a')
fig = plt.figure(figsize = [12,4])

colors = ['#2689ad','#d95f02','#7570b3']

gm = gridspec.GridSpec(20, 1000)
ax0 = plt.subplot(gm[1:-1, 1:250])
ax1 = plt.subplot(gm[1:-1, 270:740])
ax2 = plt.subplot(gm[1:-1, 770:-10])

ax0.plot(np.arange(300), ms1, color= colors[0])
ax0.fill_between(np.arange(300), ms1 -stds1, ms1 + stds1, alpha=0.3, facecolor= colors[0])

ax0.plot(np.arange(300), pms1, color= colors[1])
ax0.fill_between(np.arange(300), pms1 -pstds1, pms1 + pstds1, alpha=0.3, facecolor= colors[1])


ax1.plot(np.arange(2002), ms2, color= colors[0])
ax1.fill_between(np.arange(2002), ms2 -stds2, ms2 + stds2, alpha=0.3, facecolor= colors[0])

ax1.plot(np.arange(2002), pms2, color= colors[1])
ax1.fill_between(np.arange(2002), pms2 -pstds2, pms2 + pstds2, alpha=0.3, facecolor= colors[1])

ax0.text(195, 0.7, 'Best', fontsize = 15, color = colors[0], ha = 'center', va = 'center')
ax0.text(195, 0.3, 'Population\nAverage', fontsize = 14, color = colors[1], ha = 'center', va = 'center')

for ax in [ax0, ax1]:
    ax.set_ylim(0,1)
    ax.set_xlabel('Generation', fontsize = 18, labelpad = 2)

ax0.set_xlim(0, 300)
ax1.set_xlim(0, 2000)
ax0.set_yticks([0.0, 0.5, 1.0])
ax1.set_yticks([])
ax0.set_ylabel('Fitness', fontsize = 18, labelpad = 3)
ax0.set_xticks([50, 150, 250])
ax1.set_xticks([400, 800, 1200, 1600])

h = ax2.hist(last_gen, bins=20, orientation="vertical", color = colors[2])
ax2.set_xlim(0,1)
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
ax2.set_ylabel('Counts', fontsize = 18, rotation = 270, labelpad = 13)
#ax2.axvline(x = 0.8, linestyle = '--', color = 'red', linewidth = 2)
ax2.set_ylim(0, 110)
ax2.set_xlabel('Fitness', fontsize = 18, labelpad = 2)
#ax0.set_xticks([])
ax2.set_yticks([25, 50, 75, 100])
ax2.get_children()[19].set_facecolor('r')


plt.figtext(0.1666, 0.93, 'Stage 1', color = 'k', fontsize = 20, ha = 'center')
plt.figtext(0.5, 0.93, 'Stage 2', color = 'k', fontsize = 20, ha = 'center')
plt.figtext(0.85, 0.93, 'Ensemble', color = 'k', fontsize = 20, ha = 'center')


fig.subplots_adjust(left = 0.05, bottom = 0.13, right = 0.96, top = 0.93)

plt.savefig('fitness.png', dpi = 200)
########################################
### Selection

sel = np.where(last_gen > 0.95)[0]
os.system('mkdir ../selected')
k = 0
for i in np.delete(sel, [41, 95]): ## this solutions idd not finis the evolution, si they do not save a best gen
    os.system('cp -r ../evolution_results/%d ../selected/%d'%(i+1, k))
    k +=1








