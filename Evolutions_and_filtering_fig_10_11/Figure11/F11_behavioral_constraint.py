import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
from matplotlib_venn import venn3
nsel = 104

wavelength = np.loadtxt('Wavelength/kinogram_wavelength')
curvature = (1./np.loadtxt('Trayectory_curv/trayectory_curvature'))/1000.0
speed_ablated_B = np.loadtxt('Ablations/B_class_NMJ_ablation/speed')/0.00022
speed_ablated_A = np.loadtxt('Ablations/A_class_NMJ_ablation/speed')/0.00022
curv_profile = np.loadtxt('./curvature_profile/slope_curv_prof')

plt.close('all')
fig = plt.figure(figsize = [16,16])
gm = gridspec.GridSpec(110, 105)
ax1 = plt.subplot(gm[5:45, :40])
ax1t = plt.subplot(gm[:5, :40])
ax1r = plt.subplot(gm[5:45, 40:45])
ax2 = plt.subplot(gm[5:45, 60:100])
ax2t = plt.subplot(gm[:5, 60:100])
ax2r = plt.subplot(gm[5:45, 100:105])
ax3 = plt.subplot(gm[75:90, :40])
ax4 = plt.subplot(gm[60:, 50:])
colors = ['red', 'green', 'blue']
fz = 18

#### PLot ablation analysis
ax1.plot(speed_ablated_B, speed_ablated_A, '.k')

ax1.set_xlabel('Normalize speed (B-class ablated)', fontsize = fz)
ax1.set_xlim(0, 1.1)
ax1.axvspan(0.0, 0.2, facecolor = colors[0], alpha = 0.25, zorder = 0)

ax1.set_ylabel('Normalize speed (A-class ablated)', fontsize = fz)
ax1.set_ylim(0, 1.1)
ax1.axhspan(0.2, 1.1, facecolor = colors[0], alpha = 0.25, zorder = 0)

ax1t.hist(speed_ablated_B, bins=22, range = [0, 1.1], orientation="vertical", color = 'k')
[ax1t.get_children()[b+2].set_facecolor(colors[0]) for b in range(4)]
ax1t.set_xlim(0,1.1)
ax1t.set_ylim(0,71)

ax1r.hist(speed_ablated_A, bins=22, range = [0, 1.1], orientation="horizontal", color = 'k')
[ax1r.get_children()[b+2].set_facecolor(colors[0]) for b in range(4, 22)]
ax1r.set_ylim(0,1.1)

#### Plot body posture analysis
wavelength_threshold = 0.9
ax2.plot(wavelength, curv_profile, '.k')

ax2.set_xlabel('Wavelength', fontsize = fz)
ax2.set_xlim(0, 1.5)
ax2.axvspan(0.4, wavelength_threshold, facecolor = colors[1], alpha = 0.25, zorder = 0)

ax2.set_ylabel('Body curvature', fontsize = fz)
ax2.set_ylim(-0.02, 0.02)
ax2.axhspan(-0.02, 0.0, facecolor = colors[1], alpha = 0.25, zorder = 0)
ax2.set_yticks([-0.02, -0.01, 0.00, 0.01, 0.02])


ax2t.hist(wavelength, bins=30, range = [0, 1.5], orientation="vertical", color = '0.1')
[ax2t.get_children()[b+2].set_facecolor(colors[1]) for b in range(10, 18)]
ax2t.set_xlim(0,1.5)

ax2r.hist(curv_profile, bins=20, range = [-0.02, 0.02], orientation="horizontal", color = '0.1')
[ax2r.get_children()[b+2].set_facecolor(colors[1]) for b in range(10)]
ax2r.set_ylim(-0.02,0.02)

[ax.axis('off') for ax in [ax2t, ax1t, ax2r, ax1r]]

#### Histogram Trayectory curvature
h = ax3.hist(curvature, bins=21, range = [0, 2.1], orientation="vertical", color = '0.1')
ax3.set_xlim(0,2.1)
ax3.set_ylabel('Counts', fontsize = 18, rotation = 90, labelpad = 13)
ax3.set_ylim(0, 15)
ax3.set_xlabel('Curvature (1/mm)', fontsize = 18, labelpad = 2)
[ax3.get_children()[b].set_facecolor(colors[2]) for b in range(10)]
ax3.axvspan(0, 1, facecolor = colors[2], alpha = 0.25, zorder = 0)

###################
#### Venn diagram
#check solutions that fulfill each of the criterias
check_body = np.where((wavelength < wavelength_threshold) & (curv_profile <= 0))[0] # all of them are greater than 0.4
check_trayectory = np.where(curvature < 1)[0]
check_ablations = np.where((speed_ablated_B < 0.2) & (speed_ablated_A > 0.2))[0]

a = np.intersect1d(np.intersect1d(check_body,check_trayectory), check_ablations)
all_crit = len(a)

none = len(np.setdiff1d(np.arange(104), np.union1d(check_body, np.union1d(check_trayectory,check_ablations))))

only_b = len(np.setdiff1d(check_body, np.union1d(check_trayectory,check_ablations)))
only_t = len(np.setdiff1d(check_trayectory, np.union1d(check_body,check_ablations)))
only_a = len(np.setdiff1d(check_ablations, np.union1d(check_body,check_trayectory)))

a_t_not_b = len(np.setdiff1d(np.intersect1d(check_trayectory,check_ablations), check_body))
a_b_not_t = len(np.setdiff1d(np.intersect1d(check_ablations,check_body), check_trayectory))
b_t_not_a = len(np.setdiff1d(np.intersect1d(check_trayectory,check_body), check_ablations))

sets = (only_a, only_b, a_b_not_t, only_t, a_t_not_b, b_t_not_a, all_crit)
none = nsel - (only_a + only_b + a_b_not_t + only_t + a_t_not_b + b_t_not_a + all_crit)
###################
c = venn3(subsets = sets, set_labels = ('Ablations', 'Posture', 'Trayectory'), ax = ax4, alpha = 0.65)
for text in c.set_labels:
    text.set_fontsize(18)
for text in c.subset_labels:
    text.set_fontsize(18)
    

plt.figtext(0.09, 0.92, 'A and B class neuron ablation analysis', fontsize = 22, bbox={'facecolor':colors[0], 'alpha':0.5, 'pad':18})
plt.figtext(0.62, 0.92, 'Body posture analysis', fontsize = 22, bbox={'facecolor':colors[1], 'alpha':0.5, 'pad':18})
plt.figtext(0.132, 0.4, 'Trayectory curvature analysis', fontsize = 22, bbox={'facecolor':colors[2], 'alpha':0.5, 'pad':18})

plt.figtext(0.57, 0.13, 'None = %d'%none, color = 'k', fontsize = 18, ha = 'center')

#fig.subplots_adjust(left=0.05, bottom=0.02, right=0.98, top=0.96)
plt.savefig('behavioral_analisys.png')
plt.savefig('behavioral_analisys.pdf')

np.savetxt('list_selected', a, fmt = '%d')
