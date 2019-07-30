import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os

###################
####################
plt.close('all')
fig = plt.figure(figsize = [16,7])
gm = gridspec.GridSpec(80, 300)
ax1 = plt.subplot(gm[:50, :78])
ax11 = plt.subplot(gm[:50, 97:175])

ax2 = plt.subplot(gm[60:77, :78])
ax22 = plt.subplot(gm[60:77, 97:175])
####################
axo1 = plt.subplot(gm[:52, 200:288])
axoc = plt.subplot(gm[:52, 291:294])
axo2 = plt.subplot(gm[60:77, 200:288])
#########################################################################################

selected, wt, ns, bwm, wt2, ns2, bwm2, curv, curv_prof, speeds_mean = np.load('myfile.npy')

#######################################################################
x1, x2, x3, y1, y2 = 0.015, 0.31, 0.61, 0.96, 0.30
letters = [[x1, y1, 'A'], [x1, y2, 'B'], [x2, y1, 'C'], [x2, y2, 'D'], [x3, y1, 'E'], [x3, y2, 'F']]
for l in letters:
    plt.figtext(l[0], l[1], l[2], fontsize = 34, ha = 'center', va = 'center')
#########################################################################################
#######################################################################
control, headparal, midbodyparal, edge = '#1b9e77','#d95f02','#7570b3', 'k'
[head, neck, tail] = [0, 1, 2]

ax1.scatter(1*np.ones(15), wt[:,head], lw = 0.5, s = 100, color = control, edgecolors=edge, label = 'Control') #wild type head amplitude
ax1.scatter(2*np.ones(15), ns[:,head], lw = 0.5, s = 100, color = headparal, edgecolors=edge, label = 'Head paralyzed') #wild type head amplitude
ax1.scatter(4*np.ones(15), wt[:,tail], lw = 0.5, s = 100, color = control, edgecolors=edge) #wild type head amplitude
ax1.scatter(5*np.ones(15), ns[:,tail], lw = 0.5, s = 100, color = headparal, edgecolors=edge) #wild type head amplitude

ax1.scatter(8*np.ones(15), wt[:,head], lw = 0.5, s = 100, color = control, edgecolors=edge) #wild type head amplitude
ax1.scatter(9*np.ones(15), bwm[:,head], lw = 0.5, s = 100, color = headparal, edgecolors=edge) #wild type head amplitude

ax1.scatter(11*np.ones(15), wt[:,tail], lw = 0.5, s = 100, color = control, edgecolors=edge) #wild type head amplitude
ax1.scatter(12*np.ones(15), bwm[:,tail], lw = 0.5, s = 100, color = headparal, edgecolors=edge) #wild type head amplitude
ax1.legend(ncol = 2)

for i in range(15):
    ax1.plot([1, 2], [wt[i,head], ns[i,head]], 'k', lw = 0.3, zorder = 0) #wild type head amplitude
    ax1.plot([4, 5], [wt[i,tail], ns[i,tail]], 'k', lw = 0.3, zorder = 0) #wild type head amplitude
    ax1.plot([8, 9], [wt[i,head], bwm[i,head]], 'k', lw = 0.3, zorder = 0) #wild type head amplitude
    ax1.plot([11, 12], [wt[i,tail], bwm[i,tail]], 'k', lw = 0.3, zorder = 0) #wild type head amplitude
####################
####################

for i in range(15):
    ax11.plot([1, 2], [wt2[i,head], ns2[i,head]], 'k', lw = 0.3, zorder = 0) #wild type head amplitude
    ax11.plot([4, 5], [wt2[i,neck], ns2[i,neck]], 'k', lw = 0.3, zorder = 0) #wild type head amplitude
    ax11.plot([7, 8], [wt2[i,tail], ns2[i,tail]], 'k', lw = 0.3, zorder = 0) #wild type head amplitude
    
    ax11.plot([11, 12], [wt2[i,head], bwm2[i,head]], 'k', lw = 0.3, zorder = 0) #wild type head amplitude
    ax11.plot([14, 15], [wt2[i,neck], bwm2[i,neck]], 'k', lw = 0.3, zorder = 0) #wild type head amplitude
    ax11.plot([17, 18], [wt2[i,tail], bwm2[i,tail]], 'k', lw = 0.3, zorder = 0) #wild type head amplitude

ax11.scatter(1*np.ones(15), wt2[:,head], lw = 0.5, s = 100, color = control, edgecolors=edge, label = 'Control') #wild type head amplitude
ax11.scatter(2*np.ones(15), ns2[:,head], lw = 0.5, s = 100, color = midbodyparal, edgecolors=edge, label = 'Midbody paralyzed') #wild type head amplitude
ax11.scatter(4*np.ones(15), wt2[:,neck], lw = 0.5, s = 100, color = control, edgecolors=edge) #wild type head amplitude
ax11.scatter(5*np.ones(15), ns2[:,neck], lw = 0.5, s = 100, color = midbodyparal, edgecolors=edge) #wild type head amplitude
ax11.scatter(7*np.ones(15), wt2[:,tail], lw = 0.5, s = 100, color = control, edgecolors=edge) #wild type head amplitude
ax11.scatter(8*np.ones(15), ns2[:,tail], lw = 0.5, s = 100, color = midbodyparal, edgecolors=edge) #wild type head amplitude


ax11.scatter(11*np.ones(15), wt2[:,head], lw = 0.5, s = 100, color = control, edgecolors=edge) #wild type head amplitude
ax11.scatter(12*np.ones(15), bwm2[:,head], lw = 0.5, s = 100, color = midbodyparal, edgecolors=edge) #wild type head amplitude
ax11.scatter(14*np.ones(15), wt2[:,neck], lw = 0.5, s = 100, color = control, edgecolors=edge) #wild type head amplitude
ax11.scatter(15*np.ones(15), bwm2[:,neck], lw = 0.5, s = 100, color = midbodyparal, edgecolors=edge) #wild type head amplitude
ax11.scatter(17*np.ones(15), wt2[:,tail], lw = 0.5, s = 100, color = control, edgecolors=edge) #wild type head amplitude
ax11.scatter(18*np.ones(15), bwm2[:,tail], lw = 0.5, s = 100, color = midbodyparal, edgecolors=edge) #wild type head amplitude

################################
legend_head = ax1.legend(bbox_to_anchor=(0.4, -0.19, 0.6, .102), 
ncol=2, borderaxespad=0., fontsize = 14, handletextpad=0.0, markerscale = 1.0, borderpad = 0.0, columnspacing = 0.0, frameon=False)

legend_neck = ax11.legend(bbox_to_anchor=(0.4, -0.19, 0.6, .102), 
ncol=2, borderaxespad=0., fontsize = 14, handletextpad=0.0, markerscale = 1.0, borderpad = 0.0, columnspacing = 0.0, frameon=False)
################################
ax1.set_xlim(0, 13)
#ax1.set_xticks([1, 2, 4, 5, 8, 9, 11, 12])
ax1.set_xticks([1.5, 4.5, 8.5, 11.5])
ax1.set_xticklabels([])
ax1.set_ylabel('Amplitude', fontsize = 16)
ax1.set_ylim(0, 49.9)


ax11.set_xlim(0, 19)
#ax11.set_xticks([1, 2, 4, 5, 7, 8, 11, 12, 14, 15, 17, 18])
ax11.set_xticks([1.5, 4.5, 7.5, 11.5, 14.5, 17.5])
ax11.set_xticklabels([])
ax11.set_ylabel('Amplitude', fontsize = 16)
ax11.set_ylim(0, 49.9)
################################
xloc = [1.5, 4.5, 8.5, 11.5]
yloc = -3.2
labels = ['Head', 'Tail']
for i in range(4):
    ax1.text(xloc[i], yloc, labels[i%2], ha = 'center', fontsize = 14, rotation = 0)

xloc = [1.5, 4.5, 7.5, 11.5, 14.5, 17.5]
yloc = -3.2
labels = ['Head', 'Neck', 'Tail']
for i in range(6):
    ax11.text(xloc[i], yloc, labels[i%3], ha = 'center', fontsize = 14, rotation = 0)
##################################################################################
##############################################################################
################################
yloc = 46.7
ax1.text(1.4, yloc, 'i)', ha = 'center', fontsize = 14)
ax1.text(3.2, yloc, 'Neuron\nablation', ha = 'center', va = 'center', fontsize = 14)
ax1.text(7.1, yloc, 'ii)', ha = 'center', va = 'center', fontsize = 14)
ax1.text(10.2, yloc, 'Body wall\nmuscle ablation', ha = 'center', va = 'center', fontsize = 14)


ax11.text(2.2, yloc, 'i)', ha = 'center', fontsize = 14)
ax11.text(5.0, yloc, 'Neuron\nablation', ha = 'center', va = 'center', fontsize = 14)
ax11.text(10.2, yloc, 'ii)', ha = 'center', va = 'center', fontsize = 14)
ax11.text(15.0, yloc, 'Body wall\nmuscle ablation', ha = 'center', va = 'center', fontsize = 14)

#################################
ax1.axvspan(0, 6.3, color = '0.91', lw = 0, zorder = 0)
ax1.axvspan(6.7, 13, color = '0.91', lw = 0, zorder = 0)


ax11.axvspan(0, 9.3, color = '0.91', lw = 0, zorder = 0)
ax11.axvspan(9.7, 19, color = '0.91', lw = 0, zorder = 0)
##################################################################################

##############################################################################
############ Curvature  #######
ax2.imshow(curv.T, cmap=plt.get_cmap('seismic'), aspect='auto', vmin = -10, vmax = 10, origin = 'lower', interpolation = 'none')

############ Curvature  #######
curv = np.loadtxt('./midbody_paralysis/BWM/script/curv_%d'%23)
ax22.imshow(curv.T, cmap=plt.get_cmap('seismic'), aspect='auto', vmin = -10, vmax = 10, origin = 'lower', interpolation = 'none')
for ax in [ax2, ax22]:
    ax.set_xticks([])
    ax.set_yticks([6, 16, 41])
    ax.set_yticklabels(['Head', 'Neck', 'Tail'], fontsize = 14)
    ax.set_xticks([0, 100, 200])
    ax.set_xticklabels(['0', '5', '10'], fontsize = 16)
    ax.set_xlim(0, 210)
    ax.set_ylim(-0.5, 48.5)
    ax.set_xlabel('Time (s)', fontsize = 18)

#######################################################################################
#######################################################################################
m = axo1.imshow(curv_prof, cmap=plt.get_cmap('jet'), aspect= 'auto', origin = 'lower', vmin = 0, vmax = 1, interpolation = 'spline36')

barim = fig.colorbar(m, ax = axo1, cax = axoc, ticks=[0, 0.5, 1])
axoc.set_ylabel('Bending', labelpad = 17, rotation = -90, fontsize = 16)
axo2.errorbar(np.arange(8), np.mean(speeds_mean, axis = 0), np.std(speeds_mean, axis = 0), color = 'k', marker = 'o', capsize=4, elinewidth=1.4)
for ax in [axo1, axo2]:
    ax.set_xlim(-0.5, 7.5)
    ax.set_xticks([0, 2, 4, 6])
    ax.set_xticklabels([])
    ax.set_xticklabels(['0', '0.5', '2', '8'], fontsize = 14)
axo2.set_xlabel('Added synaptic strength', fontsize = 16)
#axo2.set_xlabel('Gap junction overexpression in B-class motorneurons', fontsize = 15)

axo1.set_ylabel('Body\ncoordenates', labelpad = -22, fontsize = 16)
#axo1.set_ylabel('Body\ncoordenates', fontsize = 15)
axo2.set_ylabel('Normalized\nSpeed', fontsize = 16)

axo1.set_yticks([4, 45])
axo1.set_yticklabels(['Head', 'Tail'], fontsize = 16)

axo2.set_ylim(-0.02, 1.05)
axo2.set_yticks([0, 0.5, 1.0])

################################################################################################################################################################
fig.subplots_adjust(left=0.042, bottom=0.06, right=0.98, top=0.97)

plt.savefig('./Consistency.png', dpi = 250)
plt.savefig('./Consistency.pdf')#, dpi = 250)


#######################################################################################
#######################################################################################
################################################################################################################################################################
#######################################################################################
################################################################################################################################################################
#######################################################################################
################################################################################################################################################################
#######################################################################################
################################################################################################################################################################
########################################################################################
##########################################################################
#np.save('myfile.npy', [selected, wt, ns, bwm, wt, ns, bwm, curv, curv_prof, speeds_mean])
### Load back in
#selected, wt, ns, bwm, wt, ns, bwm, curv, curv_prof, speeds_mean = np.load('myfile.npy')


#####################
####################
#####################
#plt.close('all')
#fig = plt.figure(figsize = [16,7])
#gm = gridspec.GridSpec(80, 300)
#ax1 = plt.subplot(gm[:50, :78])
#ax11 = plt.subplot(gm[:50, 97:175])

#ax2 = plt.subplot(gm[60:77, :78])
#ax22 = plt.subplot(gm[60:77, 97:175])
#####################
#axo1 = plt.subplot(gm[:52, 200:288])
#axoc = plt.subplot(gm[:52, 291:294])
#axo2 = plt.subplot(gm[60:77, 200:288])
##########################################################################################
########################################################################
#x1, x2, x3, y1, y2 = 0.015, 0.31, 0.61, 0.96, 0.30
#letters = [[x1, y1, 'A'], [x1, y2, 'B'], [x2, y1, 'C'], [x2, y2, 'D'], [x3, y1, 'E'], [x3, y2, 'F']]
#for l in letters:
#    plt.figtext(l[0], l[1], l[2], fontsize = 34, ha = 'center', va = 'center')
##########################################################################################
########################################################################
#selected = np.loadtxt('./figure_2_Fouad/BWM/list_selected')
#wt, ns, bwm = np.zeros((15, 3)), np.zeros((15, 3)), np.zeros((15, 3))

#i= 0
#for ind in selected:
#    wt[i] = np.loadtxt('./figure_2_Fouad/data/amplitude_fouad_wt_%d'%ind)
#    bwm[i] = np.loadtxt('./figure_2_Fouad/data/amplitude_fouad_bwm_%d'%ind)
#    ns[i] = np.loadtxt('./figure_2_Fouad/data/amplitude_fouad_ns_%d'%ind)
#    i+=1

##control, paralyzed, edge = '#1b9e77','#d95f02','#7570b3'
#control, headparal, midbodyparal, edge = '#1b9e77','#d95f02','#7570b3', 'k'

#[head, neck, tail] = [0, 1, 2]




#ax1.scatter(1*np.ones(15), wt[:,head], lw = 0.5, s = 100, color = control, edgecolors=edge, label = 'Control') #wild type head amplitude
#ax1.scatter(2*np.ones(15), ns[:,head], lw = 0.5, s = 100, color = headparal, edgecolors=edge, label = 'Head paralyzed') #wild type head amplitude
#ax1.scatter(4*np.ones(15), wt[:,tail], lw = 0.5, s = 100, color = control, edgecolors=edge) #wild type head amplitude
#ax1.scatter(5*np.ones(15), ns[:,tail], lw = 0.5, s = 100, color = headparal, edgecolors=edge) #wild type head amplitude

#ax1.scatter(8*np.ones(15), wt[:,head], lw = 0.5, s = 100, color = control, edgecolors=edge) #wild type head amplitude
#ax1.scatter(9*np.ones(15), bwm[:,head], lw = 0.5, s = 100, color = headparal, edgecolors=edge) #wild type head amplitude

#ax1.scatter(11*np.ones(15), wt[:,tail], lw = 0.5, s = 100, color = control, edgecolors=edge) #wild type head amplitude
#ax1.scatter(12*np.ones(15), bwm[:,tail], lw = 0.5, s = 100, color = headparal, edgecolors=edge) #wild type head amplitude
#ax1.legend(ncol = 2)


#for i in range(15):
#    ax1.plot([1, 2], [wt[i,head], ns[i,head]], 'k', lw = 0.3, zorder = 0) #wild type head amplitude
#    ax1.plot([4, 5], [wt[i,tail], ns[i,tail]], 'k', lw = 0.3, zorder = 0) #wild type head amplitude
#    ax1.plot([8, 9], [wt[i,head], bwm[i,head]], 'k', lw = 0.3, zorder = 0) #wild type head amplitude
#    ax1.plot([11, 12], [wt[i,tail], bwm[i,tail]], 'k', lw = 0.3, zorder = 0) #wild type head amplitude
#####################
#####################
#wt2, ns2, bwm2 = np.zeros((15, 3)), np.zeros((15, 3)), np.zeros((15, 3))

#i= 0
#for ind in selected:
#    wt2[i] = np.loadtxt('./figure_3_Fouad/data/amplitude_fouad_WT_%d'%ind)
#    bwm2[i] = np.loadtxt('./figure_3_Fouad/data/amplitude_fouad_BWM_%d'%ind)
#    ns2[i] = np.loadtxt('./figure_3_Fouad/data/amplitude_fouad_NS_%d'%ind)
#    i+=1
#[head, neck, tail] = [0, 1, 2]

#for i in range(15):
#    ax11.plot([1, 2], [wt2[i,head], ns2[i,head]], 'k', lw = 0.3, zorder = 0) #wild type head amplitude
#    ax11.plot([4, 5], [wt2[i,neck], ns2[i,neck]], 'k', lw = 0.3, zorder = 0) #wild type head amplitude
#    ax11.plot([7, 8], [wt2[i,tail], ns2[i,tail]], 'k', lw = 0.3, zorder = 0) #wild type head amplitude
#    
#    ax11.plot([11, 12], [wt2[i,head], bwm2[i,head]], 'k', lw = 0.3, zorder = 0) #wild type head amplitude
#    ax11.plot([14, 15], [wt2[i,neck], bwm2[i,neck]], 'k', lw = 0.3, zorder = 0) #wild type head amplitude
#    ax11.plot([17, 18], [wt2[i,tail], bwm2[i,tail]], 'k', lw = 0.3, zorder = 0) #wild type head amplitude

#ax11.scatter(1*np.ones(15), wt2[:,head], lw = 0.5, s = 100, color = control, edgecolors=edge, label = 'Control') #wild type head amplitude
#ax11.scatter(2*np.ones(15), ns2[:,head], lw = 0.5, s = 100, color = midbodyparal, edgecolors=edge, label = 'Midbody paralyzed') #wild type head amplitude
#ax11.scatter(4*np.ones(15), wt2[:,neck], lw = 0.5, s = 100, color = control, edgecolors=edge) #wild type head amplitude
#ax11.scatter(5*np.ones(15), ns2[:,neck], lw = 0.5, s = 100, color = midbodyparal, edgecolors=edge) #wild type head amplitude
#ax11.scatter(7*np.ones(15), wt2[:,tail], lw = 0.5, s = 100, color = control, edgecolors=edge) #wild type head amplitude
#ax11.scatter(8*np.ones(15), ns2[:,tail], lw = 0.5, s = 100, color = midbodyparal, edgecolors=edge) #wild type head amplitude


#ax11.scatter(11*np.ones(15), wt2[:,head], lw = 0.5, s = 100, color = control, edgecolors=edge) #wild type head amplitude
#ax11.scatter(12*np.ones(15), bwm2[:,head], lw = 0.5, s = 100, color = midbodyparal, edgecolors=edge) #wild type head amplitude
#ax11.scatter(14*np.ones(15), wt2[:,neck], lw = 0.5, s = 100, color = control, edgecolors=edge) #wild type head amplitude
#ax11.scatter(15*np.ones(15), bwm2[:,neck], lw = 0.5, s = 100, color = midbodyparal, edgecolors=edge) #wild type head amplitude
#ax11.scatter(17*np.ones(15), wt2[:,tail], lw = 0.5, s = 100, color = control, edgecolors=edge) #wild type head amplitude
#ax11.scatter(18*np.ones(15), bwm2[:,tail], lw = 0.5, s = 100, color = midbodyparal, edgecolors=edge) #wild type head amplitude

#################################
#legend_head = ax1.legend(bbox_to_anchor=(0.4, -0.19, 0.6, .102), 
#ncol=2, borderaxespad=0., fontsize = 14, handletextpad=0.0, markerscale = 1.0, borderpad = 0.0, columnspacing = 0.0, frameon=False)

#legend_neck = ax11.legend(bbox_to_anchor=(0.4, -0.19, 0.6, .102), 
#ncol=2, borderaxespad=0., fontsize = 14, handletextpad=0.0, markerscale = 1.0, borderpad = 0.0, columnspacing = 0.0, frameon=False)
#################################
#ax1.set_xlim(0, 13)
##ax1.set_xticks([1, 2, 4, 5, 8, 9, 11, 12])
#ax1.set_xticks([1.5, 4.5, 8.5, 11.5])
#ax1.set_xticklabels([])
#ax1.set_ylabel('Amplitude', fontsize = 16)
#ax1.set_ylim(0, 49.9)


#ax11.set_xlim(0, 19)
##ax11.set_xticks([1, 2, 4, 5, 7, 8, 11, 12, 14, 15, 17, 18])
#ax11.set_xticks([1.5, 4.5, 7.5, 11.5, 14.5, 17.5])
#ax11.set_xticklabels([])
#ax11.set_ylabel('Amplitude', fontsize = 16)
#ax11.set_ylim(0, 49.9)
#################################
#xloc = [1.5, 4.5, 8.5, 11.5]
#yloc = -3.2
#labels = ['Head', 'Tail']
#for i in range(4):
#    ax1.text(xloc[i], yloc, labels[i%2], ha = 'center', fontsize = 14, rotation = 0)

#xloc = [1.5, 4.5, 7.5, 11.5, 14.5, 17.5]
#yloc = -3.2
#labels = ['Head', 'Neck', 'Tail']
#for i in range(6):
#    ax11.text(xloc[i], yloc, labels[i%3], ha = 'center', fontsize = 14, rotation = 0)
###################################################################################
###############################################################################
#################################
#yloc = 46.7
#ax1.text(1.4, yloc, 'i)', ha = 'center', fontsize = 14)
#ax1.text(3.2, yloc, 'Neuron\nablation', ha = 'center', va = 'center', fontsize = 14)
#ax1.text(7.1, yloc, 'ii)', ha = 'center', va = 'center', fontsize = 14)
#ax1.text(10.2, yloc, 'Body wall\nmuscle ablation', ha = 'center', va = 'center', fontsize = 14)


#ax11.text(2.2, yloc, 'i)', ha = 'center', fontsize = 14)
#ax11.text(5.0, yloc, 'Neuron\nablation', ha = 'center', va = 'center', fontsize = 14)
#ax11.text(10.2, yloc, 'ii)', ha = 'center', va = 'center', fontsize = 14)
#ax11.text(15.0, yloc, 'Body wall\nmuscle ablation', ha = 'center', va = 'center', fontsize = 14)


##ax11.text(4.7, yloc, 'i)\nNeuron\nablation', ha = 'center', va = 'center', fontsize = 14)
##ax11.text(14.5, yloc, 'ii)\nBody wall muscle\nablation', ha = 'center', va = 'center', fontsize = 14)

##ax11.text(4.7, yloc, 'Neuron Ablation', ha = 'center', fontsize = 14)
##ax11.text(14.5, yloc, 'BWM ablation', ha = 'center', fontsize = 14)
##################################
#ax1.axvspan(0, 6.3, color = '0.91', lw = 0, zorder = 0)
#ax1.axvspan(6.7, 13, color = '0.91', lw = 0, zorder = 0)


#ax11.axvspan(0, 9.3, color = '0.91', lw = 0, zorder = 0)
#ax11.axvspan(9.7, 19, color = '0.91', lw = 0, zorder = 0)
###################################################################################
###############################################################################
##lstyle = "|-|, widthA=0.3, widthB=0.3"
##yloc = 46
##ax1.annotate('', xy=(0.2, yloc), xytext=(6.3, yloc),arrowprops=dict(arrowstyle=lstyle, color='k'))
##ax1.annotate('', xy=(6.7, yloc), xytext=(12.8, yloc),arrowprops=dict(arrowstyle=lstyle, color='k'))

##ax11.annotate('', xy=(0.2, yloc), xytext=(9.3, yloc),arrowprops=dict(arrowstyle=lstyle, color='k'))
##ax11.annotate('', xy=(9.7, yloc), xytext=(18.8, yloc),arrowprops=dict(arrowstyle=lstyle, color='k'))
###################################################################################
###############################################################################
############# Curvature  #######
#curv = np.loadtxt('./figure_3_Fouad/BWM/script/curv_%d'%23)
#ax22.imshow(curv.T, cmap=plt.get_cmap('seismic'), aspect='auto', vmin = -10, vmax = 10, origin = 'lower', interpolation = 'none')

############# Curvature  #######
#curv = np.loadtxt('./figure_2_Fouad/BWM/script/curv_%d'%23)
#ax2.imshow(curv.T, cmap=plt.get_cmap('seismic'), aspect='auto', vmin = -10, vmax = 10, origin = 'lower', interpolation = 'none')
#for ax in [ax2, ax22]:
#    ax.set_xticks([])
#    ax.set_yticks([6, 16, 41])
#    ax.set_yticklabels(['Head', 'Neck', 'Tail'], fontsize = 14)
#    ax.set_xticks([0, 100, 200])
#    ax.set_xticklabels(['0', '5', '10'], fontsize = 16)
#    ax.set_xlim(0, 210)
#    ax.set_ylim(-0.5, 48.5)
#    ax.set_xlabel('Time (s)', fontsize = 18)


########################################################################################
########################################################################################
########################################################################################
########################################################################################
##########################################################################

#curv_prof = np.loadtxt('./overexpression/curv_prof')
#speeds_mean = np.loadtxt('./overexpression/speed')

#m = axo1.imshow(curv_prof, cmap=plt.get_cmap('jet'), aspect= 'auto', origin = 'lower', vmin = 0, vmax = 1, interpolation = 'spline36')
##m = axo1.imshow(curv_prof, cmap=plt.get_cmap('jet'), aspect= 'auto', origin = 'lower', vmin = 0, vmax = 1, interpolation = 'None')

#barim = fig.colorbar(m, ax = axo1, cax = axoc, ticks=[0, 0.5, 1])
#axoc.set_ylabel('Bending', labelpad = 17, rotation = -90, fontsize = 16)
#axo2.errorbar(np.arange(8), np.mean(speeds_mean, axis = 0), np.std(speeds_mean, axis = 0), color = 'k', marker = 'o', capsize=4, elinewidth=1.4)
#for ax in [axo1, axo2]:
#    ax.set_xlim(-0.5, 7.5)
#    ax.set_xticks([0, 2, 4, 6])
#    ax.set_xticklabels([])
#    ax.set_xticklabels(['0', '0.5', '2', '8'], fontsize = 14)
#axo2.set_xlabel('Added synaptic strength', fontsize = 16)
##axo2.set_xlabel('Gap junction overexpression in B-class motorneurons', fontsize = 15)

#axo1.set_ylabel('Body\ncoordenates', labelpad = -22, fontsize = 16)
##axo1.set_ylabel('Body\ncoordenates', fontsize = 15)
#axo2.set_ylabel('Normalized\nSpeed', fontsize = 16)

#axo1.set_yticks([4, 45])
#axo1.set_yticklabels(['Head', 'Tail'], fontsize = 16)

#axo2.set_ylim(-0.02, 1.05)
#axo2.set_yticks([0, 0.5, 1.0])

#################################################################################################################################################################
#fig.subplots_adjust(left=0.042, bottom=0.06, right=0.98, top=0.97)

#plt.savefig('./Consistency_poster.png', dpi = 250)
##plt.savefig('./Consistency.pdf')#, dpi = 250)
##plt.savefig('../Figures/Fouad_poster.png', dpi = 250)

#np.save('myfile.npy', [selected, wt, ns, bwm, wt2, ns2, bwm2, curv, curv_prof, speeds_mean])




