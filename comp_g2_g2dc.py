import numpy as np
import matplotlib.pyplot as plt
from fit_funcs import pltCompr2Data

from fit_g2_NCS_G2M import model as md_g2
from fit_g2_NCS_G2M import data
from fit_g2dc_NCS_G2M import model as md_g2dc

ttime = 100
tvec = np.linspace(0, ttime, 4*ttime+1)
trange = (0, ttime)

s_res = md_g2.detSim(tvec)
sdc_res = md_g2dc.detSim(tvec)

# plt.plot(tvec, s_res[:,0])
# pltCompr2Data(model, data, mdxind, trange)

spname = md_g2.xname[:2] + ['RAD51+$\gamma$H2AX'] + md_g2.xname[2:]+ ['Total RAD51']

plt.rc('axes', titlesize=11, labelsize=8)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('legend', fontsize=8)

fig, axes = plt.subplots(2, 4, sharey='all', figsize=(12,6.5))

axes[-1, -1].axis('off')
for i in range(len(spname)):
    ax = axes.flatten()[i]
    plt.sca(ax)
    ax.set_title(spname[i])
    if i < 2:
        ind_dc = i
        ind = i
    elif i == 2 :
        ind_dc = 2
        ind = None
    elif i == 3:
        ind_dc = None
        ind = 2
    elif i == 4:
        ind_dc = 3
        ind = 3
    elif i == 5:
        ind_dc = 1
        ind = 4
    else:
        ind_dc = 2
        ind = 2
    

    if ind_dc is not None:
        plt.plot(tvec, sdc_res[:,ind_dc], color='cyan', label='Direct-competition model')

    if ind is not None:
        plt.plot(tvec, s_res[:,ind], color='red', label='53BP1-depend-HR model', alpha=0.7)
        
    plt.legend()
    plt.ylim(0, 30)
    plt.xlabel('Time (h)')

axes[0,0].set_ylabel('# Foci')
axes[1,0].set_ylabel('# Foci')
plt.draw()
plt.tight_layout()
plt.savefig('../figures/' + 'comp_g2_g2dc.pdf')
