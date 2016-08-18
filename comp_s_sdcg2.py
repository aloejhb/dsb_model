import numpy as np
import matplotlib.pyplot as plt
from fit_funcs import pltCompr2Data

from fit_s_NCS_SG2 import model as md_s
from fit_s_NCS_SG2 import data
from fit_sdcg2_NCS_SG2 import model as md_sdc

ttime = 100
tvec = np.linspace(0, ttime, 4*ttime+1)
trange = (0, ttime)

s_res = md_s.detSim(tvec)
sdc_res = md_sdc.detSim(tvec)

# plt.plot(tvec, s_res[:,0])
# pltCompr2Data(model, data, mdxind, trange)

s_res = np.c_[s_res, s_res[:,2] + s_res[:,3]]
sdc_res = np.c_[sdc_res, sdc_res[:,2] + sdc_res[:,3]]

spname = md_s.xname + ['total RAD51']

# plt.ion()
plt.rc('axes', titlesize=11, labelsize=8)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('legend', fontsize=8)

fig, axes = plt.subplots(2, 4, sharey='all', figsize=(12,6.5))

axes[-1, -1].axis('off')
for i in range(np.shape(s_res)[1]):
    ax = axes.flatten()[i]
    plt.sca(ax)
    ax.set_title(spname[i])
    plt.plot(tvec, sdc_res[:,i], color='cyan', label='Direct-competition model')
    plt.plot(tvec, s_res[:,i], color='red', label='53BP1-exclusion model', alpha=0.7)
    plt.legend()
    plt.ylim(0, 30)
    plt.xlabel('Time (h)')

axes[0,0].set_ylabel('# Foci')
axes[1,0].set_ylabel('# Foci')
plt.draw()
plt.tight_layout()
plt.savefig('../figures/' + 'comp_s_sdcg2.pdf')
