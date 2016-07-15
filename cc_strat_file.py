from pylab import *
from cc_strat import *
from data_funcs import *
from loadTracks import *
from scipy import signal
from scipy.interpolate import UnivariateSpline as spline
import os
import shutil


sfreq = 6.
startFrm = 200
treatFrm = 230

ttreat = (treatFrm - startFrm) / sfreq

ipath = '../data/NCSpool/'

opath = '../results/replot_data/NCS/'
if not os.path.exists(opath):
    os.makedirs(opath)
# tag_order = ['MG1', 'MG1S', 'G1', 'G1S', 'S', 'SG2', 'G2', 'G2M']
tag_order = ['G1', 'G1S', 'SG2', 'G2', 'G2M']
# for tg in ('G1', 'G1S', 'S', 'SG2', 'G2', 'G2M'):
for tg in tag_order:
    if not os.path.exists(opath + '/' + tg):
        os.makedirs(opath + '/' + tg)

tlist, fnames = Load_Tracks(ipath, fn = True)
tnames = [f.replace('.txt', '') for f in fnames]

for i,tr in enumerate(tlist):
    tag, gmnn_mn, ttrans = tag_tr(tr, ttreat)

    fig, axs = Plot_FociGemRad(tr, sfreq, startFrm)
    for ax in axs:
        ax.axvline(ttreat, color = 'grey', linestyle = '--')

    axs[0].set_ylim(0,35)
    axs[2].set_xlim(0,60)

    fnm = tag + '/' + tnames[i]
    savefig(opath + fnm + '.pdf')
    shutil.copyfile(ipath + tnames[i] + '.txt', opath  + fnm + '.txt')
    close()



