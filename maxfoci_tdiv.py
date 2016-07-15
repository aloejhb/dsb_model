from pylab import *
from data_funcs import *
from loadTracks import *
from cc_strat import *

sfreq = 6.
startFrm = 200
treatFrm = 230

ttreat = (treatFrm - startFrm) / sfreq

ipath = '../results/replot_data/'
colors = ('royalblue', 'firebrick')
# exps = ('NCS_G1S', 'NCSoff_G1S')
# labs = ('MYCNon_G1S','MYCNoff_G1S')
exps = ('NCS_G2M', 'NCSoff_G2M')
labs = ('MYCNon_G2M','MYCNoff_G2M')

figure()
track_lists = []
for i,exp in enumerate(exps):
    exp_path = exp.replace('_', '/')
    if exp_path[-1] != '/':
        exp_path +='/'
    trl = Load_Tracks(ipath + exp_path)
    track_lists.append(trl)

    ttrnl = []
    mfocil = []
    for tr in trl:
        tag, tup, gm0, tdiv = tag_tr(tr)
        gmnn = tr[:,5]
        mfoci = max(gmnn[treatFrm-startFrm:])
        if tag == 'G1S':
            ttrn = -tup
        elif tag == 'G2M':
            ttrn = tdiv
        ttrnl.append(ttrn)
        mfocil.append(mfoci)
    plot(ttrnl, mfocil, 'o', color = colors[i], label = labs[i])
legend(frameon=False, prop={'size':12})
