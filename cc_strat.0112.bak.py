from pylab import *
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


def up_time(tvec, y, th, sm = 5e6, eps = 0.1):
    spl = spline(tvec, y, s = sm)
    tup = None
    t = tvec[0]
    while t < tvec[-1]:
        if spl(t) < th and (spl(t+eps) - th) * (spl(t) - th) < 0:
            tup = t
            break 
        t = t + eps
    return tup

def deriv(y):
    dy = []
    for i in range(len(y)-1):
        dy.append(y[i+1] - y[i])
    dy.append(nan)
    return array(dy)

def drop_regn(y, drv_th, len_th = 3, smwd = 5):

    dy = deriv(y)
    win = signal.hann(smwd)
    smdy = signal.convolve(dy, win, mode = 'same')/sum(win)

    drop = []
    for i,d in enumerate(smdy):
        if d < drv_th:
            drop.append(i)

    dif_drop = deriv(drop)
    regn = []
    for i,dd in enumerate(dif_drop):
        if dd == 1:
            regn.append(drop[i])
        else:
            if len(regn) >= len_th - 1:
                regn.append(drop[i])
                break
            else:
                regn = []
    return regn

def div_time(tvec, gm, rd, add_wd=12, gmdth=-500., gmlth=5, gms=10, rddth=-1, rdlth=2, rds=10):

    gm_rg = drop_regn(gm, drv_th = gmdth, len_th = gmlth, smwd = gms)
    if len(gm_rg):
        gm_wrg = [max(0, gm_rg[0]-add_wd), min(len(tvec)-1, gm_rg[-1]+add_wd)]

        rd2 = rd[gm_wrg[0]:gm_wrg[1]]
        tvec2 = tvec[gm_wrg[0]:gm_wrg[1]]
        rd_rg = drop_regn(rd2, drv_th = rddth, len_th = rdlth, smwd = rds)

        if len(rd_rg):
            drop_typ = 2 # cell division
            tdiv = mean(tvec2[rd_rg[0]:rd_rg[1]])
        else:
            drop_typ = 1 # G2 arrest
            tdiv = None
    else:
        drop_typ = 0 # no Gmnn drop
        tdiv = None

    return tdiv, drop_typ

GMTH1 = 800.
GMTH2 = 4500.

def tag_tr(tr):
    if tr[0,0] < startFrm:
        tr = tr[int(startFrm - tr[0,0]):]
    tvec = (tr[:,0]-startFrm)/sfreq #frames to time
    gmnn = tr[:,5]
    rad = tr[:,4]

    tvec2 = tvec[treatFrm-startFrm:]
    gmnn2 = gmnn[treatFrm-startFrm:]
    rad2 = rad[treatFrm-startFrm:]
    
    gm_treat = mean(gmnn2[:int(sfreq)*3])
    tup = None
    tdiv = None
    if gm_treat < GMTH1:
        tup = up_time(tvec2, gmnn2, GMTH1)
        if tup:
            tag = 'G1S'
        else:
            tag = 'G1'
    elif gm_treat < GMTH2:
        tup = up_time(tvec2, gmnn2, GMTH2)
        if tup:
            tag = 'SG2'
        else:
            tag = 'S'
    else:
        tdiv, drop_type = div_time(tvec2, gmnn2, rad2)
        if tdiv:
            tag = 'G2M'
        else:
            tag = 'G2'
    return tag, gm_treat, tup, tdiv 



'''
ipath = '../data/NCSoffpool/'

opath = '../results/replot_data/NCSoff/'
if not os.path.exists(opath):
    os.makedirs(opath)

for tg in ('G1', 'G1S', 'S', 'SG2', 'G2', 'G2M'):
    if not os.path.exists(opath + '/' + tg):
        os.makedirs(opath + '/' + tg)
    
    

tlist, fnames = Load_Tracks(ipath, fn = True)
tnames = [f.replace('.txt', '') for f in fnames]



tdiv_g2 = []
for i,tr in enumerate(tlist):
    tag, gm_treat, tup, tdiv = tag_tr(tr)

    if tag == 'G2div':
            tdiv_g2.append(tdiv)

    fig, axs = Plot_FociGemRad(tr, sfreq, startFrm)
    for ax in axs:
        ax.axvline(ttreat, color = 'grey', linestyle = '--')

    axs[0].set_ylim(0,35)
    axs[2].set_xlim(0,60)

    fnm = tag + '/' + tnames[i]
    savefig(opath + fnm + '.pdf')
    shutil.copyfile(ipath + tnames[i] + '.txt', opath  + fnm + '.txt')
    close()

savetxt(opath + 'tdiv_g2.txt', array(tdiv_g2))

'''


