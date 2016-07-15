from numpy import *
from pylab import *
from data_funcs import *
from loadTracks import *
from stratif_by_Gmnn import *
from scipy import signal
import os
import shutil
def deriv(y):
    dy = []
    for i in range(len(y)-1):
        dy.append(y[i+1] - y[i])
    dy.append(nan)
    return array(dy)

# ion()

ipath = '../data/NCS_pool/'

opath = '../results/replot_data/NCS/'
if not os.path.exists(opath):
    os.makedirs(opath)

for tg in ('G1div', 'G1G2', 'G1', 'G2div', 'G2'):
    if not os.path.exists(opath + '/' + tg):
        os.makedirs(opath + '/' + tg)
    
    

sfreq = 6.
startFrm = 200

treatFrm = 230
ttreat = (treatFrm - startFrm) / sfreq

tlist, fnames = Load_Tracks(ipath, fn = True)

tnames = [f.replace('.txt', '') for f in fnames]

def drop_regn(y, drv_th = -1000, len_th = 3, smwd = 5):

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
add_wd = 12


GMTH = 1000.
tdiv_g1 = []
tdiv_g2 = []
# for i in range(10):
firups, tcross = getTcross(tlist, GMTH)
print firups
for i in range(len(tlist)):
    tr = tlist[i]
    if tr[0,0] < startFrm:
        tr2 = tr[int(startFrm - tr[0,0]):]
    else:
        tr2 = tr
    tvec = (tr2[:,0]-startFrm)/sfreq #frames to time
    gmnn = tr2[:,5]
    rad = tr2[:,4]

    t_cut = tvec[treatFrm-startFrm:]
    gmnn_cut = gmnn[treatFrm-startFrm:]
    rad_cut = rad[treatFrm-startFrm:]

    gmnn_rg = drop_regn(gmnn_cut, drv_th = -500, len_th = 5, smwd = 10)
    
    if len(gmnn_rg):
        gmnn_wrg = [max(0, gmnn_rg[0]-add_wd), min(len(t_cut)-1, gmnn_rg[-1]+add_wd)]
     
        rad_cut2 = rad_cut[gmnn_wrg[0]:gmnn_wrg[1]]
        t_cut2 = t_cut[gmnn_wrg[0]:gmnn_wrg[1]]
        rad_rg = drop_regn(rad_cut2, drv_th = -1, len_th = 2, smwd = 10)

        print rad_cut2
        print t_cut2
        print rad_rg

        if len(rad_rg):
            gmdrop = 2 # cell division
            tdiv = mean(t_cut2[rad_rg[0]:rad_rg[1]])
        else:
            gmdrop = 1 # G2 arrest
    else:
        gmdrop = 0 # no Gmnn drop

    if mean(gmnn_cut[:6]) < GMTH:
        if firups[i] is None:
            tag = 'G1'
        elif firups[i] > ttreat and firups[i] < ttreat + 2:
            tag = 'G1G2'
        else:
            tag = 'G1'
    else:
        if gmdrop == 2:
            tag = 'G2div'
            tdiv_g2.append(tdiv)
        else:
            tag = 'G2'



    fig, axs = Plot_FociGemRad(tr, sfreq, startFrm)
    axs[0].axvline(ttreat, color = 'grey', linestyle = '--')
    axs[1].axvline(ttreat, color = 'grey', linestyle = '--')
    axs[2].axvline(ttreat, color = 'grey', linestyle = '--')

    axs[0].set_ylim(0,35)
    axs[2].set_xlim(0,60)

    fnm = tag + '/' + tnames[i]
    savefig(opath + fnm + '.pdf')
    shutil.copyfile(ipath + tnames[i] + '.txt', opath  + fnm + '.txt')
    close()

savetxt(opath + 'tdiv_g1.txt', array(tdiv_g1))
savetxt(opath + 'tdiv_g2.txt', array(tdiv_g2))

'''
            axs[2].plot(tvec_cut[rad_rg], rad_cut[rad_rg], 'o', color = 'r')

    axs[0].axvline(ttreat, color = 'grey', linestyle = '--')
    axs[1].axvline(ttreat, color = 'grey', linestyle = '--')
    axs[2].axvline(ttreat, color = 'grey', linestyle = '--')
                
    

        axs[1].axvspan(tvec[gmnn_wrg[0]], tvec[gmnn_wrg[1]], facecolor='g', alpha=0.3)
        axs[1].plot(tgm_cut[gmnn_rg], gmnn_cut[gmnn_rg], 'o', color = 'g')
        axs[2].axvspan(tvec[gmnn_wrg[0]], tvec[gmnn_wrg[1]], facecolor='g', alpha=0.3)
'''



'''


# for i in range(5):
for i in range(len(track_list)):
    fig, axs = Plot_FociGemRad(track_list[i], sfreq, startFrm)

    axs[0].axvline(ttreat, color = 'grey', linestyle = '--')
    axs[1].axvline(ttreat, color = 'grey', linestyle = '--')
    axs[2].set_xlim(0,60)

    savefig(figPath + exp_path + tnames[i] + '.pdf')
    close()
'''
