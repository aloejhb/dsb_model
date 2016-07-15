from pylab import *
from scipy import signal
from scipy.interpolate import UnivariateSpline as spline


GM1 = 700.
GM2 = 4000.
CURPHASE_TRANGE1 = 0
CURPHASE_TRANGE2 = 4 

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

def drop_regn(tvec, y, drv_th, len_th = 3, smwd = 5):

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
    if len(regn):
        return tvec[regn[0]],tvec[regn[-1]]
    else:
        return None

# detect wether there is g1s transition in the given vector
# return transition time or None
def g1s_trans(tr):
    trans = up_time(tr[0], tr[1], GM1)
    return trans

def sg2_trans(tr):
    trans = up_time(tr[0], tr[1], GM2)
    return trans

def g2g1_trans(tr, add_rg=2., gmdth=-500., gmlth=5, gms=10, rddth=-1, rdlth=2, rds=7):
    tvec = tr[0]
    gm = tr[1]
    rd = tr[2]
    gm_rg = drop_regn(tvec, gm, drv_th = gmdth, len_th = gmlth, smwd = gms)
    rd_rg = None
    if gm_rg:
        gm_wrg = (gm_rg[0] - add_rg, gm_rg[1] + add_rg)
        rgind = (tvec > gm_wrg[0]) & (tvec < gm_wrg[1])
        rd_rg = drop_regn(tvec[rgind], rd[rgind], drv_th = rddth, len_th = rdlth, smwd = rds)
    if rd_rg:
        return mean(rd_rg)
    else:
        return None

def cur_phase(tr, tp, trange1 = CURPHASE_TRANGE1, trange2 = CURPHASE_TRANGE2):
    tr_cur = tr[:,(tr[0] > tp - trange1) & (tr[0] < tp + trange2)]
    gmnn_mn = mean(tr_cur[1])
    if g2g1_trans(tr_cur):
        cph = 'M'
    else:
        if gmnn_mn < GM1:
            cph = 'G1'
        elif gmnn_mn < GM2:
            cph = 'S'
        else:
            cph = 'G2'
    return cph, gmnn_mn
        
    
def tag_tr(trk, ttreat, startFrm = 200, treatFrm = 230, sfreq = 6):

    if trk[0,0] < startFrm:
        trk = trk[int(startFrm - trk[0,0]):]
    tvec = (trk[:,0]-startFrm)/sfreq #frames to time
    gmnn = trk[:,5]
    rad = trk[:,4]
    tr = array((tvec, gmnn, rad))

    tag, gmnn_mn = cur_phase(tr, ttreat)
        
    tr2 = tr[:,tr[0]>ttreat]
    if tag == 'M':
        ttrans = g1s_trans(tr2)
        if ttrans:
            #tag += 'G1S'
            tag = 'G1S'
        else:
            #tag +='G1'
            tag = 'G1'
    elif tag == 'G1':
        ttrans = g1s_trans(tr2)
        if ttrans:
            tag += 'S'
    elif tag == 'S':
        ttrans = sg2_trans(tr2)
        #if ttrans:
            #tag += 'G2'
        tag = 'SG2'
    else:
        ttrans = g2g1_trans(tr2)
        if ttrans:
            tag += 'M'

    return tag, gmnn_mn, ttrans
