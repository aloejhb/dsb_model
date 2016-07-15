from numpy import *
from pylab import *
from scipy.interpolate import UnivariateSpline as spline

###
sfreq = 6.
def crossThTime(spl, tstart, tend, th, eps = 0.1):
    t = tstart
    tcr = []
    direct = []
    while t < tend:
        if (spl(t+eps) - th) * (spl(t) - th) < 0:
            tcr.append(t)       
            if spl(t) - th < 0:
                direct.append('+')
            else:
                direct.append('-')
        t = t + eps
    return tcr, direct

def getTcross(track_list, gmnnTh, spl_s = 5e6, pltfig = False):
    firUps = []
    tcross = []
    for tr in track_list:
        frame = tr[:,0]
        gm = tr[:,5]
###
        tvec = (frame-200)/sfreq 
        spl = spline(tvec, gm, s = spl_s)
        tcr, direct = crossThTime(spl, tvec[0], tvec[-1], gmnnTh)
        if len(direct) and direct[0] == '+':
            firUp = tcr[0]
        else:
            firUp = None
        firUps.append(firUp)
        tcross.append((tcr, direct))
        
        if pltfig:
            figure()
            plot(tvec, gm, 'o')
            etvec = linspace(tvec[0], tvec[-1], 1000)
            plot(etvec, spl(etvec))
            axhline(y = gm_th, color = 'k', alpha = 0.7)
    return firUps, tcross 

def stratifByGmnn(track_list, gmnnTh = 1000., timeTh = 6., plthist = False):
    firUps, tcross = getTcross(track_list, gmnnTh)
    fup =  [x for x in firUps if x is not None]
    if plthist:
        binw = 1
        figure()
        hist(fup, range(0, 50 + binw, binw))
    trl = [x for i,x in enumerate(trlist_sort) if firUps[i] is not None]
    ear_trl = [trl_s[i] for i, x in enumerate(fup) if x < div]
    late_trl = [trl_s[i] for i, x in enumerate(fup) if x >= div]

    return ear_trl, late_trl


