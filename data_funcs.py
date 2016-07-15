from pylab import *
from loadTracks import Load_Tracks
from model_funcs import plotMnStd
from model_funcs import plotPerc
from scipy.stats import scoreatpercentile as perc

class FociData:
    def __init__(self, name, tvec, all, mn, st, q1, q3):
        self.name = name
        self.tvec = tvec
        self.all = all
        self.mn = mn
        self.st = st
        self.q1 = q1
        self.q3 = q3

def getDataMatrix(dtname, track_list, dtxind, sfreq, startFrm = 1):
    # header = """Frame  XPos  YPos  Foci  Radius Geminin_Int  53BP1_Int mFociArea"""    
    # indices ->   0 ,   1,     2,   3,     4,     5,           6,       7
    lastFrm = max([tr[:,0][-1] for tr in track_list])

    tvec = linspace(0, (lastFrm - startFrm)/sfreq, lastFrm - startFrm + 1)

    all = []
    for tr in track_list:
        frame = tr[:,0]
        dtvec = tr[:,dtxind]
        dtvec = dtvec.tolist()
        # pre-tach and attach NAN as place holder for missing data
        if frame[0] < startFrm:
            dtvec = dtvec[int(startFrm - frame[0]):]
        else:
            pre_nan = frame[0] - startFrm
            dtvec[:0] = [nan] * pre_nan

        ex_nan = lastFrm - frame[-1]
        dtvec.extend([nan] * ex_nan)

        all.append(dtvec)

    all = array(all)
    mn = nanmean(all, axis = 0)
    st = nanstd(all, axis = 0)

    q1 = []
    q3 = []
    for cl in all.T:
        cl2 = [c for c in cl if not isnan(c)]
        q1.append(perc(cl2, 25))
        q3.append(perc(cl2, 75))
    q1 = array(q1)
    q3 = array(q3)

    return FociData(dtname, tvec, all, mn, st, q1, q3)

def loadData(trackPath, exp, dtxind, sfreq, startFrm):
    if trackPath[-1] != '/':
        trackPath +='/'

    exp_path = exp.replace('_', '/')
    if exp_path[-1] != '/':
        exp_path +='/'

    trl = Load_Tracks(trackPath + exp_path)
    dt = getDataMatrix(exp, trl, dtxind, sfreq, startFrm)

    return trl, dt

def pltData(data_list, labs = None, colors = None, newfig = True):
    if colors is None:
        colors = cm.rainbow(linspace(0,1,len(data_list)))

    if newfig:
        figure()
        xlabel('Time a.u.')
        ylabel('#Foci')

    if labs is None:
        labs = [dt.name for dt in data_list]
    for i,dt in enumerate(data_list):
        # plotMnStd(dt.tvec, dt.mn, dt.st, label = dt.name, color = colors[i]) 
        plotPerc(dt.tvec, dt.mn, dt.q1, dt.q3, label = labs[i], color = colors[i]) 

    legend(frameon = False, prop = {'size':12})
    

