from pfuncs_ssa import *
from model_funcs import *
from data_funcs import *
from fit_funcs import *
import copy

# ion()

trackPath = '../data/'
exps = ('NCS_G1s',)
sfreq = 6.
startFrm = 200
dtxind = 3

track_lists = {}
data_list = {}
for i,exp in enumerate(exps):
    trl, dt= loadData(trackPath, exp, dtxind, sfreq, startFrm)
    track_lists[exp] = trl
    data_list[exp] = dt

ttime = 100
tvec = linspace(0, ttime, sfreq*ttime + 1)

# pltData((data_list['NCS_G1s'],))

mdname = 'g1g1'

def mk_pfuncs(p):
    if 'twash' not in p:
        p['twash'] = p['ttreat']+p['ttreat_dl']+p['gap_ttreat_twash']

    pfuncs = []
    pfuncs.append(spline_pfunc((p['b'], p['D'], p['b']), (p['ttreat'], p['twash']), (p['ttreat_dl'], p['twash_dl']), mk_hfunc(()), ttime))

    pfuncs.append(spline_pfunc((p['alph'],), (), (), mk_hfunc((0,)), ttime))

    pfuncs.append(spline_pfunc((p['beta'],), (), (), mk_hfunc((0,)), ttime))

    pfuncs.append(spline_pfunc((p['rNH'], p['rNHup'],), (p['treg'],), (p['treg_dl'],), mk_hfunc((1,)), ttime))

    pfuncs.append(spline_pfunc((p['rHR'],), (), (), mk_hfunc((2,)), ttime))

    return pfuncs

x0 = [0, 1, 0]
stoich = array(([1, 0, 0], [-1, 1, 0], [-1, 0, 1], [0, -1, 0], [0, 0, -1]))

pars = {'D':5.5, 'b':0.2, 'alph':4., 'beta':0.001, 'rNH':0.1, 'rNHup':0.1, 'rHR':0.001, 'ttreat':5., 'ttreat_dl':0.2, 'twash':7.5, 'twash_dl':1., 'treg':98., 'treg_dl':2.}
pfuncs = mk_pfuncs(pars)

xname = ['gammaH2AX', '53BP1+gammaH2AX', 'RAD51+gammaH2AX']
xcolor = ['dodgerblue', 'red', 'green']
detxcolor = ['cornflowerblue', 'rosybrown', 'lime']

obsinds = [(0,1,2)]
obsname = ['total gammaH2AX']
obscolor = ['grey']
detobscolor = ['black']

pname = ['DNA_dam', '53BP1_recrt', 'RAD51_recrt', 'NHEJ_repr', 'HR_repr']
pcolor = ['blue', 'red', 'green', 'm', 'lime']

fuzz = None
# fuzz = {'twash':('exponential', 10.)}

model = Model(mdname, x0, stoich, ttime, mk_pfuncs, pars, xname, xcolor, detxcolor, pname, pcolor, obsinds, obsname, obscolor, detobscolor, fuzz)


# model.pltSim()
# model.pltDetSim()
# model.pltMultiSim(numSim = 5)
 
dtname = 'NCS_G1s'
print 'Working on ' + dtname + 'data.\n'



mdxind = 1
trange = (0, 60)

# pltCompr2Data(model, data_list[dtname], mdxind, trange)
# pltCompr2Data(model, data_list[dtname], mdxind, trange, ssim = True, numSim = 50)

# NCS data fitting 
fixpars = {'b':1.5, 'beta':0.01, 'rHR':0.001}
fitpars = {'D':20., 'alph':0.8, 'rNH':0.4, 'rNHup':0.6, 'ttreat':12., 'ttreat_dl':0.5, 'gap_ttreat_twash':1.,'twash_dl':10., 'treg':6., 'treg_dl':2.}
bd_dict = {'D':(10., 50.), 'rNH':(0.5, 1.), 'ttreat':(10., 30.), 'ttreat_dl':(0.5, 10.), 'gap_ttreat_twash':(0.3, 10.), 'twash_dl':(1., 30.), 'treg':(1., 30.), 'treg_dl':(1., 30.)}
 
fit_met = 'Nelder-Mead'

# fit(model, data_list[dtname], mdxind, fitpars, fixpars, bd_dict, trange, fit_met)
# rndFit(model, data_list[dtname], mdxind, fitpars, fixpars, bd_dict, trange, fit_met, numRnd = 10)

par_grid = {'treg':linspace(0.3, 99, 20), 'twash':linspace(8, 50, 20)}
pnms = par_grid.keys()
# xx, yy, obj = landscape(model, data_list[dtname], mdxind, par_grid, trange, svfig = True)

# show()
