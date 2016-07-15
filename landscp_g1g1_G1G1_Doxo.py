from pfuncs_ssa import *
from model_funcs import *
from data_funcs import *
from fit_funcs import *
import copy

trackPath = '../data/'
exps = ('G1G1_Doxo','NCS_G1s')
sfreqs = (4., 6.)
startFrms = (1, 200)
dtxind = 3

track_lists = {}
data_list = {}
for i,exp in enumerate(exps):
    trl, dt= loadData(trackPath, exp, dtxind, sfreqs[i], startFrms[i])
    track_lists[exp] = trl
    data_list[exp] = dt

ttime = 100
tvec = linspace(0, ttime, sfreq*ttime + 1)

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

pars = {'D':4., 'b':1.5, 'alph':0.8, 'beta':0.001, 'rNH':0.4, 'rNHup':0.7, 'rHR':0.001, 'ttreat':2., 'ttreat_dl':0.5, 'twash':22., 'twash_dl':10., 'treg':5.5, 'treg_dl':5.}
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

model = Model(mdname, x0, stoich, ttime, mk_pfuncs, pars, xname, xcolor, detxcolor, pname, pcolor, obsinds, obsname, obscolor, detobscolor, fuzz)


mdxind = 1
trange = (0, 60)

dtname = 'G1G1_Doxo'
print 'Working on ' + dtname + ' data.\n'

par_grid = {'treg':linspace(0.3, 70, 20), 'twash':linspace(8, 35, 20)}
pnms = par_grid.keys()
# xx, yy, obj = landscape(model, data_list[dtname], mdxind, par_grid, trange, svfig = True)

# fit other par for each minima
fixpars = {'b':1.5, 'beta':0.001, 'ttreat':2., 'rHR':0.001
# fixpars.update({'twash':23., 'treg':8.})
fixpars.update({'twash':10., 'treg':40.})
fitpars = {'D':4., 'alph':0.8, 'rNH':0.4, 'rNHup':0.7, 'ttreat_dl':1., 'twash_dl':10., 'treg_dl':5.}
bd_dict = {'D':(0.01, 20.), 'rNH':(0.5, 1.), 'ttreat_dl':(1., fixpars['twash']-fixpars['ttreat']-0.1), 'twash_dl':(1., 30.), 'treg_dl':(1., 30.)}


# fit from minima
'''
fixpars = {'b':1.5, 'beta':0.001, 'ttreat':2., 'rHR':0.001, 'twash':22., 'alph':0.8, 'D':3., 'rNH':0.4, 'rNHup':0.6, 'ttreat_dl':5., 'twash_dl':10., 'treg':6., 'treg_dl':2.}}
fitpars = {'twash':23., 'treg':8.}
# fitpars = {'twash':10., 'treg':40.}
bd_dict = {'twash':(0.3, 35.), 'treg':(0.3, 70.)}
'''

fit_met = 'Nelder-Mead'
# fit(model, data_list['G1G1_Doxo'], mdxind, fitpars, fixpars, bd_dict, trange, fit_met)
# rndFit(model, data_list['G1G1_Doxo'], mdxind, fitpars, fixpars, bd_dict, trange, fit_met, numRnd = 10)
