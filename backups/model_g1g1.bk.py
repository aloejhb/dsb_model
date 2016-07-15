from pfuncs_ssa import *
from model_funcs import *
from data_funcs import *
from fit_funcs import *
import copy

# ion()
 
trackPath = '../data/'
exps = ('G1G1_Doxo', 'NCS_G1s')
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

# pltData((data_list['G1G1_Doxo'], data_list['NCS_G1s']))

mdname = 'g1g1'

def mk_pfuncs(p):
    if 'twash' not in p:
        p['twash'] = p['ttreat']+p['ttreat_dl']+p['gap_ttreat_twash']

    if 'b_to_rNH' in p:
        p['b'] = p['rNH'] * p['b_to_rNH']

    pfuncs = []
    pfuncs.append(spline_pfunc((p['b'], p['D'], p['b']), (p['ttreat'], p['twash']), (p['ttreat_dl'], p['twash_dl']), mk_hfunc(()), ttime))

    pfuncs.append(spline_pfunc((p['alph'],), (), (), mk_hfunc((0,)), ttime))

    pfuncs.append(spline_pfunc((p['beta'],), (), (), mk_hfunc((0,)), ttime))

    pfuncs.append(spline_pfunc((p['rNH'], p['rNHup'],), (p['treg'],), (p['treg_dl'],), mk_hfunc((1,)), ttime))

    pfuncs.append(spline_pfunc((p['rHR'],), (), (), mk_hfunc((2,)), ttime))

    return pfuncs

x0 = [0, 1, 0]
stoich = array(([1, 0, 0], [-1, 1, 0], [-1, 0, 1], [0, -1, 0], [0, 0, -1]))

# pars = {'D':2., 'b':0.45, 'alph':4., 'beta':0.001, 'rNH':0.15, 'rNHup':0.3, 'rHR':0.001, 'ttreat':2., 'ttreat_dl':0.15, 'twash':22., 'twash_dl':4.8, 'treg':2.5, 'treg_dl':11.}
pars = {'D':1.16, 'b':0.45, 'alph':0.6, 'beta':0.001, 'rNH':0.0001, 'rNHup':0.2, 'rHR':0.001, 'ttreat':2., 'ttreat_dl':0.15, 'twash':22., 'twash_dl':7.3, 'treg':4., 'treg_dl':7.2}
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

model_list = {}
md_pars = {}
model_list['G1G1_Doxo'] = copy.deepcopy(model)

model_list['NCS_G1s'] = copy.deepcopy(model)
# md_pars['NCS_G1s']= {'D':5.9, 'b':0.45, 'alph':4., 'beta':0.001, 'rNH':0.15, 'rNHup':0.1, 'rHR':0.001, 'ttreat':5., 'ttreat_dl':0.2, 'twash':7.5, 'twash_dl':1., 'treg':98., 'treg_dl':2.}
md_pars['NCS_G1s']= {'D':30., 'b':0.45, 'alph':0.6, 'beta':0.001, 'rNH':0.15, 'rNHup':0.1, 'rHR':0.001, 'ttreat':5., 'ttreat_dl':0.2, 'twash':5.6, 'twash_dl':0.2, 'treg':98., 'treg_dl':2.}
model_list['NCS_G1s'].setAllPars(md_pars['NCS_G1s'])

model_list['basal'] = copy.deepcopy(model)
md_pars['basal']= {'D':0.45, 'b':0.45, 'alph':4., 'beta':0.001, 'rNH':0.15, 'rNHup':0.1, 'rHR':0.001, 'ttreat':5., 'ttreat_dl':0.2, 'twash':7.5, 'twash_dl':1., 'treg':98., 'treg_dl':2.}
model_list['basal'].setAllPars(md_pars['basal'])

# model.pltSim()
# model.pltDetSim()
# model.pltMultiSim(numSim = 5)

mdxind = 1
trange = (0, 60)

dtname = 'G1G1_Doxo'
# dtname = 'NCS_G1s'
print 'Working on ' + dtname + ' data.\n'
md = model_list[dtname]



# pltCompr2Data(model_list[dtname], data_list[dtname], mdxind, trange)
# pltCompr2Data(model_list[dtname], data_list[dtname], mdxind, trange, ssim = True, numSim = 50)

par_grid = {'D':linspace(10., 50., 20), 'rNH':linspace(0.01, 0.5, 20)}
pnms = par_grid.keys()
lscp_md = copy.deepcopy(model_list[dtname])
lscp_md.pars['b_to_rNH'] = 3.
# xx, yy, obj = landscape(lscp_md, data_list[dtname], mdxind, par_grid, trange, svfig = True)

fixpars = {}
fitpars = {}
bd_dict = {}

# Doxo data fitting 
fixpars['G1G1_Doxo'] = {'b':0.45, 'beta':0.001, 'ttreat':2., 'rHR':0.001, 'treg':98., 'treg_dl':1., 'rNHup': 0.3}
fitpars['G1G1_Doxo'] = {'D':2., 'alph':4., 'rNH':0.15, 'ttreat_dl':0.15, 'twash':22., 'twash_dl':4.8}
# bd_dict['G1G1_Doxo'] = {'D':(0.01, 20.), 'rNH':(0.5, 1.), 'ttreat_dl':(1., fixpars['G1G1_Doxo']['twash']-fixpars['G1G1_Doxo']['ttreat']-0.1), 'twash_dl':(1., 30.), 'treg':(1., 30.), 'treg_dl':(1., 30.)}
bd_dict['G1G1_Doxo'] = {'D':(0.01, 20.), 'rNH':(0.5, 1.), 'ttreat_dl':(0.15, 10.), 'twash':(11., 40.), 'twash_dl':(1., 30.)}

# NCS data fitting 
fixpars['NCS_G1s'] = {'b':1.5, 'beta':0.01, 'rHR':0.001}
fitpars['NCS_G1s'] = {'D':20., 'alph':0.8, 'rNH':0.4, 'rNHup':0.6, 'ttreat':12., 'ttreat_dl':0.5, 'gap_ttreat_twash':1.,'twash_dl':10., 'treg':6., 'treg_dl':2.}
bd_dict['NCS_G1s'] = {'D':(10., 50.), 'rNH':(0.5, 1.), 'ttreat':(10., 30.), 'ttreat_dl':(0.5, 10.), 'gap_ttreat_twash':(0.3, 10.), 'twash_dl':(1., 30.), 'treg':(1., 30.), 'treg_dl':(1., 30.)}
 
fit_met = 'Nelder-Mead'

# fit(model_list[dtname], data_list[dtname], mdxind, fitpars[dtname], fixpars[dtname], bd_dict[dtname], trange, fit_met)
# rndFit(model_list[dtname], data_list[dtname], mdxind, fitpars[dtname], fixpars[dtname], bd_dict[dtname], trange, fit_met, numRnd = 10)


# show()
