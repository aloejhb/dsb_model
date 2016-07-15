from model_g1 import md_g1 as model
from data_funcs import *
from fit_funcs import *
import copy

trackPath = '../data/'
exp = 'NCSoff_G1G1s'
sfreq = 6.
startFrm = 200
dtxind = 3

track_list, data = loadData(trackPath, exp, dtxind, sfreq, startFrm)

pars = {'D':45., 'b':0.45, 'alph':0.08,'rNH':0.25, 'rNHup':0.15}
pars.update( {'beta':0.0001, 'rHR':0.0001} )
pars.update( {'ttreat':6., 'ttreat_dl':0.05, 'twash':6.6, 'twash_dl':0.05, 'treg':98., 'treg_dl':1.} )

model.setAllPars(pars)

mdxind = 1
trange = (0, 60)

# pltCompr2Data(model, data, mdxind, trange)
# pltCompr2Data(model, data, mdxind, trange, ssim = True, numSim = 50)

par_grid = {'twash':linspace(1., 40., 20), 'treg':linspace(1., 70., 20)}
pnms = par_grid.keys()
lscp_md = copy.deepcopy(model)
# lscp_md.pars['b_to_rNH'] = 3.
# xx, yy, obj = landscape(lscp_md, data, mdxind, par_grid, trange, svfig = True)


# with upreg
fixpars = {'b':0.45, 'beta':0.001, 'rHR':0.001, 'ttreat':2., 'ttreat_dl':0.15}
fitpars = {'D':2., 'alph':4., 'rNH':0.15, 'rNHup': 0.3, 'twash':23., 'twash_dl':7.3, 'treg':5.7, 'treg_dl':2.}
bd_dict = {'D':(0.01, 20.), 'alph':(0.1, 10.), 'rNH':(0.05, 1.), 'rNHup':(0.05, 1.), 'twash':(10., 40.), 'twash_dl':(1., 30.), 'treg':(3., 20.), 'treg':(1., 30.)}

# fit(model, data, mdxind, fitpars, fixpars, None, trange, 'Nelder-Mead')
# rndFit(model, data, mdxind, fitpars, fixpars, bd_dict, trange, fit_met, numRnd = 10)
