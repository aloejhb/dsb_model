from model_g2 import md_g2 as model
from data_funcs import *
from fit_funcs import *
from pylab import *
import copy
import numdifftools as nd

trackPath = '../data/'
exp = 'NCSoff_G2G1s'
sfreq = 6.
startFrm = 200
dtxind = 3

track_list, data = loadData(trackPath, exp, dtxind, sfreq, startFrm)

pars = {'D':50., 'alph2':0.25, 'beta2':0.2, 'rNH':0.15, 'rHR':0.08, 'rHR2':0.08, 'conv2':2.} 
pars.update({'b':0.45})
pars.update( {'ttreat':6., 'ttreat_dl':0.05, 'twash':6.6, 'twash_dl':0.05} )

model.setAllPars(pars)
 
mdxind = 5
trange = (0, 60)

# pltCompr2Data(model, data, mdxind, trange)
# pltCompr2Data(model, data, mdxind, trange, ssim = True, numSim = 50)

par_grid = {'beta2':linspace(0.001, 1., 20), 'rHR2':linspace(0.001, 0.2, 20)}
pnms = par_grid.keys()
lscp_md = copy.deepcopy(model)
# lscp_md.pars['b_to_rNH'] = 3.
# xx, yy, obj = landscape(lscp_md, data, mdxind, par_grid, trange, svfig = True)



# without upreg of repair
fixpars = {'b':0.45} 
fitpars = {'D':10., 'alph2':4., 'beta2':0.2, 'rNH':0.15, 'rHR':0.001, 'conv2':20., 'ttreat':5., 'twash':5.6}
# bd_dict = {'D':(0.01, 20.), 'alph':(0.1, 10.), 'rNH':(0.05, 1.), 'twash':(3., 60.), 'twash_dl':(1., 30.)}

# fit(model, data, mdxind, fitpars, fixpars, None, trange, 'Nelder-Mead')
# rndFit(model, data, mdxind, fitpars, fixpars, bd_dict, trange, fit_met, numRnd = 10)

