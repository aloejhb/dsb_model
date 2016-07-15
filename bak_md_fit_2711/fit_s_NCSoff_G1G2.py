from model_s import md_s as model
from data_funcs import *
from fit_funcs import *
from pylab import *
import copy
import numdifftools as nd

trackPath = '../results/replot_data/'
exp = 'NCSoff_G1G2'
sfreq = 6.
startFrm = 200
dtxind = 3

track_list, data = loadData(trackPath, exp, dtxind, sfreq, startFrm)

pars = {'D':60., 'alph':0.1, 'alph2':0.1, 'beta':0, 'beta2':0.07, 'rNH':0.25, 'rHR':0.05, 'rHR2':0.05, 'conv':0, 'conv2':1.} 
pars.update({'dp':0, 'dp2':0.32, 'dp_5r':0, 'dp2_5r':0.32})
pars.update({'b':0.45})
pars.update({'ttreat':5., 'ttreat_dl':0.05, 'twash':5.7, 'twash_dl':0.05, 'ts':8.5, 'ts_dl':0.2, 'tg2':15., 'tg2_dl':0.2})

model.setAllPars(pars)
 
mdxind = 5
trange = (0, 60)

# pltCompr2Data(model, data, mdxind, trange)
pltCompr2Data(model, data, mdxind, trange, ssim = True, numSim = 50, svfig = True)

par_grid = {'conv2':linspace(0.001, 30., 30), 'rHR2':linspace(0.001, 0.2, 30)}
pnms = par_grid.keys()
lscp_md = copy.deepcopy(model)
# lscp_md.pars['b_to_rNH'] = 3.
# xx, yy, obj = landscape(lscp_md, data, mdxind, par_grid, trange, svfig = True)



# without upreg of repair
fixpars = {'b':0.45} 
fitpars = {'D':10., 'alph2':4., 'beta2':0.2, 'rNH':0.15, 'rHR':0.001, 'conv2':20.}
bd_dict = {'D':(0.01, 20.), 'alph':(0.1, 10.), 'rNH':(0.05, 1.), 'twash':(3., 60.), 'twash_dl':(1., 30.)}

# fit(model, data, mdxind, fitpars, fixpars, None, trange, 'Nelder-Mead')
# rndFit(model, data, mdxind, fitpars, fixpars, bd_dict, trange, fit_met, numRnd = 10)

