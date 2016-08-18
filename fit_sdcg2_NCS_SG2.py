from model_sdcg2 import md_sdcg2 as model
from data_funcs import *
from fit_funcs import *
from pylab import *
import copy
import numdifftools as nd

trackPath = '../results/replot_data/'
exp = 'NCS_SG2'
sfreq = 6.
startFrm = 200
dtxind = 3

track_list, data = loadData(trackPath, exp, dtxind, sfreq, startFrm)
# pars = {'D':23.89, 'alph':0.8, 'alph_s':0.42, 'rNH':0.052, 'beta':0.02, 'rHR':0.02, 'conv':0.024, 'rHR2':0.02}
# fitted parameters
pars = {'D':23.71, 'alph':0.8, 'alph_s':0.0745, 'rNH':0.043, 'beta':0.064, 'rHR':0.02, 'conv':0.00001, 'rHR2':0.02}
pars.update({'b':0.18})
pars.update({'ttreat':5., 'ttreat_dl':0.05, 'twash':6., 'twash_dl':0.05})
pars.update({'ts':6.5, 'ts_dl':1., 'tg2':12., 'tg2_dl':2.})

model.setAllPars(pars)
 
mdxind = 5
# trange = (0, 55)
trange = (0, 40)

# pltCompr2Data(model, data, mdxind, trange)
# pltCompr2Data(model, data, mdxind, trange, ssim = True, numSim = 50)
'''
figname = 'MYCNon_SG2'
data.name = figname
pltCompr2Data(model, data, mdxind, trange, ssim = True, numSim = 50, svfig = True)
xlim(trange)
fig, ax1, ax2 = model.pltMultiSim(numSim=50)
ax1.set_ylim(0,30)
fig.suptitle(figname + ' model simulation', fontsize=18)
fig.savefig('../results/manual_fit/' + figname + '_model_sim.png')

par_grid = [('conv',linspace(0.001, 0.2, 20)), ('rHR2',linspace(0.001, 0.1, 20))]
lscp_md = copy.deepcopy(model)
# xx, yy, obj = landscape(lscp_md, data, mdxind, par_grid, trange, svfig = True)
'''

fitpars = {}
fixpars = {} 
fitpars_nm = ['D', 'alph_s', 'beta', 'rNH', 'conv']
for nm in fitpars_nm:
    fitpars[nm] = pars[nm]

for nm in pars.keys():
    if nm not in fitpars_nm:
        fixpars[nm] = pars[nm]

bd_dict = {}


# fit(model, data, mdxind, fitpars, fixpars, None, trange, 'Nelder-Mead')
# rndFit(model, data, mdxind, fitpars, fixpars, bd_dict, trange, fit_met, numRnd = 10)

