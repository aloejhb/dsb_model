from model_g2 import md_g2 as model
from data_funcs import *
from fit_funcs import *
from pylab import *
import copy
import numdifftools as nd

trackPath = '../results/replot_data/'
exp = 'NCSoff_G2M'
sfreq = 6.
startFrm = 200
dtxind = 3

track_list, data = loadData(trackPath, exp, dtxind, sfreq, startFrm)
# pars = {'D':30., 'alph':0.12, 'rNH':0.24, 'conv':0.0001, 'rHR2':0.02} 
pars = {'D':30., 'alph':0.13, 'rNH':0.24, 'conv':0., 'rHR2':0.02} 
pars1 = {'b':0.18}
pars2 = {'ttreat':5., 'ttreat_dl':0.05, 'twash':6., 'twash_dl':0.05} 
pars.update(pars1)
pars.update(pars2)

model.setAllPars(pars)

mdxind = 4
trange = (0, 55)

# pltCompr2Data(model, data, mdxind, trange)
# pltCompr2Data(model, data, mdxind, trange, ssim = True, numSim = 50)
figname = 'MYCNoff_G2M'
data.name = figname
pltCompr2Data(model, data, mdxind, trange, ssim = True, numSim = 50, svfig = True)
xlim(trange)
fig, ax1, ax2 = model.pltMultiSim(numSim=50)
ax1.set_ylim(0,45)
fig.suptitle(figname + ' model simulation', fontsize=18)
fig.savefig('../results/manual_fit/' + figname + '_model_sim.png')

par_grid = [('conv',linspace(0.001, 0.2, 20)), ('rHR2',linspace(0.001, 0.1, 20))]
lscp_md = copy.deepcopy(model)
# xx, yy, obj = landscape(lscp_md, data, mdxind, par_grid, trange, svfig = True)



fitpars = {}
fixpars = {} 
fitpars_nm = ['alph', 'conv']
for nm in fitpars_nm:
    fitpars[nm] = pars[nm]

for nm in pars.keys():
    if nm not in fitpars_nm:
        fixpars[nm] = pars[nm]

bd_dict = {}

# fit(model, data, mdxind, fitpars, fixpars, None, trange, 'Nelder-Mead')
# rndFit(model, data, mdxind, fitpars, fixpars, bd_dict, trange, fit_met, numRnd = 10)

