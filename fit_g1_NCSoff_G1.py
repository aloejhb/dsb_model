from model_g1 import md_g1 as model
from data_funcs import *
from fit_funcs import *
import copy
import numdifftools as nd

trackPath = '../results/replot_data/'
exp = 'NCSoff_G1'
sfreq = 6.
startFrm = 200
dtxind = 3

track_list, data = loadData(trackPath, exp, dtxind, sfreq, startFrm)

# pars = {'D':30., 'alph':0.066,'rNH':0.24}
pars = {'D':30., 'alph':0.065,'rNH':0.235}
pars.update({'b':0.18})
pars.update({'ttreat':5., 'ttreat_dl':0.05, 'twash':6., 'twash_dl':0.05})

model.setAllPars(pars)

mdxind = 1
trange = (0, 30)
trange = (0, 55)

# pltCompr2Data(model, data, mdxind, trange)
# pltCompr2Data(model, data, mdxind, trange, ssim = True, numSim = 50)
figname = 'MYCNoff_G1'
data.name = figname
pltCompr2Data(model, data, mdxind, trange, ssim = True, numSim = 50, svfig = True)
xlim(trange)
fig, ax1, ax2 = model.pltMultiSim(numSim=50)
ax1.set_ylim(0,45)
fig.suptitle(figname + ' model simulation', fontsize=18)
fig.savefig('../results/manual_fit/' + figname + '_model_sim.pdf')



par_grid = [('twash',linspace(1., 40., 20)), ('treg',linspace(1., 70., 20))]
lscp_md = copy.deepcopy(model)
# xx, yy, obj = landscape(lscp_md, data, mdxind, par_grid, trange, svfig = True)


fitpars = {}
fixpars = {} 
fitpars_nm = ['alph', 'rNH']
for nm in fitpars_nm:
    fitpars[nm] = pars[nm]

for nm in pars.keys():
    if nm not in fitpars_nm:
        fixpars[nm] = pars[nm]

bd_dict = {}


# fit(model, data, mdxind, fitpars, fixpars, None, trange, 'Nelder-Mead')
# rndFit(model, data, mdxind, fitpars, fixpars, bd_dict, trange, fit_met, numRnd = 10)

