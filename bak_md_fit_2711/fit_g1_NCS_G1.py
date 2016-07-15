from model_g1 import md_g1 as model
from data_funcs import *
from fit_funcs import *
import copy
import numdifftools as nd

trackPath = '../results/replot_data/'
exp = 'NCS_G1'
sfreq = 6.
startFrm = 200
dtxind = 3

track_list, data = loadData(trackPath, exp, dtxind, sfreq, startFrm)

pars1 = {'D':30., 'alph':0.5,'rNH':0.13, 'rNHup':0.15}
pars2 = {'b':0.45}
pars3 = {'ttreat':5., 'ttreat_dl':0.05, 'twash':5.7, 'twash_dl':0.05, 'treg':98., 'treg_dl': 1.}
pars = {}
pars.update(pars1)
pars.update(pars2)
pars.update(pars3)

model.setAllPars(pars)

mdxind = 1
trange = (0, 60)

# pltCompr2Data(model, data, mdxind, trange)
# pltCompr2Data(model, data, mdxind, trange, ssim = True, numSim = 50)

par_grid = [('twash',linspace(1., 40., 20)), ('treg',linspace(1., 70., 20))]
lscp_md = copy.deepcopy(model)
# lscp_md.pars['b_to_rNH'] = 3.
# xx, yy, obj = landscape(lscp_md, data, mdxind, par_grid, trange, svfig = True)


# with upreg
fixpars = {'b':0.45, 'beta':0.001, 'rHR':0.001, 'ttreat':2., 'ttreat_dl':0.15}
fitpars = {'D':2., 'alph':4., 'rNH':0.15, 'rNHup': 0.3, 'twash':23., 'twash_dl':7.3, 'treg':5.7, 'treg_dl':2.}
bd_dict = {'D':(0.01, 20.), 'alph':(0.1, 10.), 'rNH':(0.05, 1.), 'rNHup':(0.05, 1.), 'twash':(10., 40.), 'twash_dl':(1., 30.), 'treg':(3., 20.), 'treg':(1., 30.)}

# fit(model, data, mdxind, fitpars, fixpars, None, trange, 'Nelder-Mead')
# rndFit(model, data, mdxind, fitpars, fixpars, bd_dict, trange, fit_met, numRnd = 10)


'''
pnames = ['D', 'alph', 'rNH']
p0 = [pars[pn] for pn in pnames]
fixpars2 = dict([p for p in pars.items() if p[0] not in pnames])

def objfunc2(p):
        tvec = data.tvec
        data_mn = data.mn
        return objFunc(p, model, tvec, data_mn, mdxind, pnames, fixpars2)


Jacob_objfunc = nd.Jacobian(objfunc2)
jfunc = Jacob_objfunc(p0)

Hess_objfunc = nd.Hessian(objfunc2)
hsfunc = Hess_objfunc(p0)

ion()
figure()
imshow(hsfunc, interpolation='none', vmin=hsfunc.min(), vmax=hsfunc.max(), aspect='auto', cmap=cm.coolwarm)
colorbar()
xticks(range(len(pnames)), pnames, rotation = 45)
yticks(range(len(pnames)), pnames, rotation = 45)

savefig('../results/others/Hess_g1_vs_NCS_G1.png')
'''
