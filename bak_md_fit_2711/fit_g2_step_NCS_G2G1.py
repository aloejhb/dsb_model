from model_g2 import md_g2_step as model
from data_funcs import *
from fit_funcs import *
from pylab import *
import copy
import numdifftools as nd

trackPath = '../data/'
exp = 'NCS_G2G1s'
sfreq = 6.
startFrm = 200
dtxind = 3

track_list, data = loadData(trackPath, exp, dtxind, sfreq, startFrm)

pars = {'D':50., 'alph2':0.6, 'beta2':0.2, 'rNH':0.15, 'rHR':0.08, 'rHR2':0.08, 'conv2':2.} 
pars.update({'b':0.45})
pars.update({'ttreat':5., 'twash':5.6})

model.setAllPars(pars)
 
mdxind = 5
trange = (0, 60)

# pltCompr2Data(model, data, mdxind, trange)
# pltCompr2Data(model, data, mdxind, trange, ssim = True, numSim = 50)

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


'''
pnames = ['D', 'alph2', 'beta2', 'rNH', 'rHR', 'rHR2', 'conv2']
fixpars2 = {'b':0.45, 'ttreat':5., 'twash':5.6}
tvec = data.tvec
data_mn = data.mn

def objfunc2(p):
        return objFunc(p, model, tvec, data_mn, mdxind, pnames, fixpars2)

p0 = [50., 0.6, 0.2, 0.15, 0.08, 0.08, 2.]

Jacob_objfunc = nd.Jacobian(objfunc2)
jfunc = Jacob_objfunc(p0)

Hess_objfunc = nd.Hessian(objfunc2)
hsfunc = Hess_objfunc(p0)

ion()
figure()
imshow(hsfunc, interpolation='none', vmin=-2e5, vmax=2e5, aspect='auto', cmap=cm.coolwarm)
colorbar()
xticks(range(len(pnames)), pnames, rotation = 45)
yticks(range(len(pnames)), pnames, rotation = 45)

# savefig('../results/others/Hess_g2_NCS_G2G1.pdf')
'''
